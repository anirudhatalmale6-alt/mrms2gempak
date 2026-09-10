#!/usr/bin/env python3
"""
mrms_s3.py -- anonymous listing and fetching of MRMS GRIB2 from AWS Open Data.

The MRMS archive lives in the public bucket `noaa-mrms-pds` (NOAA Open Data
Dissemination).  It is readable without an AWS account, without credentials and
without the aws CLI -- plain anonymous HTTPS is enough, which is why this module
uses nothing but the Python standard library.

Layout of the bucket (verified 2026-09-10):

    s3://noaa-mrms-pds/CONUS/<PRODUCT>/<YYYYMMDD>/MRMS_<PRODUCT>_<YYYYMMDD>-<HHMMSS>.grib2.gz

e.g. the 24-hour multi-sensor Pass-2 QPE valid 2026-09-09 12Z:

    CONUS/MultiSensor_QPE_24H_Pass2_00.00/20260909/
        MRMS_MultiSensor_QPE_24H_Pass2_00.00_20260909-120000.grib2.gz

Note the product name appears TWICE: once as the directory, once inside the
file name.  The trailing "_00.00" is the MRMS level string (surface), it is part
of the name and must be kept.

Cadence differs per product.  Measured on 2026-09-09 in this bucket: every QPE
accumulation product (01H through 72H, Pass 1 and Pass 2) is written hourly on
the hour -- 24 files per day, ~4.6 MB gzipped for 24H, ~0.6 MB for 01H.  The
instantaneous fields such as PrecipRate are 2-minute (717 files that day).
Rather than hard-coding any of that, this module LISTS the day's prefix and then
picks the file closest to the time you asked for, so a change of cadence
upstream cannot break it.
"""

from __future__ import annotations

import datetime as dt
import gzip
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

BUCKET = "noaa-mrms-pds"
HOST = f"https://{BUCKET}.s3.amazonaws.com"
S3_NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"

# Products this tool knows how to name.  "<dur>" is substituted.
#   MultiSensor_QPE_<dur>_Pass<n>_00.00   radar + gauge + model, Pass 1 / Pass 2
#   RadarOnly_QPE_<dur>_00.00             radar only
# Pass 2 is the later, gauge-corrected pass -- the one you want for verification.
DURATIONS = ("01H", "03H", "06H", "12H", "24H", "48H", "72H")

USER_AGENT = "mrms2gem/1.0 (+stdlib urllib)"


def product_name(duration: str, source: str = "MultiSensor", qpe_pass: int = 2) -> str:
    """Build the MRMS product string, e.g. MultiSensor_QPE_24H_Pass2_00.00."""
    duration = duration.upper()
    if duration not in DURATIONS:
        raise ValueError(f"duration {duration!r} not one of {DURATIONS}")
    if source == "MultiSensor":
        if qpe_pass not in (1, 2):
            raise ValueError("qpe_pass must be 1 or 2")
        return f"MultiSensor_QPE_{duration}_Pass{qpe_pass}_00.00"
    if source == "RadarOnly":
        # RadarOnly has no Pass distinction.
        return f"RadarOnly_QPE_{duration}_00.00"
    raise ValueError(f"unknown source {source!r} (MultiSensor or RadarOnly)")


def _http_get(url: str, retries: int = 4, timeout: int = 60) -> bytes:
    """GET with a few retries.  S3 occasionally 503s under load; that is normal."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            # A 404 is a real answer, not a transient failure -- do not retry it.
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 404:
                raise
            last = exc
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def list_day(product: str, day: dt.date, region: str = "CONUS") -> list[tuple[dt.datetime, str, int]]:
    """
    List every MRMS file for one product on one UTC day.

    Returns [(valid_time_utc, s3_key, size_bytes), ...] sorted by time.
    Uses the anonymous ListObjectsV2 REST call and follows continuation tokens,
    so a 2-minute product (720 keys/day) lists completely.
    """
    prefix = f"{region}/{product}/{day:%Y%m%d}/"
    out: list[tuple[dt.datetime, str, int]] = []
    token = None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        url = f"{HOST}/?{urllib.parse.urlencode(params)}"
        root = ET.fromstring(_http_get(url))
        for c in root.findall(f"{S3_NS}Contents"):
            key = c.findtext(f"{S3_NS}Key") or ""
            size = int(c.findtext(f"{S3_NS}Size") or 0)
            stamp = _time_from_key(key)
            if stamp is not None:
                out.append((stamp, key, size))
        if (root.findtext(f"{S3_NS}IsTruncated") or "false").lower() != "true":
            break
        token = root.findtext(f"{S3_NS}NextContinuationToken")
        if not token:
            break
    out.sort()
    return out


def _time_from_key(key: str) -> dt.datetime | None:
    """Pull the valid time out of ..._YYYYMMDD-HHMMSS.grib2.gz."""
    base = key.rsplit("/", 1)[-1]
    for suffix in (".grib2.gz", ".grib2"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    else:
        return None
    stamp = base.rsplit("_", 1)[-1]
    try:
        return dt.datetime.strptime(stamp, "%Y%m%d-%H%M%S").replace(tzinfo=dt.timezone.utc)
    except ValueError:
        return None


def find_nearest(product: str, want: dt.datetime, tolerance_min: int = 30,
                 region: str = "CONUS") -> tuple[dt.datetime, str, int] | None:
    """
    Find the file whose valid time is closest to `want`, within `tolerance_min`.

    Looks at the requested day and, if `want` is near midnight, the neighbouring
    day too -- otherwise a 00:00Z request could miss a 23:58Z file.
    """
    want = _as_utc(want)
    days = {want.date()}
    if want.hour == 0:
        days.add((want - dt.timedelta(days=1)).date())
    if want.hour == 23:
        days.add((want + dt.timedelta(days=1)).date())

    candidates: list[tuple[dt.datetime, str, int]] = []
    for day in sorted(days):
        try:
            candidates.extend(list_day(product, day, region=region))
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                raise
    if not candidates:
        return None
    best = min(candidates, key=lambda row: abs((row[0] - want).total_seconds()))
    if abs((best[0] - want).total_seconds()) > tolerance_min * 60:
        return None
    return best


def download(key: str, dest_dir: str, decompress: bool = True,
             overwrite: bool = False) -> str:
    """
    Download one S3 key into dest_dir.  Returns the path to the usable GRIB2.

    The archive stores .grib2.gz; with decompress=True the .gz is expanded and
    removed, leaving a plain .grib2 that wgrib2 and nagrib2 can both read.
    Writes to a .part file first and renames, so an interrupted run never leaves
    a truncated GRIB2 behind that a later run would happily treat as valid.
    """
    os.makedirs(dest_dir, exist_ok=True)
    base = key.rsplit("/", 1)[-1]
    gz_path = os.path.join(dest_dir, base)
    final = gz_path[:-3] if (decompress and base.endswith(".gz")) else gz_path

    if os.path.exists(final) and not overwrite and os.path.getsize(final) > 0:
        return final

    part = gz_path + ".part"
    data = _http_get(f"{HOST}/{urllib.parse.quote(key)}")
    with open(part, "wb") as fh:
        fh.write(data)

    if decompress and base.endswith(".gz"):
        with gzip.open(part, "rb") as src, open(final + ".part", "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.replace(final + ".part", final)
        os.unlink(part)
    else:
        os.replace(part, final)

    # Cheap integrity check: every GRIB2 message starts "GRIB" and ends "7777".
    with open(final, "rb") as fh:
        head = fh.read(4)
        fh.seek(-4, os.SEEK_END)
        tail = fh.read(4)
    if head != b"GRIB" or tail != b"7777":
        raise RuntimeError(f"{final} does not look like a complete GRIB2 file")
    return final


def _as_utc(when: dt.datetime) -> dt.datetime:
    return when.replace(tzinfo=dt.timezone.utc) if when.tzinfo is None else when.astimezone(dt.timezone.utc)


if __name__ == "__main__":
    # Quick self-test / browse helper:
    #   python3 mrms_s3.py 24H 2 2026-09-09
    dur = sys.argv[1] if len(sys.argv) > 1 else "24H"
    pss = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    day = dt.date.fromisoformat(sys.argv[3]) if len(sys.argv) > 3 else dt.date.today()
    prod = product_name(dur, qpe_pass=pss)
    rows = list_day(prod, day)
    print(f"{prod} {day}: {len(rows)} file(s)")
    for stamp, key, size in rows:
        print(f"  {stamp:%Y-%m-%d %H:%M:%SZ}  {size/1e6:6.2f} MB  {key}")

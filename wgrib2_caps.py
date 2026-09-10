#!/usr/bin/env python3
"""
wgrib2_caps.py -- find out what THIS wgrib2 build can do, and how it spells it.

Two things vary between wgrib2 builds and versions, and both of them matter here:

1. Interpolation (-new_grid) is only compiled in when the build had a Fortran
   compiler and USE_IPOLATES != 0.  A build without it still decodes and writes
   GRIB2 perfectly well, it just cannot regrid -- so we check before relying on
   it rather than producing a confusing error half way through a batch.

2. The spelling of "turn these sentinel values into missing data" has changed
   across versions (-rpn with a mask operator, -undefine_val, ...).  MRMS uses
   -3 for "no radar coverage" and -999 for "missing"; if those survive into an
   interpolation they are averaged into neighbouring grid boxes and you get
   negative rainfall.  So instead of trusting one spelling, we TRY the candidates
   on a real file and keep the first one that demonstrably works, where
   "demonstrably works" means: wgrib2 exits 0, the number of defined points goes
   DOWN (so something was actually masked) and the minimum value is no longer
   negative (so the sentinels are gone).  The answer is cached in
   ~/.cache/mrms2gem/caps.json so the probe runs once, not once per file.

That last check is the important one.  A mask command that silently does nothing
looks exactly like success if you only test the exit status.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys

CACHE = os.path.join(os.path.expanduser("~"), ".cache", "mrms2gem", "caps.json")

# Candidate ways of flagging MRMS sentinels as undefined.  {lo} is substituted
# with the threshold: every value below it becomes undefined.  Ordered most
# modern / most explicit first.
MASK_CANDIDATES: list[tuple[str, list[str]]] = [
    # Verified working on wgrib2 v3.1.3: "-undefine_val X" takes val or low:high.
    ("undefine_val_range", ["-undefine_val", "-1e30:{lo}"]),
    # Fallbacks for builds that predate it.
    ("rpn_mask", ["-rpn", "sto_1:{lo}:lt:mask:rcl_1:swap:merge"]),
    ("rpn_mask_simple", ["-rpn", "{lo}:lt:mask"]),
]


def wgrib2_path(explicit: str | None = None) -> str:
    """Locate wgrib2: --wgrib2 flag, then $WGRIB2, then PATH."""
    for cand in (explicit, os.environ.get("WGRIB2"), "wgrib2"):
        if not cand:
            continue
        found = shutil.which(cand) if os.path.sep not in cand else (cand if os.access(cand, os.X_OK) else None)
        if found:
            return found
    raise RuntimeError(
        "wgrib2 not found.  Install it (package `wgrib2` on Debian/Ubuntu/EPEL, "
        "or build from https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/wgrib2.tgz) "
        "and either put it on PATH or set WGRIB2=/path/to/wgrib2."
    )


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise RuntimeError(
            "command failed: " + " ".join(cmd) + "\n" + (proc.stderr or proc.stdout).strip()
        )
    return proc


def version(wgrib2: str) -> str:
    out = run([wgrib2, "--version"], check=False).stdout.strip()
    return out.splitlines()[0] if out else "unknown"


def has_option(wgrib2: str, option: str) -> bool:
    """True if this build lists `option` in its own -help output."""
    out = run([wgrib2, "-help", option.lstrip("-")], check=False)
    text = (out.stdout + out.stderr)
    return option in text and "not found" not in text.lower()


def has_interpolation(wgrib2: str) -> bool:
    """
    -new_grid needs the ipolates library, which needs a Fortran compiler at
    build time.  `wgrib2 -config` prints either "interpolation package is not
    installed" or the list of vector pairs it will interpolate.
    """
    text = run([wgrib2, "-config"], check=False).stdout
    for line in text.splitlines():
        if "interpolation" in line.lower():
            return "not installed" not in line.lower()
    return has_option(wgrib2, "-new_grid")


def field_stats(wgrib2: str, path: str, record: int = 1) -> dict:
    """
    Return {'npts','ndata','undef','min','max','mean'} for one GRIB2 record.

    -stats already reports ndata and undef alongside min/max/mean, which is what
    makes it possible to show that a mask changed something.  (There is no
    -ndata option; it is part of the -stats output.)
    """
    text = run([wgrib2, "-d", str(record), "-npts", "-stats", path]).stdout.strip()
    out: dict[str, float] = {}
    for key in ("npts", "ndata", "undef", "min", "max", "mean"):
        m = re.search(rf"\b{key}=([-\d.eE+]+)", text)
        if m:
            out[key] = float(m.group(1))
    if "npts" not in out or "min" not in out:
        raise RuntimeError(f"could not parse wgrib2 stats output:\n{text}")
    return out


def probe_mask(wgrib2: str, sample: str, threshold: float = -0.001,
               tmpdir: str | None = None) -> list[str]:
    """
    Work out how to mask MRMS sentinels with this wgrib2, proving it worked.

    `sample` must be a real MRMS GRIB2 file containing at least some -3 or -999
    values -- i.e. any CONUS file, since the grid always extends over ocean and
    off-radar terrain.  Returns the wgrib2 argument list to use.
    """
    import tempfile

    before = field_stats(wgrib2, sample)
    if before["min"] >= threshold:
        raise RuntimeError(
            f"{sample} has min={before['min']} -- no negative sentinels in it, so it "
            "cannot be used to prove the mask works.  Pass a plain CONUS MRMS file."
        )

    tmpdir = tmpdir or tempfile.mkdtemp(prefix="mrms2gem-probe-")
    os.makedirs(tmpdir, exist_ok=True)
    tried: list[str] = []
    for name, template in MASK_CANDIDATES:
        args = [a.format(lo=threshold) for a in template]
        out = os.path.join(tmpdir, f"probe_{name}.grib2")
        if os.path.exists(out):
            os.unlink(out)
        proc = run([wgrib2, sample, *args, "-grib_out", out], check=False)
        if proc.returncode != 0 or not os.path.exists(out) or os.path.getsize(out) == 0:
            tried.append(f"{name}: wgrib2 rejected it ({(proc.stderr or proc.stdout).strip().splitlines()[-1] if (proc.stderr or proc.stdout).strip() else 'no output'})")
            continue
        after = field_stats(wgrib2, out)
        # NB: wgrib2's "ndata" is the total number of data points and does not
        # change when points are masked -- "undef" is the one that moves.
        masked = after.get("undef", 0) - before.get("undef", 0)
        if masked <= 0:
            tried.append(f"{name}: ran but masked nothing (undef still {after.get('undef')})")
            continue
        if after["min"] < threshold:
            tried.append(f"{name}: masked {masked} points but min is still {after['min']}")
            continue
        return args
    raise RuntimeError(
        "none of the known ways to mask MRMS sentinel values worked with this "
        f"wgrib2 ({version(wgrib2)}).  Tried:\n  " + "\n  ".join(tried)
    )


def capabilities(wgrib2: str | None = None, sample: str | None = None,
                 refresh: bool = False) -> dict:
    """Load cached capabilities, probing if needed.  Cache key is the version string."""
    exe = wgrib2_path(wgrib2)
    ver = version(exe)
    cached: dict = {}
    if os.path.exists(CACHE) and not refresh:
        try:
            with open(CACHE) as fh:
                cached = json.load(fh)
        except (OSError, json.JSONDecodeError):
            cached = {}
    if cached.get("version") == ver and cached.get("exe") == exe and cached.get("mask_args"):
        return cached

    caps = {
        "exe": exe,
        "version": ver,
        "interpolation": has_interpolation(exe),
        "mask_args": None,
    }
    if sample:
        caps["mask_args"] = probe_mask(exe, sample)
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        with open(CACHE, "w") as fh:
            json.dump(caps, fh, indent=2)
    return caps


if __name__ == "__main__":
    exe = wgrib2_path(sys.argv[1] if len(sys.argv) > 1 else None)
    print("wgrib2        :", exe)
    print("version       :", version(exe))
    print("interpolation :", "yes" if has_interpolation(exe) else "NO (-new_grid unavailable)")
    if len(sys.argv) > 2:
        print("mask args     :", " ".join(probe_mask(exe, sys.argv[2])))

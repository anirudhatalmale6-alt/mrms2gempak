#!/usr/bin/env python3
"""
mrms2gem.py -- MRMS QPE (GRIB2, AWS) -> GEMPAK .gem, on a grid GEMPAK can use.

    ./mrms2gem.py --duration 24H --start 2026-09-09T00 --end 2026-09-09T23 \
                  --grid-template /data/nbm/blend.t12z.core.f001.co.grib2 \
                  --gemfile /data/gempak/mrms_qpe24_%Y%m%d.gem

What it does, per valid time:

  1. lists the product's day prefix in the public S3 bucket noaa-mrms-pds and
     picks the file nearest the requested time (no credentials needed);
  2. downloads and gunzips it, skipping anything already on disk;
  3. flags MRMS's sentinel values as missing   (-3 = no radar coverage,
     -999 = missing; left alone they get averaged into neighbouring boxes by the
     interpolation and come out as negative rainfall);
  4. interpolates onto your target grid with wgrib2 -new_grid, budget
     interpolation by default, which is the mode meant for precipitation
     accumulations -- it conserves the areal total instead of sampling a point;
  5. optionally re-tags the record as APCP so GEMPAK's standard GRIB2 tables
     recognise it (see --retag below);
  6. runs nagrib2 to produce/append to a .gem file.

WHY INTERPOLATE AND NOT JUST CLIP THE SAME BOX
MRMS CONUS is a 0.01 degree lat/lon grid, 7000 x 3500 = 24,500,000 points.
That is (a) far past GEMPAK's compiled-in maximum grid size and (b) a different
projection and spacing from the NBM's 2.5 km Lambert conformal grid.  Cutting
the same lat/lon box out of both with wgrib2 -small_grib leaves you with two
grids over the same ground that still have different navigation, and GEMPAK
cannot do arithmetic across different navigation -- nagrib2 will not even write
the second grid into the same .gem file.  Interpolating onto the NBM grid fixes
both problems at once: the point count drops to 3,744,965 and the navigation
block is identical, so GDDIAG/GDPLOT maths between MRMS and NBM just works.

DEPENDENCIES
  * python3 >= 3.9, standard library only (no boto3, no numpy, no aws CLI)
  * wgrib2, built WITH the interpolation (ipolates) library -- check with
    `wgrib2 -config | grep -i interpolation`.  A build that says "interpolation
    package is not installed" cannot do step 4; this script checks up front and
    tells you instead of failing half way through a batch.
  * GEMPAK 7.x with nagrib2 on PATH (source Gemenviron / Gemenviron.profile).
    Only step 6 needs GEMPAK -- use --no-gempak to stop after the GRIB2 stage.

See README.md for the full walk-through.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
import textwrap
import threading

import grid_spec
import mrms_s3
import wgrib2_caps

# MRMS sentinels: -3 means "outside radar coverage", -999 means "missing".
# Everything below this threshold is turned into GRIB2 missing data.
SENTINEL_THRESHOLD = -0.001

# The MRMS CONUS domain, read off the files themselves: 7000 x 3500 at 0.01 deg,
# lon 230.005 to 299.995 (i.e. 129.995W to 60.005W), lat 20.005 to 54.995.
MRMS_CONUS_DOMAIN = (230.005, 299.995, 20.005, 54.995)
MRMS_CONUS_POINTS = 7000 * 3500

# Point count above which a grid is worth warning about.  This is the size of the
# full NBM CONUS grid, which this client's GEMPAK demonstrably holds, so it is a
# floor on his real limit rather than a guess at it -- the limit itself is LLMXGD
# in GEMPRM.PRM and has to be read off the install.
GEMPAK_GRID_WARN = 2345 * 1597

DURATION_HOURS = {"01H": 1, "03H": 3, "06H": 6, "12H": 12, "24H": 24, "48H": 48, "72H": 72}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def log(msg: str) -> None:
    print(f"[mrms2gem] {msg}", flush=True)


def run(cmd: list[str], what: str) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{what} failed (exit {proc.returncode})\n"
            f"  command: {' '.join(cmd)}\n"
            f"  output : {(proc.stderr or proc.stdout).strip()}"
        )
    return proc


def find_exe(name: str, env_var: str | None = None, hint: str = "") -> str:
    for cand in (os.environ.get(env_var) if env_var else None, name):
        if not cand:
            continue
        found = shutil.which(cand) if os.path.sep not in cand else (cand if os.access(cand, os.X_OK) else None)
        if found:
            return found
    raise RuntimeError(f"{name} not found on PATH. {hint}")


def parse_duration(text: str) -> str:
    """
    Accept an accumulation length written any reasonable way.

    "6", "6h", "6H", "06", "06H", "6hr", "6 hour" all mean the same thing, and
    there is no reason to make someone remember that the bucket spells it 06H.
    """
    cleaned = text.strip().lower().replace("hour", "").replace("hr", "").replace("h", "").strip()
    try:
        hours = int(cleaned)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"cannot read {text!r} as an accumulation length; "
            f"use one of {', '.join(sorted(DURATION_HOURS))}"
        ) from None
    key = f"{hours:02d}H"
    if key not in DURATION_HOURS:
        raise argparse.ArgumentTypeError(
            f"MRMS has no {hours}-hour QPE product; available: "
            f"{', '.join(sorted(DURATION_HOURS))}"
        )
    return key


def normalise_clip(text: str) -> str:
    """
    Accept a clip box written any sensible way and return it as
    lonW:lonE:latS:latN with longitudes in 0-360, matching how MRMS stores them.

    Western longitudes are negative in ordinary use (-105) but 255 in the MRMS
    files, and the two have to be compared against each other, so normalise once
    here rather than in three places.
    """
    parts = text.split(":")
    if len(parts) != 4:
        raise ValueError(f"--clip needs lonW:lonE:latS:latN, got {text!r}")
    try:
        lon_w, lon_e, lat_s, lat_n = (float(v) for v in parts)
    except ValueError:
        raise ValueError(f"--clip values must be numbers, got {text!r}") from None
    lon_w = lon_w + 360.0 if lon_w < 0 else lon_w
    lon_e = lon_e + 360.0 if lon_e < 0 else lon_e
    if lon_e <= lon_w or lat_n <= lat_s:
        raise ValueError(
            f"--clip {text!r} is not a box: expected lonW:lonE:latS:latN with "
            f"west < east and south < north (got W={lon_w} E={lon_e} S={lat_s} N={lat_n})"
        )
    dw, de, ds, dn = MRMS_CONUS_DOMAIN
    if lon_e <= dw or lon_w >= de or lat_n <= ds or lat_s >= dn:
        raise ValueError(
            f"--clip {text!r} lies outside the MRMS CONUS domain "
            f"({dw}-{de}E, {ds}-{dn}N, i.e. 129.995W-60.005W)"
        )
    return f"{lon_w:.4f}:{lon_e:.4f}:{lat_s:.4f}:{lat_n:.4f}"


def glue_negative_values(argv: list[str]) -> list[str]:
    """
    Let "--clip -105:-90:35:45" work.

    argparse treats any token starting with "-" as an option, so a clip box over
    the western US -- which every clip box over the US is -- fails with "expected
    one argument" unless written --clip=-105:... This rewrites the spaced form
    into the joined form before argparse sees it.
    """
    out: list[str] = []
    skip = False
    for i, tok in enumerate(argv):
        if skip:
            skip = False
            continue
        if tok == "--clip" and i + 1 < len(argv) and argv[i + 1].startswith("-") \
                and ":" in argv[i + 1]:
            out.append(f"--clip={argv[i + 1]}")
            skip = True
        else:
            out.append(tok)
    return out


def parse_time(text: str) -> dt.datetime:
    """Accept 2026-09-09T12, 2026-09-09T12:00, 2026090912, 20260909-1200, ..."""
    text = text.strip()
    for fmt in ("%Y-%m-%dT%H:%M", "%Y-%m-%dT%H", "%Y-%m-%d %H:%M", "%Y-%m-%d %H",
                "%Y%m%d-%H%M%S", "%Y%m%d-%H%M", "%Y%m%d%H%M", "%Y%m%d%H", "%Y-%m-%d", "%Y%m%d"):
        try:
            return dt.datetime.strptime(text, fmt).replace(tzinfo=dt.timezone.utc)
        except ValueError:
            continue
    raise argparse.ArgumentTypeError(f"cannot parse time {text!r}")


def wgrib2_has_interpolation(wgrib2: str) -> bool:
    out = subprocess.run([wgrib2, "-config"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "interpolation" in line.lower():
            # wgrib2 prints either "... is not installed" or lists the vectors
            return "not installed" not in line.lower()
    return True   # older builds do not report it; let wgrib2 complain if absent


def field_stats(wgrib2: str, path: str) -> dict:
    """{'ndata','undef','min','max','mean'} for record 1, straight from -stats."""
    out = run([wgrib2, "-d", "1", "-npts", "-stats", path], "wgrib2 -stats").stdout
    stats: dict[str, float] = {}
    for token in out.strip().split(":"):
        if "=" in token:
            key, _, val = token.partition("=")
            try:
                stats[key.strip()] = float(val)
            except ValueError:
                pass
    missing = {"ndata", "undef", "min"} - stats.keys()
    if missing:
        raise RuntimeError(f"could not read {sorted(missing)} from wgrib2 output:\n{out}")
    return stats


# --------------------------------------------------------------------------- #
# the GRIB2 stage
# --------------------------------------------------------------------------- #
#: wgrib2 options this script would like to use.  Not every build has all of
#: them -- older wgrib2 predates several -- so the capabilities are probed once
#: and the command is assembled from what is actually there.
WANTED_OPTIONS = ("-set", "-set_lev", "-set_date", "-set_ftime2", "-set_ftime",
                  "-set_ftime_mode", "-set_grib_type", "-undefine_val",
                  "-small_grib", "-new_grid")


def probe_capabilities(wgrib2: str) -> dict[str, bool]:
    """Ask the wgrib2 binary which of the options we want it actually has."""
    caps = {opt: wgrib2_caps.supports(wgrib2, opt) for opt in WANTED_OPTIONS}
    missing = [opt for opt, ok in caps.items() if not ok]
    log(f"wgrib2: {wgrib2_caps.version(wgrib2).split(',')[0]}")
    if missing:
        log("wgrib2 does not have: " + " ".join(missing)
            + " -- working around where possible")
    return caps


def resolve_retag(mode: str, caps: dict[str, bool]) -> str:
    """
    Downgrade the re-tagging mode to what this wgrib2 can actually do, loudly.

    Silently producing a differently-tagged record than asked for would send us
    chasing the wrong thing in GEMPAK later, so every downgrade is logged with
    the reason and what it means.
    """
    if mode == "none":
        return mode
    if not caps["-set"]:
        log("WARNING: this wgrib2 has no -set, so the record cannot be re-tagged "
            "as APCP at all.  Falling back to --retag none: the MRMS identifiers "
            "are kept and GEMPAK needs the parameter table "
            "(see make_mrms_table.py in the README).")
        return "none"
    if mode == "apcp-acc" and not (caps["-set_ftime2"] or caps["-set_ftime"]):
        log("WARNING: this wgrib2 has neither -set_ftime2 nor -set_ftime, so the "
            "record cannot be turned into an N-hour accumulation.  Falling back to "
            "--retag apcp: still APCP at the surface, but tagged as an analysis at "
            "its valid time rather than an accumulation, so GEMPAK will not name it "
            "P06M/P12M/P24M.")
        return "apcp"
    if mode == "apcp-acc" and not caps["-set_date"]:
        # Without the reference-time shift an "0-N hour acc" record claims to be
        # valid N hours in the future.  That is a wrong date on every grid, which
        # is worse than not having the accumulation tag at all, so refuse it
        # rather than work around it.
        log("WARNING: this wgrib2 has no -set_date, so the reference time cannot be "
            "shifted back.  An accumulation record without that shift would be dated "
            f"N hours too late on every file.  Falling back to --retag apcp, which "
            "keeps the correct valid time.")
        return "apcp"
    return mode


def retag_args(mode: str, hours: int, caps: dict[str, bool] | None = None) -> list[str]:
    """
    GEMPAK's shipped GRIB2 tables have no entry for MRMS: the records carry
    discipline 209, centre 161, local table 1, category 6, parameter 41, and
    nagrib2 has nothing to call that.  Two ways out, and this picks between them.

      apcp      re-tag as the standard WMO APCP (discipline 0, category 1,
                parameter 8) at the surface, leaving it an analysis at its own
                valid time.  No GEMPAK table edits needed at all.
      apcp-acc  as above, and additionally rebuild the product definition as an
                N-hour accumulation (PDT 8) with the reference time moved back N
                hours, so the record reads "0-N hour acc fcst" ending at the MRMS
                valid time.  This is what makes GEMPAK name the grid P06M/P12M/
                P24M, matching how your NBM precip grids are already named.
      none      leave the MRMS identifiers alone.  Needs the parameter table in
                tables/ passed to nagrib2 via G2TBLS (see README).

    MRMS accumulations are valid AT the file's timestamp and cover the preceding
    N hours, which is why apcp-acc shifts the reference time back by N hours --
    without that shift the record would claim to be valid N hours in the future.
    """
    caps = caps or {opt: True for opt in WANTED_OPTIONS}
    if mode == "none":
        return []
    args = [
        "-set", "discipline", "0",
        "-set", "center", "7",
        "-set", "master_table", "2",
        "-set", "local_table", "1",
        "-set", "table_4.1", "1",
        "-set", "table_4.2", "8",
    ]
    if caps["-set_lev"]:
        args += ["-set_lev", "surface"]
    if mode == "apcp-acc":
        # resolve_retag() has already guaranteed -set_date and one of the ftime
        # options exist by the time this mode survives; -set_ftime2 is the modern
        # spelling and -set_ftime the older one, taking the same string.
        args += ["-set", "table_1.2", "1"]     # reference time = start of forecast
        args += ["-set_date", f"-{hours}hr"]
        args += ["-set_ftime2" if caps["-set_ftime2"] else "-set_ftime",
                 f"0-{hours} hour acc fcst"]
    return args


_MASK_ARGS_CACHE: dict[str, list[str]] = {}


def mask_args(wgrib2: str, caps: dict[str, bool], probe_sample: str) -> list[str]:
    """
    The wgrib2 arguments that turn MRMS sentinels into missing data.

    `-undefine_val low:high` is the clean way and exists in current wgrib2.  An
    older build without it falls back to -rpn, and the fallback has to
    demonstrably work on a real file before it is accepted -- see
    wgrib2_caps.probe_mask.  If nothing works we stop, because averaging or
    interpolating unmasked -3s produces negative rainfall.

    `probe_sample` must be the FULL CONUS download, not a clipped subset: the
    proof that a mask fired is that the undefined count went up, and a clipped
    box over well-covered land legitimately contains no sentinels to mask.  The
    probe result is cached, since each attempt rewrites a 24.5M-point grid.
    """
    if caps["-undefine_val"]:
        return ["-undefine_val", f"-1e30:{SENTINEL_THRESHOLD}"]
    if wgrib2 in _MASK_ARGS_CACHE:
        return _MASK_ARGS_CACHE[wgrib2]

    log("this wgrib2 has no -undefine_val; working out an older way to mask the sentinels")
    if field_stats(wgrib2, probe_sample)["min"] >= SENTINEL_THRESHOLD:
        # Nothing in this file to mask, so the mask cannot be proven on it -- and
        # equally cannot do any harm.  The per-file check after conversion still
        # enforces that no negative values survive.
        args = [a.format(lo=SENTINEL_THRESHOLD) for a in wgrib2_caps.MASK_CANDIDATES[1][1]]
        log(f"  WARNING: {os.path.basename(probe_sample)} contains no sentinel values, "
            f"so the fallback mask could not be proven on it.  Using {' '.join(args)} "
            "unverified; the per-file check still rejects any negative output.")
    else:
        args = wgrib2_caps.probe_mask(wgrib2, probe_sample, threshold=SENTINEL_THRESHOLD)
        log("  using: " + " ".join(args) + " (verified on this file)")
    _MASK_ARGS_CACHE[wgrib2] = args
    return args


def convert_grib(wgrib2: str, src: str, dst: str, new_grid: list[str] | None,
                 interp: str, retag: str, hours: int, clip: str | None,
                 caps: dict[str, bool] | None = None,
                 grib_type: str = "c3") -> dict:
    """
    Mask sentinels, optionally clip, optionally regrid, write `dst`.

    Returns the before/after statistics so the caller can prove the mask and the
    interpolation actually did something.

    Every wgrib2 option used here is checked against `caps` first.  `caps` comes
    from probe_capabilities(), because older wgrib2 builds are missing several of
    these and a missing option is a fatal error, not a warning -- an old build
    stops at the first one it does not recognise with "unknown option" and nothing
    gets converted at all.
    """
    caps = caps or {opt: True for opt in WANTED_OPTIONS}
    work = src
    tmp_clip = None
    if clip:
        if not caps["-small_grib"]:
            raise RuntimeError(
                "this wgrib2 has no -small_grib, so it cannot cut a box out of the "
                "MRMS grid.  Without a box the native grid is 24,500,000 points, "
                "which GEMPAK will not take -- a newer wgrib2 is needed here."
            )
        # Trim the source to a lat/lon box.  With --native this is what keeps the
        # grid inside GEMPAK's limit; when regridding it is only a speed-up, since
        # budget interpolation over all 24.5M MRMS points is the slow part.
        lon_w, lon_e, lat_s, lat_n = [float(v) for v in clip.split(":")]
        tmp_clip = dst + ".clip.grib2"
        clip_cmd = [wgrib2, src]
        if caps["-set_grib_type"]:
            clip_cmd += ["-set_grib_type", "same"]
        clip_cmd += ["-small_grib", f"{lon_w}:{lon_e}", f"{lat_s}:{lat_n}", tmp_clip]
        run(clip_cmd, "wgrib2 -small_grib")
        work = tmp_clip

    # Read the "before" numbers from whatever is actually going into the
    # conversion.  Taking them from the unclipped source would compare a CONUS
    # minimum of -3 against a clipped region that legitimately has no sentinels
    # in it, and the mask check below would fail a perfectly good file.
    before = field_stats(wgrib2, work)

    cmd = [wgrib2]
    if caps["-set_ftime_mode"]:
        # Stops wgrib2 rewriting "0-24 hour acc" as "0-1 day acc", which would
        # change the name GEMPAK derives for the grid.  It is an init option, so
        # it has to come before the input file.  Builds without it still produce
        # correct data, just possibly a "1 day" time code -- see the note logged
        # by probe_capabilities().
        cmd += ["-set_ftime_mode", "1"]
    cmd += [work]
    # Probe against `src`, the full CONUS download, not the clipped `work`.
    cmd += mask_args(wgrib2, caps, src)
    cmd += retag_args(retag, hours, caps)
    if caps["-set_grib_type"]:
        cmd += ["-set_grib_type", grib_type]
    if new_grid:
        cmd += ["-new_grid_winds", "earth",
                "-new_grid_interpolation", interp,
                "-new_grid", *new_grid, dst]
    else:
        cmd += ["-grib_out", dst]
    run(cmd, "wgrib2 conversion")
    if tmp_clip and os.path.exists(tmp_clip):
        os.unlink(tmp_clip)

    after = field_stats(wgrib2, dst)

    # Positive control on the mask.  Two separate things have to be true, because
    # either one alone can pass while the mask does nothing:
    #   * the minimum must no longer be negative, and
    #   * if the source DID contain sentinels, the undefined count must have gone
    #     up.  (wgrib2's "ndata" is the total point count and does not move when
    #     points are masked -- "undef" is the one that does.)
    # A source with min >= 0 genuinely has nothing to mask, which happens on a
    # small clip over an area with full radar coverage; that is not a failure.
    if after["min"] < SENTINEL_THRESHOLD:
        raise RuntimeError(
            f"{dst}: minimum is still {after['min']} after masking -- the sentinel "
            "mask did not take effect.  Do not use this file; report it."
        )
    if before["min"] < SENTINEL_THRESHOLD and after["undef"] <= before["undef"]:
        raise RuntimeError(
            f"{dst}: the source had values down to {before['min']} but the masked "
            f"output has no more undefined points than the input ({after['undef']}). "
            "The mask did not fire -- do not use this file; report it."
        )
    return {"before": before, "after": after}


# --------------------------------------------------------------------------- #
# the GEMPAK stage
# --------------------------------------------------------------------------- #
# Written with the full words "run" and "exit" rather than the "r"/"e"
# abbreviations, and with no leading whitespace, which is how Unidata's own GEMPAK
# driver scripts do it.  The blank line after "run" answers the prompt GEMPAK
# redisplays when the run finishes; "exit" then ends the program rather than
# leaving it sitting at that prompt waiting for input that will never come.
NAGRIB2_TEMPLATE = """\
GBFILE   = {gbfile}
INDXFL   =
GDOUTF   = {gdoutf}
PROJ     =
GRDAREA  =
KXKY     =
MAXGRD   = {maxgrd}
CPYFIL   = gds
GAREA    = grid
OUTPUT   = T
G2TBLS   = {g2tbls}
G2DIAG   =
OVERWR   = {overwr}
PDSEXT   = NO
run

exit
"""


def nagrib2_deck(grib: str, gemfile: str, maxgrd: int, g2tbls: str = "",
                 overwrite: bool = False) -> str:
    """The exact input nagrib2 is fed.  Printable, so it can be run by hand."""
    return NAGRIB2_TEMPLATE.format(
        gbfile=grib, gdoutf=gemfile, maxgrd=maxgrd,
        g2tbls=g2tbls, overwr="YES" if overwrite else "NO",
    )


def run_nagrib2(nagrib2: str, grib: str, gemfile: str, maxgrd: int,
                g2tbls: str = "", overwrite: bool = False,
                gemenviron: str | None = None, show_deck: bool = False,
                timeout: float = 180.0) -> str:
    """
    Drive nagrib2 non-interactively.

    CPYFIL=gds tells nagrib2 to take the navigation straight from the GRIB2 grid
    definition section, which is what you want here: the grid already IS the
    target grid, so nothing should be re-navigated.  Appending to an existing
    .gem file only works when that file's navigation matches.

    `gemenviron` is the path to a GEMPAK environment script (Gemenviron.profile).
    When given, nagrib2 is run inside a shell that sources it first.  A GEMPAK
    install that only works after sourcing its environment will otherwise fail
    here with something unhelpful about tables or $GEMTBL, even though running
    the same command by hand in a sourced shell works fine.

    The deck is echoed on any failure (and always with show_deck) so it can be
    pasted into nagrib2 by hand -- the conversion up to this point has already
    produced a GEMPAK-ready GRIB2 file, and a broken GEMPAK environment should
    not make that work unreachable.
    """
    script = nagrib2_deck(grib, gemfile, maxgrd, g2tbls, overwrite)
    # Always leave the deck on disk next to the GRIB2.  If anything goes wrong with
    # GEMPAK the useful next step is "nagrib2 < thatfile" in a shell where GEMPAK
    # works, and retyping fifteen lines from an error message is nobody's idea of
    # a good time.
    deck_path = grib + ".nagrib2.deck"
    try:
        with open(deck_path, "w") as fh:
            fh.write(script)
    except OSError:
        deck_path = ""
    if show_deck:
        log("nagrib2 input deck:\n" + script)

    # Remember the state of the output file.  "Did the file appear?" is not a
    # sufficient test when appending into an existing .gem -- the file is already
    # there, so a run that writes nothing would look like a success.
    was = os.stat(gemfile) if os.path.exists(gemfile) else None

    if gemenviron:
        # "." not "source": /bin/sh is not bash.  Sourcing noise is discarded so
        # it cannot be mistaken for nagrib2's own output.
        cmd = ["/bin/sh", "-c",
               f'. "{gemenviron}" >/dev/null 2>&1; exec "{nagrib2}"']
    else:
        cmd = [nagrib2]

    # Stream nagrib2's output as it arrives rather than capturing it silently.
    # GEMPAK programs are interactive: if one stops to ask something the deck does
    # not answer, a captured-output run looks like a hang with nothing on screen
    # and no way to tell what it is waiting for.  Echoing every line means the
    # last thing printed is the question it is stuck on.
    lines: list[str] = []
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    assert proc.stdin and proc.stdout
    try:
        proc.stdin.write(script)
        proc.stdin.close()          # EOF, so nagrib2 knows there is no more input
    except BrokenPipeError:
        pass                        # it exited before reading the deck; output says why

    killed: list[bool] = []

    def give_up() -> None:
        killed.append(True)
        proc.kill()

    timer = threading.Timer(timeout, give_up)
    timer.start()
    try:
        for line in proc.stdout:
            lines.append(line.rstrip("\n"))
            print(f"  [nagrib2] {line.rstrip()}", flush=True)
        proc.wait()
    finally:
        timer.cancel()

    text = "\n".join(lines).strip()
    if killed:
        raise RuntimeError(
            f"nagrib2 produced no more output for {timeout}s and was stopped.  It is "
            "almost certainly sitting at a prompt waiting for an answer the input "
            "deck does not give it -- the last [nagrib2] line above is the question.\n"
            "  Run the GEMPAK step by hand, in a shell where GEMPAK works:\n"
            f"    {nagrib2} < {deck_path or '(deck could not be written)'}\n"
            "  That way you see the prompts and can answer them.  The converted\n"
            f"  GRIB2 it is reading is ready at:\n    {grib}\n"
        )

    def fail(reason: str) -> RuntimeError:
        return RuntimeError(
            f"{reason}\n"
            f"  nagrib2 output:\n    " + (text.replace("\n", "\n    ") or "(nothing)") +
            "\n  The regridded GRIB2 is still there and is ready for GEMPAK:\n"
            f"    {grib}\n"
            "  To run the GEMPAK step by hand, source your GEMPAK environment and:\n"
            f"    {nagrib2} < {deck_path or '(deck could not be written)'}\n"
        )

    # nagrib2 can exit 0 having written nothing, so the file is the real test.
    if proc.returncode != 0:
        raise fail(f"nagrib2 exited {proc.returncode} for {grib}")
    if not os.path.exists(gemfile):
        raise fail(f"nagrib2 reported no error but {gemfile} was not created")
    now = os.stat(gemfile)
    if now.st_size == 0:
        raise fail(f"nagrib2 created {gemfile} but it is empty")
    if was is not None and (now.st_size, now.st_mtime_ns) == (was.st_size, was.st_mtime_ns):
        raise fail(f"nagrib2 exited cleanly but did not modify {gemfile} -- no grid "
                   "was written.  A navigation mismatch with the existing file, or a "
                   "parameter nagrib2 could not name, will do this quietly")
    return text


# --------------------------------------------------------------------------- #
# clip box helpers
# --------------------------------------------------------------------------- #
def auto_clip_box(wgrib2: str, template: str, margin: float = 0.5) -> str:
    """
    Lat/lon bounding box of a template grid, with a margin, as lonW:lonE:latS:latN.

    A Lambert grid's edges bow in lat/lon, so the four corners are not the
    bounding box -- this walks the whole perimeter (coarsely) and takes the
    extremes, then pads by `margin` degrees so the interpolation never reaches
    for a source point that was clipped away.
    """
    grid = grid_spec.read_grid(template)
    nx, ny = grid.fields["nx"], grid.fields["ny"]
    step_x = max(1, nx // 40)
    step_y = max(1, ny // 40)
    points: list[tuple[int, int]] = []
    for i in range(1, nx + 1, step_x):
        points += [(i, 1), (i, ny)]
    for j in range(1, ny + 1, step_y):
        points += [(1, j), (nx, j)]
    points += [(nx, ny)]

    args = [wgrib2, "-d", "1"]
    for i, j in points:
        args += ["-ijlat", str(i), str(j)]
    args += [template]
    out = run(args, "wgrib2 -ijlat").stdout

    lons = [float(v) for v in re.findall(r"lon=([-\d.]+)", out)]
    lats = [float(v) for v in re.findall(r"lat=([-\d.]+)", out)]
    if not lons or not lats:
        raise RuntimeError(f"could not read grid point lat/lons from wgrib2:\n{out}")

    # Intersect with the MRMS CONUS domain -- a target grid that reaches past
    # MRMS (the NBM CONUS grid does, at both ends) would otherwise produce a
    # "clip" box wider than the data, which clips nothing and just costs a pass.
    lon_w = max(min(lons) - margin, MRMS_CONUS_DOMAIN[0])
    lon_e = min(max(lons) + margin, MRMS_CONUS_DOMAIN[1])
    lat_s = max(min(lats) - margin, MRMS_CONUS_DOMAIN[2])
    lat_n = min(max(lats) + margin, MRMS_CONUS_DOMAIN[3])
    covers_all = (lon_w <= MRMS_CONUS_DOMAIN[0] and lon_e >= MRMS_CONUS_DOMAIN[1]
                  and lat_s <= MRMS_CONUS_DOMAIN[2] and lat_n >= MRMS_CONUS_DOMAIN[3])
    if covers_all:
        return ""      # nothing to gain; caller skips the clip
    return f"{lon_w:.4f}:{lon_e:.4f}:{lat_s:.4f}:{lat_n:.4f}"


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def build_times(args) -> list[dt.datetime]:
    if args.valid:
        return [parse_time(v) for v in args.valid]
    if args.times_file:
        with open(args.times_file) as fh:
            return [parse_time(line) for line in fh if line.strip() and not line.startswith("#")]
    start, end = args.start, args.end or args.start
    out, cur = [], start
    while cur <= end:
        out.append(cur)
        cur += dt.timedelta(hours=args.step)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mrms2gem.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Convert MRMS QPE GRIB2 from AWS into GEMPAK .gem grids.",
        epilog=textwrap.dedent("""\
            examples
              # one day of 24-hour Pass-2 QPE onto your NBM grid, one .gem per day
              ./mrms2gem.py --duration 24H --start 2026-09-09T00 --end 2026-09-09T23 \\
                  --grid-template nbm_co.grib2 --gemfile mrms_qpe24_%Y%m%d.gem

              # 6, 12 and 24 hour in one run, 12Z only
              ./mrms2gem.py --duration 06H --duration 12H --duration 24H \\
                  --valid 2026-09-09T12 --grid-template nbm_co.grib2 \\
                  --gemfile mrms_%Y%m%d_{dur}.gem

              # what is on S3 for that day, download nothing
              ./mrms2gem.py --duration 24H --start 2026-09-09 --list-only
            """),
    )
    p.add_argument("--duration", action="append", default=[], type=parse_duration,
                   metavar="06H|12H|24H|...",
                   help="accumulation length: 01H 03H 06H 12H 24H 48H 72H. "
                        "Plain numbers are fine too (--duration 6). Repeat the option "
                        "for several in one run (default 24H)")
    p.add_argument("--pass", dest="qpe_pass", type=int, default=2, choices=(1, 2),
                   help="MRMS QPE pass (default 2, the gauge-corrected one)")
    p.add_argument("--source", default="MultiSensor", choices=("MultiSensor", "RadarOnly"))
    p.add_argument("--start", type=parse_time, help="first valid time, UTC")
    p.add_argument("--end", type=parse_time, help="last valid time, UTC (default = start)")
    p.add_argument("--step", type=int, default=1, help="hours between valid times (default 1)")
    p.add_argument("--valid", action="append", help="explicit valid time; repeatable")
    p.add_argument("--times-file", help="file of valid times, one per line")
    p.add_argument("--tolerance", type=int, default=30,
                   help="minutes either side of the requested time to accept (default 30)")

    g = p.add_argument_group("target grid")
    g.add_argument("--grid-template", help="GRIB2 file whose grid MRMS should be put on "
                                          "(e.g. one of your NBM files)")
    g.add_argument("--grid", help="explicit wgrib2 -new_grid spec, three space-separated fields")
    g.add_argument("--native", action="store_true",
                   help="do not regrid.  Note: the native 24.5M-point MRMS grid will "
                        "almost certainly exceed your GEMPAK's maximum grid size")
    g.add_argument("--interp", default="budget",
                   choices=("budget", "bilinear", "bicubic", "neighbor", "neighbor-budget"),
                   help="wgrib2 interpolation (default budget, correct for accumulations)")
    g.add_argument("--clip", default="auto",
                   help="lonW:lonE:latS:latN to cut MRMS down to.  With --native this is "
                        "how you get a grid small enough for GEMPAK; with a target grid it "
                        "is only a speed-up.  'auto' derives it from the template "
                        "(default), 'none' disables it")

    o = p.add_argument_group("output")
    o.add_argument("--gemfile", default="mrms_{dur}_%Y%m%d.gem",
                   help="output .gem path; strftime codes and {dur} are substituted "
                        "(default mrms_{dur}_%%Y%%m%%d.gem)")
    o.add_argument("--retag", default="apcp-acc", choices=("apcp-acc", "apcp", "none"),
                   help="how to make GEMPAK recognise the parameter (default apcp-acc)")
    o.add_argument("--maxgrd", type=int, default=2000, help="nagrib2 MAXGRD (default 2000)")
    o.add_argument("--g2tbls", default="", help="nagrib2 G2TBLS value, for --retag none")
    o.add_argument("--workdir", default="./mrms_work", help="scratch for downloads and GRIB2")
    o.add_argument("--keep-grib", action="store_true", help="keep the regridded GRIB2 files")
    o.add_argument("--no-gempak", action="store_true", help="stop after the GRIB2 stage")
    o.add_argument("--overwrite-gem", action="store_true",
                   help="pass OVERWR=YES to nagrib2 (replaces matching grids)")
    o.add_argument("--list-only", action="store_true", help="list what would be fetched, then stop")
    o.add_argument("--wgrib2", help="path to wgrib2 (default: $WGRIB2 or PATH)")
    o.add_argument("--nagrib2", help="path to nagrib2 (default: $GEMEXE/nagrib2 or PATH)")
    o.add_argument("--gemenviron",
                   help="path to your GEMPAK environment script (Gemenviron.profile). "
                        "nagrib2 is then run inside a shell that sources it first -- use "
                        "this if GEMPAK only works in a shell where you have sourced it")
    o.add_argument("--gempak-timeout", type=float, default=180.0, metavar="SECONDS",
                   help="give up if nagrib2 goes this long without output (default 180). "
                        "A GEMPAK program stuck at a prompt would otherwise hang forever")
    o.add_argument("--show-deck", action="store_true",
                   help="print the nagrib2 input deck, so you can run the GEMPAK step "
                        "by hand if your GEMPAK environment is awkward")

    args = p.parse_args(glue_negative_values(list(sys.argv[1:] if argv is None else argv)))
    if not args.duration:
        args.duration = ["24H"]
    if not (args.start or args.valid or args.times_file):
        p.error("give --start (optionally with --end) or --valid or --times-file")
    if not (args.grid_template or args.grid or args.native or args.list_only):
        p.error("give --grid-template FILE, or --grid SPEC, or --native (see --help)")

    times = build_times(args)
    log(f"{len(times)} valid time(s) x {len(args.duration)} duration(s)")

    # ---- resolve the target grid once -------------------------------------- #
    wgrib2 = None
    new_grid = None
    clip = None
    caps: dict[str, bool] = {}
    retag = args.retag
    if not args.list_only:
        wgrib2 = find_exe(args.wgrib2 or "wgrib2", "WGRIB2",
                          "Install the wgrib2 package or set WGRIB2=/path/to/wgrib2.")
        caps = probe_capabilities(wgrib2)
        retag = resolve_retag(args.retag, caps)
        if args.grid:
            new_grid = args.grid.split()
            if len(new_grid) != 3:
                p.error("--grid needs exactly three fields, e.g. "
                        "'lambert:265:25:25:25 233.7234:2345:2539.703 19.229:1597:2539.703'")
        elif args.grid_template:
            grid = grid_spec.read_grid(args.grid_template)
            new_grid = grid.new_grid_args()
            log(f"target grid from {os.path.basename(args.grid_template)}: {grid.describe()}")
            log(f"  -new_grid {' '.join(new_grid)}")
        if new_grid and not wgrib2_has_interpolation(wgrib2):
            raise SystemExit(
                f"{wgrib2} was built without the interpolation (ipolates) library, so it "
                "cannot do -new_grid.\nRebuild it with a Fortran compiler present and "
                "USE_IPOLATES=3 in the makefile, or install a distro wgrib2 package that "
                "has it.  Check with: wgrib2 -config | grep -i interpolation"
            )
        # --clip is independent of regridding.  With --native it is the whole
        # point: it cuts a native-resolution MRMS subgrid small enough for GEMPAK
        # to hold, with no interpolation and no template file involved.
        if args.clip not in (None, "none"):
            if args.clip == "auto":
                clip = auto_clip_box(wgrib2, args.grid_template) if args.grid_template else None
                if clip == "":
                    log("clip auto: target grid spans the whole MRMS domain, not clipping")
                    clip = None
            else:
                try:
                    clip = normalise_clip(args.clip)
                except ValueError as exc:
                    p.error(str(exc))
            if clip:
                log(f"clipping MRMS to {clip}"
                    + (" before interpolating" if new_grid else " (native resolution)"))
                lon_w, lon_e, lat_s, lat_n = [float(v) for v in clip.split(":")]
                nx = round((lon_e - lon_w) / 0.01) + 1
                ny = round((lat_n - lat_s) / 0.01) + 1
                # Approximate: the box edges rarely land exactly on grid points,
                # so this is for the size warning below, not a promise.  The exact
                # count appears in the per-file line once the clip has run.
                log(f"  that box is roughly {nx} x {ny}, about {nx * ny / 1e6:.2f}M MRMS points")
                if not new_grid and nx * ny > GEMPAK_GRID_WARN:
                    log(f"  WARNING: that is more than {GEMPAK_GRID_WARN:,} points. "
                        "GEMPAK has a compiled-in maximum grid size (LLMXGD); if this "
                        "is over it, nagrib2 will refuse the grid.  Check yours with: "
                        "grep LLMXGD $GEMINC/GEMPRM.PRM")
        if args.native and not clip:
            log(f"WARNING: --native with no --clip means the full {MRMS_CONUS_POINTS:,}-point "
                "MRMS grid.  That is far past the grid size any stock GEMPAK is built for; "
                "expect nagrib2 to refuse it.  Add --clip lonW:lonE:latS:latN.")

    nagrib2 = None
    if not (args.list_only or args.no_gempak):
        gemexe = os.environ.get("GEMEXE")
        cand = args.nagrib2 or (os.path.join(gemexe, "nagrib2") if gemexe else "nagrib2")
        if args.gemenviron:
            # Do not insist on finding it now: the whole point of --gemenviron is
            # that nagrib2 only becomes visible once that script has been sourced.
            if not os.path.exists(args.gemenviron):
                p.error(f"--gemenviron {args.gemenviron} does not exist")
            nagrib2 = cand
            log(f"nagrib2 will be run as '{cand}' in a shell that sources "
                f"{args.gemenviron}")
        else:
            nagrib2 = find_exe(cand, None,
                               "Source your GEMPAK environment first (e.g. "
                               ". $NAWIPS/Gemenviron.profile), or pass "
                               "--gemenviron /path/to/Gemenviron.profile to have this "
                               "script source it for you, or --nagrib2 /path/to/nagrib2, "
                               "or --no-gempak to stop after the GRIB2 stage.")

    os.makedirs(args.workdir, exist_ok=True)
    made: list[str] = []
    failures: list[str] = []

    for duration in args.duration:
        product = mrms_s3.product_name(duration, args.source, args.qpe_pass)
        hours = DURATION_HOURS[duration]
        for want in times:
            tag = f"{product} {want:%Y-%m-%d %H:%M}Z"
            try:
                hit = mrms_s3.find_nearest(product, want, tolerance_min=args.tolerance)
                if not hit:
                    log(f"MISS  {tag}: nothing within {args.tolerance} min on S3")
                    failures.append(tag)
                    continue
                stamp, key, size = hit
                if args.list_only:
                    log(f"would fetch {stamp:%Y-%m-%d %H:%M}Z  {size/1e6:.2f} MB  s3://{mrms_s3.BUCKET}/{key}")
                    continue

                raw = mrms_s3.download(key, os.path.join(args.workdir, "raw"))
                regridded = os.path.join(
                    args.workdir, f"{product}_{stamp:%Y%m%d-%H%M%S}_{'native' if args.native else 'target'}.grib2")
                stats = convert_grib(wgrib2, raw, regridded,
                                     None if args.native else new_grid,
                                     args.interp, retag, hours, clip, caps)
                b, a = stats["before"], stats["after"]
                log(f"OK    {tag}: {int(b['ndata']):,} pts -> {int(a['ndata']):,} pts, "
                    f"masked {int(a['undef']):,}, mean {b['mean']:.3f} -> {a['mean']:.3f} mm, "
                    f"max {b['max']:.1f} -> {a['max']:.1f} mm")

                if not args.no_gempak:
                    gemfile = stamp.strftime(args.gemfile).replace("{dur}", duration.lower())
                    os.makedirs(os.path.dirname(os.path.abspath(gemfile)), exist_ok=True)
                    run_nagrib2(nagrib2, regridded, gemfile, args.maxgrd,
                                args.g2tbls, args.overwrite_gem,
                                gemenviron=args.gemenviron, show_deck=args.show_deck,
                                timeout=args.gempak_timeout)
                    log(f"      -> {gemfile}")
                    if gemfile not in made:
                        made.append(gemfile)
                if not args.keep_grib and not args.no_gempak:
                    os.unlink(regridded)
            except Exception as exc:           # one bad time must not kill the batch
                log(f"FAIL  {tag}: {exc}")
                failures.append(tag)

    if made:
        log("gem file(s): " + ", ".join(made))
    if failures:
        log(f"{len(failures)} of {len(times) * len(args.duration)} failed: " + "; ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

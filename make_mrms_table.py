#!/usr/bin/env python3
"""
make_mrms_table.py -- build a GEMPAK GRIB2 parameter table that knows about MRMS.

Only needed if you run mrms2gem.py with --retag none, i.e. you want the MRMS
records to keep their own identifiers instead of being re-tagged as APCP.

MRMS records are stamped:

    discipline 209, centre 161, local table 1, category 6, parameter 41

None of GEMPAK's shipped tables has a row for that, which is why nagrib2 either
skips the field or gives it no name.  The fix is a row in a GRIB2 parameter table
handed to nagrib2 through its G2TBLS keyword.

The awkward part is that these tables are FIXED-WIDTH, and the column positions
differ slightly between GEMPAK releases.  Rather than guess, this script reads
one of your installed tables, works out where each column actually starts from a
real data row in it, and writes the new table using exactly those offsets.  Then
it reads its own output back and checks every field lands in the right column --
a table whose columns are off by one parses as garbage, silently.

    # typical use
    python3 make_mrms_table.py --source $GEMTBL/grid/g2varsncep1.tbl \\
                               --out tables/g2varsmrms1.tbl

    # then in mrms2gem.py
    --retag none --g2tbls "g2varswmo2.tbl;/abs/path/tables/g2varsmrms1.tbl"

nagrib2's G2TBLS takes up to four table names separated by semicolons:
WMO parameter table ; local (centre) parameter table ; WMO vertical coordinate
table ; local vertical coordinate table.  Leave a slot empty to keep the default.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

# The rows we want to add.  GNAM is the name GEMPAK will show and that you use in
# GFUNC.  Keep them 4 characters or fewer plus a digit pair so they cannot collide
# with the NBM's own P06M/P12M/P24M.
MRMS_ROWS = [
    # All of these identifiers were read off real files in s3://noaa-mrms-pds on
    # 2026-09-10 -- they are not guesses.  Every one is discipline 209, centre
    # 161, local table 1, category 6; only the parameter number changes with the
    # accumulation length, which is why one table row per duration is needed.
    # disc cat parm pdt  description                               units      gnam  scale  missing  hzremap direction
    (209, 6, 37, 0, "MRMS MultiSensor QPE 1 hour Pass 2", "kg/m**2", "QP01", 0, -9999.00, 0, 0),
    (209, 6, 38, 0, "MRMS MultiSensor QPE 3 hour Pass 2", "kg/m**2", "QP03", 0, -9999.00, 0, 0),
    (209, 6, 39, 0, "MRMS MultiSensor QPE 6 hour Pass 2", "kg/m**2", "QP06", 0, -9999.00, 0, 0),
    (209, 6, 40, 0, "MRMS MultiSensor QPE 12 hour Pass 2", "kg/m**2", "QP12", 0, -9999.00, 0, 0),
    (209, 6, 41, 0, "MRMS MultiSensor QPE 24 hour Pass 2", "kg/m**2", "QP24", 0, -9999.00, 0, 0),
    (209, 6, 42, 0, "MRMS MultiSensor QPE 48 hour Pass 2", "kg/m**2", "QP48", 0, -9999.00, 0, 0),
    (209, 6, 43, 0, "MRMS MultiSensor QPE 72 hour Pass 2", "kg/m**2", "QP72", 0, -9999.00, 0, 0),
    (209, 6, 34, 0, "MRMS MultiSensor QPE 24 hour Pass 1", "kg/m**2", "QA24", 0, -9999.00, 0, 0),
    (209, 6, 6, 0, "MRMS RadarOnly QPE 24 hour", "kg/m**2", "QR24", 0, -9999.00, 0, 0),
]

NUM_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s")


NUMERIC_FIELDS = ("disc", "cat", "parm", "pdt", "scale", "missing", "hzremap", "direction")
TEXT_FIELDS = ("desc", "units", "gnam")
FIELD_ORDER = ("disc", "cat", "parm", "pdt", "desc", "units", "gnam",
               "scale", "missing", "hzremap", "direction")


def find_sample_row(lines: list[str]) -> tuple[str, dict[str, tuple[int, int]]]:
    """
    Find a real data row in an existing table and learn each column's span.

    Rows start with four integers (discipline, category, parameter, PDT); the
    description may contain spaces, so the trailing fields are identified from
    the right: direction, hzremap, missing, scale, GNAM, units.  Returns
    {field: (start, end)} character spans taken from that row, so the rows we
    write line up with whatever layout this GEMPAK release actually uses.
    """
    for raw in lines:
        line = raw.rstrip("\n")
        if line.startswith("!") or not line.strip():
            continue
        m = NUM_RE.match(line)
        if not m:
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        tail = parts[-6:]                     # units gnam scale missing hzremap direction
        try:
            float(tail[2]); float(tail[3]); int(tail[4]); int(tail[5])
        except ValueError:
            continue

        spans: dict[str, tuple[int, int]] = {}
        cursor = len(line)
        ok = True
        for name, value in (("direction", tail[5]), ("hzremap", tail[4]),
                            ("missing", tail[3]), ("scale", tail[2]),
                            ("gnam", tail[1]), ("units", tail[0])):
            idx = line.rfind(value, 0, cursor)
            if idx < 0:
                ok = False
                break
            spans[name] = (idx, idx + len(value))
            cursor = idx
        if not ok:
            continue

        for i, name in enumerate(("disc", "cat", "parm", "pdt")):
            spans[name] = (m.start(i + 1), m.end(i + 1))
        desc_start = m.end(4)
        while desc_start < len(line) and line[desc_start] == " ":
            desc_start += 1
        # A text column's WIDTH cannot be taken from the sample value -- the
        # sample's units might be "K", one character, while the column is twelve.
        # Each text column therefore runs up to the start of the next column.
        spans["desc"] = (desc_start, spans["units"][0] - 1)
        spans["units"] = (spans["units"][0], spans["gnam"][0] - 1)
        spans["gnam"] = (spans["gnam"][0], spans["scale"][0] - 1)
        return line, spans
    raise SystemExit(
        "could not find a usable data row in the source table -- check that you "
        "pointed --source at a GEMPAK g2vars*.tbl file."
    )


def place(spans: dict[str, tuple[int, int]], fields: dict[str, str]) -> str:
    """
    Write one row at the learned column spans.

    Numbers are right-aligned on the span's end column and text is left-aligned
    on its start column, which is how the shipped tables are laid out.  Text
    longer than the learned column is truncated rather than allowed to push the
    following columns sideways -- a row whose columns have shifted reads back as
    nonsense, and the check at the end of this script would reject it anyway.
    """
    width = max(end for _, end in spans.values()) + 4
    row = [" "] * width

    def put(start: int, text: str) -> None:
        row[start:start + len(text)] = list(text)

    for name in FIELD_ORDER:
        start, end = spans[name]
        text = fields[name]
        if name in NUMERIC_FIELDS:
            put(max(0, end - len(text)), text)
        else:
            limit = max(1, end - start)
            put(start, text[:limit])
    return "".join(row).rstrip()


def parse_row(line: str) -> tuple:
    """Read a row back the same way GEMPAK-style whitespace parsing would."""
    parts = line.split()
    ints = [int(v) for v in parts[:4]]
    tail = parts[-6:]
    desc = " ".join(parts[4:-6])
    return (*ints, desc, tail[0], tail[1], float(tail[2]), float(tail[3]),
            int(tail[4]), int(tail[5]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    help="an existing GEMPAK table to copy the column layout from, "
                         "e.g. $GEMTBL/grid/g2varsncep1.tbl")
    ap.add_argument("--out", default="tables/g2varsmrms1.tbl", help="table to write")
    args = ap.parse_args(argv)

    if not os.path.exists(args.source):
        raise SystemExit(f"{args.source} does not exist -- is GEMTBL set?")
    with open(args.source, errors="replace") as fh:
        lines = fh.readlines()

    sample, spans = find_sample_row(lines)
    print("layout learned from:")
    print("  " + sample)
    print("  columns: " + ", ".join(f"{k}@{v[0]}-{v[1]}" for k, v in
                                    sorted(spans.items(), key=lambda kv: kv[1][0])))

    header = [ln for ln in lines if ln.startswith("!")]
    out_lines = header + ["!\n", "! MRMS entries added by make_mrms_table.py\n", "!\n"]
    expected: list[tuple] = []
    truncated: list[str] = []
    for disc, cat, parm, pdt, desc, units, gnam, scale, missing, hz, direction in MRMS_ROWS:
        # What actually fits in this layout's text columns.  Anything wider would
        # push the following columns sideways, so it is cut instead -- and the
        # comparison below is made against the cut value, not the wish.
        widths = {name: max(1, spans[name][1] - spans[name][0]) for name in TEXT_FIELDS}
        fitted = {"desc": desc[:widths["desc"]].rstrip(),
                  "units": units[:widths["units"]].rstrip(),
                  "gnam": gnam[:widths["gnam"]].rstrip()}
        for name, full in (("desc", desc), ("units", units), ("gnam", gnam)):
            if fitted[name] != full:
                truncated.append(f"{gnam}: {name} {full!r} -> {fitted[name]!r} "
                                 f"({widths[name]} char column)")
        row = place(spans, {
            "disc": str(disc), "cat": str(cat), "parm": str(parm), "pdt": str(pdt),
            "scale": str(scale), "missing": f"{missing:.2f}",
            "hzremap": str(hz), "direction": str(direction), **fitted,
        })
        out_lines.append(row + "\n")
        expected.append((disc, cat, parm, pdt, fitted["desc"], fitted["units"],
                         fitted["gnam"], float(scale), float(missing), hz, direction))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.writelines(out_lines)

    # Read it back and check every field survived the column placement.
    if truncated:
        print("\nfields cut to fit this layout's columns:")
        for note in truncated:
            print("  " + note)

    print(f"\nwrote {args.out}; checking it parses back:")
    ok = True
    written = [ln for ln in open(args.out) if not ln.startswith("!") and ln.strip()]
    if len(written) != len(expected):
        raise SystemExit(f"wrote {len(expected)} rows but read back {len(written)}")
    for want, line in zip(expected, written):
        got = parse_row(line)
        if got != want:
            ok = False
            print(f"  BAD {line.rstrip()}\n      parsed {got}\n      wanted {want}")
        else:
            print(f"  ok  {line.rstrip()}")
    if not ok:
        raise SystemExit("the generated table does not read back correctly -- do not use it")
    print("\nall rows read back correctly.")
    print("use it with:  --retag none --g2tbls \"g2varswmo2.tbl;"
          f"{os.path.abspath(args.out)}\"")
    return 0


if __name__ == "__main__":
    sys.exit(main())

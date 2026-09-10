#!/usr/bin/env python3
"""
grid_spec.py -- read a GRIB2 file's grid definition and emit the wgrib2
`-new_grid` arguments that reproduce it exactly.

Why this exists: to do GEMPAK maths between MRMS and NBM the two grids must
share one navigation block, which means MRMS has to be interpolated onto the NBM
grid.  wgrib2 ships a helper (grid_defn.pl) that prints the -new_grid arguments
for a file, but it is not installed by every package and it works by scraping
wgrib2's human-readable -grid output.  This module reads Section 3 of the GRIB2
message directly instead, so the numbers come out at full stored precision
(angles are 1e-6 degree, grid lengths are 1e-3 m) with no text parsing and no
dependency beyond the standard library.

Supported grid templates: 0 (lat/lon), 10 (Mercator), 20 (polar stereographic),
30 (Lambert conformal).  Those cover MRMS (0), NBM CONUS/OCONUS (30 and 20) and
the HRRR (30).  Anything else raises, rather than guessing.

Run it directly to see the spec for a file:

    python3 grid_spec.py blend.t12z.core.f001.co.grib2
    lambert:265.000000:25.000000:25.000000:25.000000 \
        233.723400:2345:2539.703000 19.229000:1597:2539.703000
"""

from __future__ import annotations

import struct
import sys

MICRO = 1e-6   # GRIB2 stores angles in 1e-6 degree
MILLI = 1e-3   # ... and grid lengths in 1e-3 m


class GribGrid:
    """The parsed Section 3 of the first GRIB2 message in a file."""

    def __init__(self, template: int, npts: int, fields: dict):
        self.template = template
        self.npts = npts
        self.fields = fields

    def __repr__(self) -> str:
        return f"GribGrid(template={self.template}, npts={self.npts}, {self.fields})"

    # -- wgrib2 -new_grid argument triple -------------------------------------
    def new_grid_args(self) -> list[str]:
        f = self.fields
        if self.template == 0:
            # latlon lon0:nx:dlon lat0:ny:dlat   (dlon/dlat in degrees)
            return [
                "latlon",
                f"{f['lo1']:.6f}:{f['nx']}:{f['dlon']:.6f}",
                f"{f['la1']:.6f}:{f['ny']}:{f['dlat']:.6f}",
            ]
        if self.template == 30:
            # lambert:lov:latin1:latin2:lad lon0:nx:dx lat0:ny:dy   (dx/dy in m)
            return [
                f"lambert:{f['lov']:.6f}:{f['latin1']:.6f}:{f['latin2']:.6f}:{f['lad']:.6f}",
                f"{f['lo1']:.6f}:{f['nx']}:{f['dx']:.6f}",
                f"{f['la1']:.6f}:{f['ny']}:{f['dy']:.6f}",
            ]
        if self.template == 20:
            # nps/sps:lov:lad -- bit 1 of the projection centre flag set = South
            proj = "sps" if (f["proj_flag"] & 0x80) else "nps"
            return [
                f"{proj}:{f['lov']:.6f}:{f['lad']:.6f}",
                f"{f['lo1']:.6f}:{f['nx']}:{f['dx']:.6f}",
                f"{f['la1']:.6f}:{f['ny']}:{f['dy']:.6f}",
            ]
        if self.template == 10:
            return [
                f"mercator:{f['lad']:.6f}",
                f"{f['lo1']:.6f}:{f['nx']}:{f['dx']:.6f}",
                f"{f['la1']:.6f}:{f['ny']}:{f['dy']:.6f}",
            ]
        raise ValueError(
            f"grid template {self.template} is not one this tool can turn into a "
            "-new_grid spec.  Pass the target grid yourself with --grid, or use "
            "wgrib2's grid_defn.pl on the template file."
        )

    def describe(self) -> str:
        f = self.fields
        names = {0: "lat/lon", 10: "Mercator", 20: "polar stereographic",
                 30: "Lambert conformal"}
        head = f"{names.get(self.template, 'template ' + str(self.template))}" \
               f" {f.get('nx')} x {f.get('ny')} = {self.npts:,} points"
        if self.template == 0:
            return head + f", {f['dlon']:.4f} deg, first point {f['la1']:.4f}N {_lon180(f['lo1']):.4f}"
        return head + f", {f['dx']/1000:.3f} km, first point {f['la1']:.4f}N {_lon180(f['lo1']):.4f}"


def _lon180(lon: float) -> float:
    return lon - 360.0 if lon > 180.0 else lon


def _u4(b: bytes, i: int) -> int:
    return struct.unpack(">I", b[i:i + 4])[0]


def _s4(b: bytes, i: int) -> int:
    """
    GRIB2 signed values are sign-and-magnitude, not two's complement: the top
    bit is the sign.  Python's struct '>i' would get negative latitudes wrong.
    """
    v = struct.unpack(">I", b[i:i + 4])[0]
    return -(v & 0x7FFFFFFF) if v & 0x80000000 else v


def read_grid(path: str) -> GribGrid:
    """Parse Section 3 of the first GRIB2 message in `path`."""
    with open(path, "rb") as fh:
        data = fh.read(4_000_000)   # Section 3 is always near the start
    if data[:4] != b"GRIB":
        raise ValueError(f"{path} is not a GRIB2 file (no GRIB magic)")
    if data[7] != 2:
        raise ValueError(f"{path} is GRIB edition {data[7]}, this tool needs edition 2")

    pos = 16
    while pos < len(data) - 4:
        if data[pos:pos + 4] == b"7777":
            break
        seclen = _u4(data, pos)
        secnum = data[pos + 4]
        if seclen <= 0:
            raise ValueError(f"{path}: zero-length section at byte {pos}")
        if secnum == 3:
            return _parse_sec3(data[pos:pos + seclen])
        pos += seclen
    raise ValueError(f"{path}: no Section 3 (grid definition) found")


def _parse_sec3(s: bytes) -> GribGrid:
    npts = _u4(s, 6)
    template = struct.unpack(">H", s[12:14])[0]
    f: dict = {}
    if template == 0:
        f["nx"], f["ny"] = _u4(s, 30), _u4(s, 34)
        f["la1"], f["lo1"] = _s4(s, 46) * MICRO, _s4(s, 50) * MICRO
        f["la2"], f["lo2"] = _s4(s, 55) * MICRO, _s4(s, 59) * MICRO
        f["dlon"], f["dlat"] = _u4(s, 63) * MICRO, _u4(s, 67) * MICRO
        f["scan"] = s[71]
        # A lat/lon grid stored north-to-south has a positive Dj but the first
        # point is the NORTH edge.  wgrib2's -new_grid wants the SOUTH edge with
        # a positive increment, so flip when bit 2 of the scanning mode is clear.
        if not (f["scan"] & 0x40):
            f["la1"] = min(f["la1"], f["la2"])
    elif template == 30:
        f["nx"], f["ny"] = _u4(s, 30), _u4(s, 34)
        f["la1"], f["lo1"] = _s4(s, 38) * MICRO, _s4(s, 42) * MICRO
        f["lad"], f["lov"] = _s4(s, 47) * MICRO, _s4(s, 51) * MICRO
        f["dx"], f["dy"] = _u4(s, 55) * MILLI, _u4(s, 59) * MILLI
        f["proj_flag"], f["scan"] = s[63], s[64]
        f["latin1"], f["latin2"] = _s4(s, 65) * MICRO, _s4(s, 69) * MICRO
    elif template == 20:
        f["nx"], f["ny"] = _u4(s, 30), _u4(s, 34)
        f["la1"], f["lo1"] = _s4(s, 38) * MICRO, _s4(s, 42) * MICRO
        f["lad"], f["lov"] = _s4(s, 47) * MICRO, _s4(s, 51) * MICRO
        f["dx"], f["dy"] = _u4(s, 55) * MILLI, _u4(s, 59) * MILLI
        f["proj_flag"], f["scan"] = s[63], s[64]
    elif template == 10:
        f["nx"], f["ny"] = _u4(s, 30), _u4(s, 34)
        f["la1"], f["lo1"] = _s4(s, 38) * MICRO, _s4(s, 42) * MICRO
        f["lad"] = _s4(s, 47) * MICRO
        f["la2"], f["lo2"] = _s4(s, 51) * MICRO, _s4(s, 55) * MICRO
        f["scan"] = s[59]
        f["dx"], f["dy"] = _u4(s, 64) * MILLI, _u4(s, 68) * MILLI
    else:
        return GribGrid(template, npts, {})
    # Sanity: nx*ny should equal the point count in the section header.
    if f.get("nx") and f.get("ny") and f["nx"] * f["ny"] != npts:
        raise ValueError(
            f"grid parse disagrees with the file: nx*ny={f['nx']}*{f['ny']}"
            f"={f['nx']*f['ny']} but Section 3 says {npts} points"
        )
    return GribGrid(template, npts, f)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(f"usage: {sys.argv[0]} FILE.grib2 [--describe]")
    grid = read_grid(sys.argv[1])
    if "--describe" in sys.argv:
        print(grid.describe())
    else:
        print(" ".join(grid.new_grid_args()))

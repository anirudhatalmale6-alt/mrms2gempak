# MRMS QPE (AWS GRIB2) → GEMPAK .gem

Turns MRMS Quantitative Precipitation Estimate files from the public AWS archive
into GEMPAK grids **on a grid you can actually do arithmetic with** — specifically,
on the same navigation as your existing NBM grids, so `GDDIAG` maths between MRMS
and NBM works.

Runs unattended, takes a date/time range or a list of times, batches a whole day.

---

## 1. The two problems this solves

### 1.1 GEMPAK has no idea what MRMS is

Read a raw MRMS file and the parameter has no name:

```
$ wgrib2 MRMS_MultiSensor_QPE_24H_Pass2_00.00_20260909-120000.grib2
1:0:d=2026090912:var discipline=209 center=161 local_table=1 parmcat=6 parm=41: ...
```

MRMS uses **discipline 209, centre 161, local table 1, category 6**, and a
parameter number that changes with the accumulation length. GEMPAK's shipped
GRIB2 tables have no row for any of it, so `nagrib2` has nothing to call the
field. Parameter numbers, read off real files in the bucket on 2026-09-10:

| product | parameter | GNAM used by `make_mrms_table.py` |
| --- | --- | --- |
| `MultiSensor_QPE_01H_Pass2` | 37 | `QP01` |
| `MultiSensor_QPE_03H_Pass2` | 38 | `QP03` |
| `MultiSensor_QPE_06H_Pass2` | 39 | `QP06` |
| `MultiSensor_QPE_12H_Pass2` | 40 | `QP12` |
| `MultiSensor_QPE_24H_Pass2` | 41 | `QP24` |
| `MultiSensor_QPE_48H_Pass2` | 42 | `QP48` |
| `MultiSensor_QPE_72H_Pass2` | 43 | `QP72` |
| `MultiSensor_QPE_24H_Pass1` | 34 | `QA24` |
| `RadarOnly_QPE_24H` | 6 | `QR24` |

Two ways round it, both supported — see §5.

### 1.2 The grid is too big, and it is the wrong grid

```
$ wgrib2 -d 1 -grid MRMS_MultiSensor_QPE_24H_Pass2_00.00_20260909-120000.grib2
	lat-lon grid:(7000 x 3500) units 1e-06 ... #points=24500000
	lat 54.995000 to 20.005001 by 0.010000
	lon 230.004999 to 299.994997 by 0.010000
```

24,500,000 points on a 0.01° lat/lon grid. That is past the grid size GEMPAK is
compiled for, and it is a different projection and spacing from the NBM's 2.5 km
Lambert conformal grid:

```
$ wgrib2 -d 1 -grid blend.t12z.core.f001.co.grib2
	Lambert Conformal: (2345 x 1597) ... Dx 2539.703000 m Dy 2539.703000 m
	Lat1 19.229000 Lon1 233.723400 LoV 265.000000 LatD 25.000000 Latin1 25.0 Latin2 25.0
```

**Cutting the same lat/lon box out of both with `wgrib2 -small_grib` does not make
them compatible.** You get two grids over the same ground with different
navigation. GEMPAK cannot do arithmetic across different navigation, and
`nagrib2` will not write the second grid into a `.gem` file whose navigation does
not match.

What does work is interpolating MRMS onto the NBM grid. Point count drops from
24,500,000 to 3,744,965, the navigation block becomes identical, and the two
grids can live in one `.gem` file. `--interp budget` (the default) is wgrib2's
interpolation mode for accumulations: it conserves the areal total rather than
sampling a single point, which is what you want going from ~1 km to 2.5 km.

---

## 2. Dependencies

| what | why | check |
| --- | --- | --- |
| python3 ≥ 3.9 | the scripts. **Standard library only** — no boto3, no numpy, no aws CLI | `python3 -V` |
| wgrib2, **built with interpolation** | masking, re-tagging, `-new_grid` | `wgrib2 -config \| grep -i interpolation` |
| GEMPAK 7.x (`nagrib2` on `PATH`) | the `.gem` stage only | `which nagrib2` |

The wgrib2 check matters. `-new_grid` needs the `ipolates` library, which needs a
Fortran compiler at build time. A build that reports

```
interpolation package is not installed, default vectors: ...
```

decodes and writes GRIB2 perfectly but **cannot regrid**. `mrms2gem.py` tests for
this before it downloads anything and tells you, rather than failing part-way
through a batch. If your wgrib2 lacks it, rebuild from
<https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/wgrib2.tgz> with `gfortran`
installed and `USE_IPOLATES=3` in the makefile.

No AWS credentials are needed anywhere. `noaa-mrms-pds` is a public Open Data
bucket and the scripts read it over plain anonymous HTTPS.

---

## 3. Quick start

```bash
# what is on S3 for a day — downloads nothing
./mrms2gem.py --duration 24H --start 2026-09-09 --list-only

# one day of 24-hour Pass-2 QPE onto your NBM grid, hourly, one .gem per day
./mrms2gem.py --duration 24H \
              --start 2026-09-09T00 --end 2026-09-09T23 \
              --grid-template /data/nbm/blend.t12z.core.f001.co.grib2 \
              --gemfile /data/gempak/mrms_{dur}_%Y%m%d.gem

# 6, 12 and 24 hour at 12Z only
./mrms2gem.py --duration 06H --duration 12H --duration 24H \
              --valid 2026-09-09T12 \
              --grid-template /data/nbm/blend.t12z.core.f001.co.grib2 \
              --gemfile /data/gempak/mrms_{dur}_%Y%m%d.gem

# stop after the GRIB2 stage (no GEMPAK needed) and keep the files
./mrms2gem.py --duration 24H --valid 2026-09-09T12 \
              --grid-template nbm_co.grib2 --no-gempak --keep-grib
```

`--gemfile` takes `strftime` codes and `{dur}`. Times are **UTC**, and are
accepted as `2026-09-09T12`, `2026-09-09 12:00`, `2026090912`,
`20260909-1200`, …

Per run you get one line per file:

```
[mrms2gem] target grid from blend.t12z.core.f001.co.grib2: Lambert conformal 2345 x 1597 = 3,744,965 points, 2.540 km, first point 19.2290N -126.2766
[mrms2gem]   -new_grid lambert:265.000000:25.000000:25.000000:25.000000 233.723400:2345:2539.703000 19.229000:1597:2539.703000
[mrms2gem] OK    MultiSensor_QPE_24H_Pass2_00.00 2026-09-09 12:00Z: 24,500,000 pts -> 3,744,965 pts, masked 1,512,918, mean 1.686 -> 1.994 mm, max 339.8 -> 339.8 mm
```

The numbers are there to be read, not decoration. `masked` is how many points
were sentinels (§4); the mean moving 1.686 → 1.994 mm is those sentinels no
longer dragging the average down. A run that masked 0 points on a CONUS file
means the mask did not fire, and the script stops rather than hand you a grid
with negative rainfall in it.

### Useful options

| option | default | notes |
| --- | --- | --- |
| `--duration` | `24H` | repeatable: `01H 03H 06H 12H 24H 48H 72H` |
| `--pass` | `2` | Pass 2 is the gauge-corrected one |
| `--source` | `MultiSensor` | or `RadarOnly` |
| `--step` | `1` | hours between valid times |
| `--tolerance` | `30` | minutes either side of the requested time to accept |
| `--grid-template` | — | GRIB2 file whose grid to copy (one of your NBM files) |
| `--grid` | — | explicit wgrib2 `-new_grid` spec instead of a template |
| `--native` | off | no regridding. Will almost certainly blow GEMPAK's grid limit |
| `--interp` | `budget` | `budget` for accumulations; `neighbor` to keep exact values |
| `--clip` | `auto` | trim MRMS to the target's bounding box first, purely for speed |
| `--retag` | `apcp-acc` | see §5 |
| `--maxgrd` | `2000` | `nagrib2` MAXGRD |
| `--overwrite-gem` | off | `OVERWR=YES`: replaces matching grids |
| `--no-gempak` | off | stop after GRIB2 |
| `--keep-grib` | off | keep the regridded GRIB2 |

`--clip auto` derives the target grid's lat/lon bounding box by walking its
perimeter with `wgrib2 -ijlat` (a Lambert grid's edges bow, so the four corners
are not the bounding box) and trims MRMS to it before interpolating. For the full
NBM CONUS grid the box is wider than MRMS at both ends, so the script says so and
skips the clip instead of doing a pointless pass. It pays off for regional runs.

---

## 4. The sentinel values — why this is not optional

MRMS encodes **−3 as "outside radar coverage"** and **−999 as "missing"**. They
are ordinary values in the file, not flagged missing:

```
before:  ndata=24500000:undef=0:mean=1.68587:min=-3:max=339.8
```

Interpolate that as-is and every −3 gets averaged into the neighbouring output
boxes. A 2.5 km cell that straddles a coverage edge comes out negative, and a
verification statistic computed over it is quietly wrong. So before anything else
the script runs

```
wgrib2 … -undefine_val -1e30:-0.001 …
```

which turns them into proper GRIB2 missing data:

```
after:   ndata=24500000:undef=1512918:mean=1.99427:min=0:max=339.8
```

1,512,918 points flagged, minimum back to 0, mean up by 0.3 mm. The script then
checks the minimum is no longer negative and **aborts that file** if it is —
a mask that silently fails looks exactly like a mask that worked if you only
check the exit status.

---

## 5. Making GEMPAK recognise the parameter

`--retag` picks between three approaches. Default is `apcp-acc`.

### `--retag apcp-acc` (default, no GEMPAK table edits)

Re-stamps the record as standard WMO `APCP` at the surface **and** rebuilds the
product definition as an N-hour accumulation:

```
wgrib2 -set_ftime_mode 1 in.grib2 \
    -undefine_val -1e30:-0.001 \
    -set discipline 0 -set center 7 -set master_table 2 -set local_table 1 \
    -set table_4.1 1 -set table_4.2 8 -set_lev surface \
    -set table_1.2 1 -set_date -24hr -set_ftime2 "0-24 hour acc fcst" \
    -set_grib_type c3 … -new_grid … out.grib2
```

Result:

```
1:0:12Z08sep2026:APCP Total Precipitation [kg/m^2]:surface:0-24 hour acc fcst
```

Two details in there that are easy to get wrong:

* **`-set_date -24hr`.** An MRMS 24-hour accumulation is valid *at* the file's
  timestamp and covers the preceding 24 hours. Turning it into a "0-24 hour acc"
  record without moving the reference time back 24 hours would make it claim to
  be valid 24 hours in the *future* — a whole day of data misdated, with nothing
  in the output to show it.
* **`-set_ftime_mode 1`.** Without it wgrib2 helpfully rewrites
  `0-24 hour acc fcst` as `0-1 day acc fcst`, which changes the name GEMPAK
  derives for the grid.

This is the mode that should make GEMPAK name the grids `P06M` / `P12M` / `P24M`,
matching how your NBM precip grids are already named, distinguished from them by
reference time.

### `--retag apcp`

Same re-stamp to `APCP` at the surface, but left as an analysis at its own valid
time (PDT 0). Use this if you would rather the grid carry an analysis time tag
than look like a forecast.

### `--retag none` + a parameter table

Keeps MRMS's own identifiers and teaches GEMPAK about them. `make_mrms_table.py`
builds the table:

```bash
python3 make_mrms_table.py --source $GEMTBL/grid/g2varsncep1.tbl \
                           --out tables/g2varsmrms1.tbl

./mrms2gem.py --duration 24H --valid 2026-09-09T12 \
              --grid-template nbm_co.grib2 \
              --retag none \
              --g2tbls "g2varswmo2.tbl;$PWD/tables/g2varsmrms1.tbl"
```

These tables are fixed-width and the column positions differ between GEMPAK
releases, so the script does not guess: it reads a real data row out of *your*
table, works out where each column starts and ends, writes the MRMS rows at those
same columns, then reads its own output back and checks every field landed in the
right column. A table that is one column out parses as garbage without
complaining.

`nagrib2`'s `G2TBLS` takes up to four table names separated by `;` — WMO
parameter table, local parameter table, WMO vertical coordinate table, local
vertical coordinate table. Empty slots keep the defaults.

---

## 6. Doing the maths in GEMPAK

The reliable route is **both grids in one `.gem` file**, which is exactly what the
matching navigation buys you. Point `--gemfile` at a *copy* of your NBM file:

```bash
cp /data/gempak/nbm_20260908.gem /data/gempak/compare_20260909.gem
./mrms2gem.py --duration 24H --valid 2026-09-09T12 \
              --grid-template /data/nbm/blend.t12z.core.f001.co.grib2 \
              --gemfile /data/gempak/compare_20260909.gem
gdinfo GDFILE=/data/gempak/compare_20260909.gem LSTALL=YES OUTPUT=T
```

**Work on a copy.** `nagrib2` appends into the file you name; do not point it at
the only copy of a set of NBM grids until you have seen a run come out right.
If the navigation does not match, `nagrib2` refuses — which is the failure mode
you want, but read its output rather than assuming.

Then difference them:

```
gddiag
 GDFILE  = /data/gempak/compare_20260909.gem
 GDOUTF  = /data/gempak/compare_20260909.gem
 GFUNC   = SUB(P24M^20260909/1200,P24M^20260908/1200F024)
 GDATTIM = 20260909/1200
 GLEVEL  = 0
 GVCORD  = none
 GRDNAM  = QPEERR
 r
```

Exact grid identifiers come from `gdinfo`; the `^time` tags are what keep the
MRMS analysis and the NBM forecast apart when both are named `P24M`. With
`--retag none` the MRMS grid is `QP24` instead and no `^time` disambiguation is
needed.

---

## 7. Files

| file | what it does |
| --- | --- |
| `mrms2gem.py` | the driver: fetch → mask → regrid → `nagrib2` |
| `mrms_s3.py` | anonymous listing/fetching from `noaa-mrms-pds`; run it directly to browse a day |
| `grid_spec.py` | reads Section 3 of a GRIB2 file and emits the `-new_grid` spec; run it directly on a file |
| `make_mrms_table.py` | builds a GEMPAK GRIB2 parameter table with MRMS rows, matched to your layout |
| `wgrib2_caps.py` | reports what a wgrib2 build can do and proves the sentinel mask works on a real file |

Each is standalone and documents itself:

```bash
python3 mrms_s3.py 24H 2 2026-09-09          # list a day
python3 grid_spec.py nbm_co.grib2            # print the -new_grid spec
python3 grid_spec.py nbm_co.grib2 --describe
python3 wgrib2_caps.py wgrib2 sample_mrms.grib2
```

---

## 8. What is verified, and what is not

Verified by running it, against real files from the bucket:

* bucket layout, product names, and cadence — every QPE accumulation product
  (01H…72H, Pass 1 and Pass 2) is hourly, 24 files a day, ~4.6 MB gzipped at 24H;
  `PrecipRate` is 2-minute (717 files on 2026-09-09)
* anonymous listing and download, continuation tokens, gunzip, truncation check
* the MRMS identifiers in the table in §1.1
* the sentinel mask, with before/after counts (−3 → undefined, 1,512,918 points)
* the re-tagging to `APCP`, the reference-time shift, and the `0-24 hour acc
  fcst` product definition
* `grid_spec.py` against a real NBM CONUS file — its output matches wgrib2's own
  `-grid` report field for field
* `--clip` with `-small_grib`, and the `auto` box derivation
* `make_mrms_table.py` against two deliberately different table layouts,
  including the read-back check catching a mis-aligned column

**Not yet run end to end here:**

* `-new_grid` itself. The wgrib2 I built for testing has no Fortran compiler
  available, so it reports `interpolation package is not installed`. The command
  is assembled and the guard that detects this is tested; the interpolation needs
  a wgrib2 that has `ipolates`.
* the `nagrib2` stage, for the obvious reason that GEMPAK is not installed here.

So the first run on your machine is the one that settles the GEMPAK half — which
`--retag` mode gives the grid name you want, and whether `MAXGRD` and the
navigation match your existing files. Send me the `nagrib2` output from that run
and I will adjust.

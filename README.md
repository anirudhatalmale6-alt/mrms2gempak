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

Two ways round it, both supported — see §6.

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
| wgrib2 | masking, re-tagging, clipping | `wgrib2 --version` |
| …**built with interpolation** | only for `-new_grid`, i.e. regridding | `wgrib2 -config \| grep -i interpolation` |
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

### Older wgrib2 builds

wgrib2 has gained options over the years and a build that lacks one stops dead
with `*** FATAL ERROR: unknown option -set_ftime_mode ***` — it does not warn and
carry on, so one missing option converts nothing. Rather than require a
particular version, `mrms2gem.py` asks your binary what it has, reports it, and
assembles the command from what is there:

```
[mrms2gem] wgrib2: v3.1.3 10/2023 Wesley Ebisuzaki
[mrms2gem] wgrib2 does not have: -set_ftime_mode -- working around where possible
```

What each missing option costs:

| missing | consequence |
| --- | --- |
| `-set_ftime_mode` | time code may read `0-1 day` instead of `0-24 hour`; data unaffected |
| `-set_grib_type` | output GRIB2 is uncompressed, so larger (4.6 MB → 37 MB for CONUS) |
| `-set_lev` | level stays as MRMS wrote it rather than being set to surface |
| `-set_ftime2`/`-set_ftime` | cannot build an accumulation record → falls back to `--retag apcp` |
| `-set_date` | **blocks** `--retag apcp-acc` rather than working around it: without the reference-time shift every grid would be dated N hours late |
| `-set` | cannot re-tag at all → falls back to `--retag none` + parameter table |
| `-undefine_val` | falls back to `-rpn dup:-0.001:>=:mask`, proven on the file before use |
| `-small_grib` | **fatal** — no way to cut the grid down to a size GEMPAK accepts |

Every downgrade is logged with its reason. The one deliberate refusal is
`-set_date`: a wrong date on every grid is worse than a missing accumulation tag,
so that case loses the tag instead of the date.

---

## 3. Straight conversion, no NBM involved

If you would rather get MRMS into GEMPAK first and deal with remapping yourself,
skip the template entirely:

```bash
python3 mrms2gem.py --duration 24H --valid 2026-09-09T12 \
                    --native --clip=-105:-90:35:45 \
                    --gemfile ~/mrms_{dur}_%Y%m%d.gem
```

No interpolation, no NBM file, no `ipolates` — MRMS stays on its own 0.01°
lat/lon grid and only gets cut down to the box you asked for.

**The clip is not optional here, and that is the one real constraint.** Full
CONUS native is 24,500,000 points, far past the grid size any stock GEMPAK is
built for (`LLMXGD` in `GEMPRM.PRM`, fixed at compile time). If you cannot
rebuild GEMPAK, a box is the only way to keep native resolution. The script warns
before it wastes a download.

At 0.01° the arithmetic is easy — degrees × 100 per side:

| box | points |
| --- | --- |
| 10° × 5° | 500,000 |
| 15° × 10° | 1,500,000 |
| 20° × 15° | 3,000,000 |
| 30° × 20° | 6,000,000 |

To find what your build actually allows:

```bash
grep LLMXGD $GEMINC/GEMPRM.PRM
```

Note `--clip=-105:-90:35:45` with the **equals sign**. A value starting with `-`
looks like an option to the argument parser; the script rewrites the spaced form
for you, but the joined form is the habit worth having. Longitudes may be given
either as −105 or as 255.

## 4. Quick start

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
were sentinels (§5); the mean moving 1.686 → 1.994 mm is those sentinels no
longer dragging the average down. A run that masked 0 points on a CONUS file
means the mask did not fire, and the script stops rather than hand you a grid
with negative rainfall in it.

### Useful options

| option | default | notes |
| --- | --- | --- |
| `--duration` | `24H` | `01H 03H 06H 12H 24H 48H 72H`; plain numbers work too (`--duration 6`). Repeatable |
| `--pass` | `2` | Pass 2 is the gauge-corrected one |
| `--source` | `MultiSensor` | or `RadarOnly` |
| `--step` | `1` | hours between valid times |
| `--tolerance` | `30` | minutes either side of the requested time to accept |
| `--grid-template` | — | GRIB2 file whose grid to copy (one of your NBM files) |
| `--grid` | — | explicit wgrib2 `-new_grid` spec instead of a template |
| `--native` | off | no regridding — use with `--clip`, see §3 |
| `--interp` | `budget` | `budget` for accumulations; `neighbor` to keep exact values |
| `--clip` | `auto` | `lonW:lonE:latS:latN`. Required size control with `--native`; only a speed-up when regridding |
| `--retag` | `apcp-acc` | see §6 |
| `--maxgrd` | `2000` | `nagrib2` MAXGRD |
| `--gemenviron` | — | path to `Gemenviron.profile`; `nagrib2` is run in a shell that sources it |
| `--show-deck` | off | print the `nagrib2` input deck |
| `--overwrite-gem` | off | `OVERWR=YES`: replaces matching grids |
| `--no-gempak` | off | stop after GRIB2 |
| `--keep-grib` | off | keep the regridded GRIB2 |

`--clip auto` derives the target grid's lat/lon bounding box by walking its
perimeter with `wgrib2 -ijlat` (a Lambert grid's edges bow, so the four corners
are not the bounding box) and trims MRMS to it before interpolating. For the full
NBM CONUS grid the box is wider than MRMS at both ends, so the script says so and
skips the clip instead of doing a pointless pass. It pays off for regional runs.

### If your GEMPAK environment is awkward

A GEMPAK that only works in a shell where you have sourced its environment will
fail here with something unhelpful about tables or `$GEMTBL`. Two ways out:

```bash
# let the script source it for you before calling nagrib2
python3 mrms2gem.py ... --gemenviron $NAWIPS/Gemenviron.profile

# or see the deck and run the GEMPAK step yourself
python3 mrms2gem.py ... --show-deck
```

With `--gemenviron` the script does not insist on finding `nagrib2` up front —
that is the whole point, since it only becomes visible after sourcing.

Either way, **any** `nagrib2` failure prints the full input deck and the path of
the GRIB2 file that is already converted and waiting. Worst case the script's job
is "hand me a GEMPAK-ready GRIB2" and you drive `nagrib2` yourself; a fighting
GEMPAK environment does not put the rest of the work out of reach.

**nagrib2 does not exit when fed from a pipe.** Observed on GEMPAK 7 with
nagrib2 v3.0.2-era tables: it processes the deck, prints its report, returns to
`GEMPAK-NAGRIB2>` and sits there — `exit` at the end of the deck does not end it,
and without intervention the run never returns. So the script does not wait for
it to exit. It watches for

```
         1 grids were written to the GEMPAK file:
```

which is nagrib2 announcing the work is done, allows three seconds for trailing
output, then closes it down. That line is also how success is judged, rather than
by exit status or file timestamps:

| nagrib2 says | treated as |
| --- | --- |
| `N grids were written`, N > 0 | success, count logged |
| `0 grids were written` + `[GD -10] Grid already exists.` | success — the grid is already in the file; `--overwrite-gem` to replace it |
| `0 grids were written`, nothing about existing | failure, with the output and the deck |
| no report at all | falls back to checking the file changed |

That second row is the one that looks alarming and is not: re-running the same
time is a no-op because `OVERWR=NO`, and the desired grid is present either way.

Output is echoed live (prefixed `[nagrib2]`) so if it does stop somewhere
unexpected, the last line on screen is the prompt it is stuck on, and
`--gempak-timeout` (default 180s) is the backstop. The deck is always written
next to the GRIB2 as `<grib2>.nagrib2.deck`, so the manual route is one command:

```bash
nagrib2 < mrms_work/MultiSensor_QPE_24H_Pass2_00.00_20260909-120000_native.grib2.nagrib2.deck
```

`nagrib2` can exit 0 having written nothing — a navigation mismatch or a
parameter it cannot name both do that quietly. The script stats the output file
before and after and treats "exited cleanly but did not change the file" as a
failure, because otherwise an empty append reads as success.

---

## 5. The sentinel values — why this is not optional

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

## 6. Making GEMPAK recognise the parameter

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

## 7. Doing the maths in GEMPAK

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

## 8. Files

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

## 9. What is verified, and what is not

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
* the straight-conversion path of §3: `--native --clip=-105:-90:35:45` cuts
  1,500,000 points out of CONUS at native resolution, no interpolation involved
* `--clip` validation: a reversed box, a box outside the MRMS domain, and a
  spaced negative longitude (`--clip -105:-90:35:45`) are all handled
* `make_mrms_table.py` against two deliberately different table layouts,
  including the read-back check catching a mis-aligned column
* the older-wgrib2 handling, against a stand-in that rejects named options the
  same way a build without them does: each of `-set_ftime_mode`, `-set_ftime2`,
  `-set_ftime`, `-set_date`, `-set`, `-set_lev`, `-set_grib_type` and
  `-undefine_val` removed in turn, and several at once
* the `-rpn dup:-0.001:>=:mask` fallback, which flags the same 1,512,918 points
  and reaches the same minimum as `-undefine_val` on a real CONUS file
* the GEMPAK stage against stand-ins built from a real nagrib2 transcript,
  including the one that writes nothing because the grid already exists and then
  sits at its prompt: the run now finishes in ~8s instead of hanging, and the
  three zero/non-zero/no-report cases are told apart
* the GEMPAK stage against a stand-in `nagrib2`: `--gemenviron` really does
  source the environment (and the same run without it really does not — both
  directions checked), sourcing noise stays out of the reported output, and a
  `nagrib2` that exits 0 without touching the file is caught rather than
  reported as success

**Not yet run end to end here:**

* `-new_grid` itself. The wgrib2 I built for testing has no Fortran compiler
  available, so it reports `interpolation package is not installed`. The command
  is assembled and the guard that detects this is tested; the interpolation needs
  a wgrib2 that has `ipolates`.
* the **real** `nagrib2`, for the obvious reason that GEMPAK is not installed
  here. The stand-in above exercises the plumbing around it, not GEMPAK itself.

So the first run on your machine is the one that settles the GEMPAK half — which
`--retag` mode gives the grid name you want, and whether `MAXGRD` and the
navigation match your existing files. Send me the `nagrib2` output from that run
and I will adjust.

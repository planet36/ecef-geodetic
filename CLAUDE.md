# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A test suite that measures the accuracy and speed of about 75 ECEF-to-geodetic conversion
functions on the WGS 84 ellipsoid.  It backs the 2020 SIW paper in `docs/`, whose
recommendation is `olson_1996` (custom height) with 2 iterations for iterative algorithms.

## Commands

```sh
make -j $(nproc) input   # generate the input files in input/
make -j $(nproc)         # build both test binaries (and the inputs)
make full                # run both tests and combine them, about 1 minute
make acc1                # quick single-point accuracy check on the speed inputs
make lint                # clang-tidy with the repo's .clang-tidy
make clean-all           # remove binaries, .d/.opts files, and input files
```

There is no unit-test framework.  To exercise only some algorithms, pass their names (the
keys in `include/map_func_name_to_func_info.hpp`) as trailing arguments:

```sh
./test-ecef_to_geodetic-acc -v -a -g olson_1996 < input/geod.2d.region-all.txt
NUM_THREADS=1 ./test-ecef_to_geodetic-speed olson_1996 sedris < input/ecef.2d.speed.txt
```

The accuracy binary's flags are `-a` for the accuracy test, `-1` for the single-point test,
`-g` when stdin holds geodetic rather than ECEF points, `-t` for multiple threads, `-c` to
collect every error for precise statistics, `-m N` to skip algorithms whose
`ilog10_mean_dist_err` exceeds N, and `-s N` for N rounds of a built-in speed test, up to
10000.  In its JSON, an `ilog10_mean_dist_err` of 99 marks an inaccurate algorithm and -99 an
exact one.

Both binaries exit with an error message on a malformed input line, when they read no input
coordinates, and on an algorithm name that is unknown or given more than once.  The accuracy
binary also exits on an option value that is not an integer in range.

The speed binary is a Google Benchmark program that accepts the usual `--benchmark_*` flags and
requires ECEF input.  `NUM_THREADS` sets how many threads run each benchmark.  It defaults to 1,
0 means one per hardware thread, and a value that is not an integer from 0 to the number of
hardware threads is an error.

Each `make` run writes its test output into its own folder, `results/TIMESTAMP/`, named for
the time the run started.  `make acc`, `make acc1`, and `make speed` write `acc.json`,
`acc1.json`, and `speed.json` there and embed the compile options, which the build extracts
from the binary with `readelf` into `*.opts`.  Runs from before 2026-10-10 still sit directly in
`results/` under timestamped file names.

`make full` runs both tests into one folder, then runs
`bash scripts/process-results.bash results/TIMESTAMP`.  That writes `acc-speed.csv` and
`acc-speed.filtered.csv` next to the JSON, and `make full` then plots the filtered CSV into
`acc-speed.png` in the same folder.  To plot one by hand, run
`python3 scripts/plot-results.py -o OUT.png results/TIMESTAMP/acc-speed.filtered.csv`, or leave
out `-o` to show the plot in a window.

## Build constraints

- The compiler must be g++ with `-std=c++26`.  The Makefile states that clang++ is not
  supported.
- Required libraries are Google Benchmark, nlohmann-json, and oneTBB.  The Makefile aborts if
  any program in `REQUIRED_BINS` (including `jq`, `sponge`, and `numfmt`) is missing.
- Never enable `-ffinite-math-only`, `-ffast-math`, or `-Ofast`.  The other floating-point
  flags left commented in the Makefile trade accuracy for speed and are off on purpose.

## Architecture

`include/ecef_to_geodetic-funcs.hpp` holds every algorithm under test.  Each one lives in its
own namespace with a `double`-only `ecef_to_geodetic(x, y, z, lat_rad, lon_rad, ht)` and a
`func_info` object.  An `_xN` suffix is the iteration count, and `_customht` marks an
algorithm that uses its own height formula instead of `ell.get_ht`.

The `func_info_t` metadata includes a code-size figure computed from `__LINE__` markers placed
around the function body.  The constructor adds the lines of whichever prologue macro the
function uses, `COMMON_FIRST_DECLS` or `COMMON_FIRST_DECLS_CHECKED`, and `lines_extra` counts
shared helpers.  Keep those markers and constants accurate when editing a function, because
the line counts are reported as results.

A hand-counted constant runs from a helper's `template` line through its closing brace,
without its doc comments.  A constant that covers several helpers is written as a sum with one
term per helper, such as `6 + 5 + 5 + 10 + 6 + 9`.

To add an algorithm, write a new namespace in that file and add a matching entry to
`include/map_func_name_to_func_info.hpp`.  A comment in that header gives the `grep | awk`
command that regenerates the map from the namespaces.  Both test programs discover
algorithms only through that map.

The accuracy test measures round-trip error.  It converts each input point to ECEF exactly,
runs the algorithm under test, converts the result back exactly with
`include/geodetic_to_ecef.hpp`, and reports the Euclidean distance between the two ECEF
points.  Every algorithm returns the same longitude, so inputs vary only latitude and
height.

`include/ecef_to_geodetic.hpp` is separate from the test harness.  It is a standalone,
templated Olson 1996 implementation built on the library types in `include/` (`angle`,
`ECEF`, `Geodetic`, and the ellipsoids), and it does not depend on the test code.

Two headers are generated.  `generate/generate.bash` runs `calc` scripts to write
`include/aux-lat-conv.hpp` and `include/utm-ups-const.hpp`, so edit the `.cal` sources rather
than those headers.

The Python scripts in `scripts/` generate and plot input data.  `Nd-arange.py` piped into
`polar-to-cartesian.py` produces the ECEF grids, and `plot-points.py` draws them (see the
`plot-ecef` and `plot-geod` targets).  `ellipsoid.py` holds the scalar ellipsoid math that the
converter and plotting scripts share.  `process-results.bash`, `filter-benchmark-results.jq`,
and `plot-results.py` turn the test results into CSV files and a plot.  `olson_1996/` is the
original C version with its own Makefile.

## Style

All prose, including code comments, doc blocks, commit messages, and documentation files,
follows `COMMENT-STYLE.md`.  Read it before writing any of those.  Headers use Doxygen with
`///` briefs above `/** */` detail blocks, and each file carries SPDX license headers.

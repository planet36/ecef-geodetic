# ecef-geodetic

A test suite of ECEF-to-geodetic coordinate conversion functions

## Requirements

### C++ libraries

- [Google Benchmark](https://github.com/google/benchmark)
- [fmt](https://github.com/fmtlib/fmt)
- [nlohmann-json](https://github.com/nlohmann/json)
- [oneTBB](https://github.com/oneapi-src/oneTBB)

### Python modules

- [gmpy2](https://github.com/aleaxit/gmpy)
- [matplotlib](https://matplotlib.org/)
- [more-itertools](https://github.com/more-itertools/more-itertools)
- [numpy](https://numpy.org/)

### Programs

See [Makefile](Makefile) for detailed list of programs that are required to run.

## Results

Only algorithms with a mean distance error less than 10nm were included in the figure below.  All iterative algorithms did 2 iterations.  "(c.h.)" means the algorithm used a "custom height" formula to calculate ellipsoid height instead of the standard formula.

![Scatter plot of accuracy vs speed](results/acc-speed.20240112T204815.png)
*Scatter plot of ***accuracy***, measured by mean distance error (nm), versus ***speed***, measured in millions of conversions per second*

Results closer to the upper-left corner of the figure are better.

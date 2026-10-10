# ecef-geodetic

A test suite of ECEF-to-geodetic coordinate conversion functions

## Usage

To create the input data files
> `make -j $(nproc) input`

They are text files with ECEF and Geodetic coordinates.

To build the accuracy and speed tests
> `make -j $(nproc)`

To run the accuracy and speed tests
> `make full`
* Takes about 1 minute to finish
* Runs the two tests one after the other, even under `make -j`

The output of the accuracy and speed tests are json files which are put in a timestamped folder
within the `results` folder.
`make full` then combines them into a CSV file and a filtered CSV file in that same folder, and
plots the filtered one into a PNG file there.

See [Makefile](Makefile) for all possible targets.

## Requirements

### C++ libraries

- [Google Benchmark](https://github.com/google/benchmark)
- [nlohmann-json](https://github.com/nlohmann/json)
- [oneTBB](https://github.com/oneapi-src/oneTBB)

### Python modules

- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)

### Programs

See [Makefile](Makefile) for detailed list of programs that are required to run.

## Results

Only algorithms with a mean distance error less than 10nm were included in the figure below.  All iterative algorithms did 2 iterations.  "(c.h.)" means the algorithm used a "custom height" formula to calculate ellipsoid height instead of the standard formula.

![Scatter plot of accuracy vs speed](results/20261010T150508/acc-speed.png)
*Scatter plot of ***accuracy***, measured by mean distance error (nm), versus ***speed***, measured in millions of conversions per second*

Results closer to the upper-left corner of the figure are better.

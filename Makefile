# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

export LC_ALL := C

# https://how.wtf/check-if-a-program-exists-from-a-makefile.html
REQUIRED_BINS := \
awk \
bash \
g++ \
grep \
join \
jq \
mkdir \
python3 \
readelf \
rm \
sed \
sort \
sponge \
tr

$(foreach bin,$(REQUIRED_BINS),\
    $(if $(shell command -v $(bin) 2> /dev/null),,$(error Please install `$(bin)`)))

# clang++ not supported
CXX = g++

CPPFLAGS = -MMD -MP
CPPFLAGS += -Iinclude

CXXFLAGS = -pipe -Wall -Wextra -Wpedantic -Wfatal-errors
CXXFLAGS += -std=c++23
# -frecord-gcc-switches is used by readelf
CXXFLAGS += -frecord-gcc-switches
CXXFLAGS += -O3 -flto=auto -march=native -fno-math-errno
# https://gcc.gnu.org/wiki/FloatingPointMath
#CXXFLAGS += -freciprocal-math # slightly decreased accuracy, slightly increased speed
#CXXFLAGS += -fno-signed-zeros # slightly decreased accuracy, slightly increased speed
#CXXFLAGS += -fno-trapping-math # slightly decreased accuracy, same speed
# NOTE: -fassociative-math requires -fno-signed-zeros and -fno-trapping-math
#CXXFLAGS += -fassociative-math -fno-signed-zeros -fno-trapping-math # almost same accuracy, almost same speed

# Do not use -ffinite-math-only (enabled with -ffast-math (enabled with -Ofast))

#LDFLAGS +=

LDLIBS += -lbenchmark
LDLIBS += -lfmt
LDLIBS += -ltbb

ALL_INFILES_GEOD := \
geod.2d.region-0.txt \
geod.2d.region-1.txt \
geod.2d.region-2.txt \
geod.2d.region-3.txt \
geod.2d.region-4.txt \
geod.2d.region-all.txt \
geod.2d.neg-ht-1.txt \
geod.2d.neg-ht-2.txt \

ALL_INFILES_ECEF := \
ecef.2d.region-0.txt \
ecef.2d.region-1.txt \
ecef.2d.region-2.txt \
ecef.2d.region-3.txt \
ecef.2d.region-4.txt \
ecef.2d.region-all.txt \
ecef.2d.speed.txt \

# Use N-1 threads in the speed test
export NUM_THREADS := $(shell nproc --ignore 1)

# Should be an odd number for simpler median
NUM_SPEED_TESTS := 11

# Used by plot
DPI := 180

DATETIME := $(shell date -u +'%Y%m%dT%H%M%S')

OUTPUT_DIR := results

SRC_ACC := test-ecef_to_geodetic-acc.cpp
#BIN_ACC = $(addsuffix .out, $(basename $(SRC_ACC)))
BIN_ACC = $(basename $(SRC_ACC))

SRC_SPEED := test-ecef_to_geodetic-speed.cpp
#BIN_SPEED = $(addsuffix .out, $(basename $(SRC_SPEED)))
BIN_SPEED = $(basename $(SRC_SPEED))

SRCS := $(SRC_ACC) $(SRC_SPEED)
DEPS := $(SRC_ACC:.cpp=.d) $(SRC_SPEED:.cpp=.d)
#OBJS := $(SRC_ACC:.cpp=.o) $(SRC_SPEED:.cpp=.o)
BINS := $(BIN_ACC) $(BIN_SPEED)

all: $(BINS) input | $(OUTPUT_DIR)

# The built-in recipe for the implicit rule uses $^ instead of $<
%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)
	@# Extract compile options
	readelf -p .GCC.command.line $@ | grep -F 'GNU GIMPLE' | sed -E -e 's/^\s*\[\s*[0-9]+\]\s*//' | tr -d '\n' > $@.opts

input: $(ALL_INFILES_ECEF) $(ALL_INFILES_GEOD)

# ECEF points
# NOTE: They can be generated 2 ways:
# 1) vary W (meters) and Z (meters)
# 2) vary r (meters) and theta (degrees), and convert from polar to cartesian
#    r is distance from center of earth
#    theta is geocentric latitude (not geodetic)
#    polar-to-cartesian.py is good enough for this case, even though it's input is geocentric latitude.

ecef.2d.region-0.txt:
	python3 Nd-arange.py 0 50_000 1_000 0 90 30 | python3 polar-to-cartesian.py > $@

ecef.2d.region-1.txt:
	python3 Nd-arange.py 0 7_000_000 100_000 0 90 30 | python3 polar-to-cartesian.py > $@

ecef.2d.region-2.txt:
	python3 Nd-arange.py 6_300_000 6_500_000 1_000 0 90 30 | python3 polar-to-cartesian.py > $@

ecef.2d.region-3.txt:
	python3 Nd-arange.py 6_350_000 6_400_000 100 0 90 30 | python3 polar-to-cartesian.py > $@

ecef.2d.region-4.txt:
	python3 Nd-arange.py 0 100_000_000 1_000_000 0 90 30 | python3 polar-to-cartesian.py > $@

ecef.2d.region-all.txt: ecef.2d.region-0.txt ecef.2d.region-1.txt ecef.2d.region-2.txt ecef.2d.region-3.txt ecef.2d.region-4.txt
	LC_ALL=C sort -u -- $^ > $@

ecef.2d.speed.txt: create-speed-points.py
	python3 create-speed-points.py > $@

# Geodetic points
# vary geodetic latitude (degrees) and ellipsoid height (meters)

geod.2d.region-0.txt:
	python3 Nd-arange.py -90 90 0.0001 0 0 1 > $@

geod.2d.region-1.txt:
	python3 Nd-arange.py -90 90 0.001 -100 1000 100 > $@

geod.2d.region-2.txt:
	python3 Nd-arange.py -90 90 0.01 -10_000 100_000 1_000 > $@

geod.2d.region-3.txt:
	python3 Nd-arange.py -90 90 0.1 -1_000_000 10_000_000 10_000 > $@

geod.2d.region-4.txt:
	python3 Nd-arange.py -90 90 1 -5_000_000 500_000_000 100_000 > $@

geod.2d.region-all.txt: geod.2d.region-0.txt geod.2d.region-1.txt geod.2d.region-2.txt geod.2d.region-3.txt geod.2d.region-4.txt
	LC_ALL=C sort -u -- $^ > $@

geod.2d.neg-ht-1.txt:
	python3 Nd-arange.py -90 90 5 -6_383_000 0 1_000 > $@

geod.2d.neg-ht-2.txt:
	python3 Nd-arange.py -90 90 1 -6_383_000 0 10_000 > $@

# https://www.gnu.org/software/make/manual/html_node/Double_002dColon.html

plot-ecef:: ecef.2d.region-all.txt
	for F in $^; do python3 plot-points.py -v --ell --evo --lim --km --dpi=$(DPI) < $$F; done

plot-ecef:: ecef.2d.speed.txt
	for F in $^; do python3 plot-points.py -v --ell --evo       --km --dpi=$(DPI) < $$F; done

plot-geod:: geod.2d.region-0.txt geod.2d.region-1.txt geod.2d.region-2.txt geod.2d.region-3.txt geod.2d.region-4.txt
	for F in $^; do python3 plot-points.py -v -g --ell --evo --lim --km --dpi=$(DPI) < $$F; done

plot-geod:: geod.2d.neg-ht-1.txt geod.2d.neg-ht-2.txt
	for F in $^; do python3 plot-points.py -v -g --ell --evo       --km --dpi=$(DPI) < $$F; done

acc: $(BIN_ACC) input | $(OUTPUT_DIR)
	./$< -v -t -g -a < geod.2d.region-all.txt > $(OUTPUT_DIR)/$@.$(DATETIME).json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' < $(OUTPUT_DIR)/$@.$(DATETIME).json | sponge $(OUTPUT_DIR)/$@.$(DATETIME).json

acc1: $(BIN_ACC) input | $(OUTPUT_DIR)
	@# NOTE: Only run this test with a few input points
	./$< -v -t -1 < ecef.2d.speed.txt > $(OUTPUT_DIR)/$@.$(DATETIME).json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' < $(OUTPUT_DIR)/$@.$(DATETIME).json | sponge $(OUTPUT_DIR)/$@.$(DATETIME).json

speed: $(BIN_SPEED) input | $(OUTPUT_DIR)
	@# NOTE: The input data format must be ECEF, not Geodetic
	./$< \
		--benchmark_enable_random_interleaving=true \
		--benchmark_repetitions=$(NUM_SPEED_TESTS) \
		--benchmark_report_aggregates_only=true \
		--benchmark_out_format=json \
		--benchmark_out=$(OUTPUT_DIR)/$@.$(DATETIME).json \
		< ecef.2d.speed.txt

	@# Preserve the given order because --benchmark_enable_random_interleaving=true shuffles the order of the tests.
	jq '.benchmarks |= sort_by(.family_index)' < $(OUTPUT_DIR)/$@.$(DATETIME).json | sponge $(OUTPUT_DIR)/$@.$(DATETIME).json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' < $(OUTPUT_DIR)/$@.$(DATETIME).json | sponge $(OUTPUT_DIR)/$@.$(DATETIME).json

$(OUTPUT_DIR):
	mkdir --verbose --parents -- $@

clean:
	@$(RM) --verbose -- $(DEPS) $(BINS) *.opts

clean-input:
	@$(RM) --verbose -- \
		$(ALL_INFILES_ECEF) $(ALL_INFILES_GEOD)

clean-output:
	git clean -i $(OUTPUT_DIR)

clean-all: clean clean-input clean-output

lint:
	-clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: all input plot-ecef plot-geod acc acc1 speed clean clean-input clean-output clean-all lint

-include $(DEPS)

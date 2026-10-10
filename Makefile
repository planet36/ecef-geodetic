# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

export LC_ALL = C

# https://how.wtf/check-if-a-program-exists-from-a-makefile.html
REQUIRED_BINS := \
awk \
bash \
basename \
cat \
dirname \
g++ \
grep \
join \
jq \
mkdir \
numfmt \
python3 \
readelf \
rm \
rmdir \
sed \
sort \
sponge \
tr

$(foreach bin,$(REQUIRED_BINS),\
    $(if $(shell command -v $(bin) 2> /dev/null),,$(error Please install `$(bin)`)))

# clang++ not supported
CXX = g++

CPPFLAGS = -MMD -MP
CPPFLAGS += -I include

CXXFLAGS = -std=c++26
CXXFLAGS += -pipe -Wall -Wextra -Wpedantic -Wfatal-errors
CXXFLAGS += -O3 -flto=auto -march=native
CXXFLAGS += -fno-math-errno
#CXXFLAGS += -march=raptorlake
# -frecord-gcc-switches is used by readelf
CXXFLAGS += -frecord-gcc-switches
# https://gcc.gnu.org/wiki/FloatingPointMath
#CXXFLAGS += -freciprocal-math # slightly decreased accuracy, slightly increased speed
#CXXFLAGS += -fno-signed-zeros # slightly decreased accuracy, slightly increased speed
#CXXFLAGS += -fno-trapping-math # slightly decreased accuracy, same speed
# NOTE: -fassociative-math requires -fno-signed-zeros and -fno-trapping-math
#CXXFLAGS += -fassociative-math -fno-signed-zeros -fno-trapping-math # almost same accuracy, almost same speed

# Do not use -ffinite-math-only (enabled with -ffast-math (enabled with -Ofast))

#LDFLAGS =

LDLIBS = -lbenchmark -ltbb

INPUT_DIR = input

ALL_INFILES_GEOD := \
$(INPUT_DIR)/geod.2d.region-0.txt \
$(INPUT_DIR)/geod.2d.region-1.txt \
$(INPUT_DIR)/geod.2d.region-2.txt \
$(INPUT_DIR)/geod.2d.region-3.txt \
$(INPUT_DIR)/geod.2d.region-4.txt \
$(INPUT_DIR)/geod.2d.region-all.txt \
$(INPUT_DIR)/geod.2d.neg-ht-1.txt \
$(INPUT_DIR)/geod.2d.neg-ht-2.txt \

ALL_INFILES_ECEF := \
$(INPUT_DIR)/ecef.2d.region-0.txt \
$(INPUT_DIR)/ecef.2d.region-1.txt \
$(INPUT_DIR)/ecef.2d.region-2.txt \
$(INPUT_DIR)/ecef.2d.region-3.txt \
$(INPUT_DIR)/ecef.2d.region-4.txt \
$(INPUT_DIR)/ecef.2d.region-all.txt \
$(INPUT_DIR)/ecef.2d.speed.txt \

# Use N-1 threads in the speed test
export NUM_THREADS := $(shell nproc --ignore 1)

# Should be an odd number for simpler median
BENCHMARK_REPS = 5

# Used by plot
DPI = 180

DATETIME := $(shell date -u +'%Y%m%dT%H%M%S')

OUTPUT_DIR = results/$(DATETIME)

SCRIPTS_DIR = scripts

SRC_ACC = test-ecef_to_geodetic-acc.cpp
#BIN_ACC = $(addsuffix .out, $(basename $(SRC_ACC)))
BIN_ACC = $(basename $(SRC_ACC))

SRC_SPEED = test-ecef_to_geodetic-speed.cpp
#BIN_SPEED = $(addsuffix .out, $(basename $(SRC_SPEED)))
BIN_SPEED = $(basename $(SRC_SPEED))

SRCS = $(SRC_ACC) $(SRC_SPEED)
DEPS = $(SRC_ACC:.cpp=.d) $(SRC_SPEED:.cpp=.d)
#OBJS = $(SRC_ACC:.cpp=.o) $(SRC_SPEED:.cpp=.o)
BINS = $(BIN_ACC) $(BIN_SPEED)

all: $(BINS) input

# The built-in recipe for the implicit rule uses $^ instead of $<
%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)
	@# Extract compile options
	readelf -p .GCC.command.line $@ | grep -F 'GNU GIMPLE' | \
		sed -E -e 's/^\s*\[\s*[0-9]+\]\s*//' | tr -d '\n' > $@.opts

input: $(ALL_INFILES_ECEF) $(ALL_INFILES_GEOD)

# The trailing /. keeps the directory distinct from the phony target of the same name.
$(ALL_INFILES_ECEF) $(ALL_INFILES_GEOD): | $(INPUT_DIR)/.

# ECEF points
# NOTE: They can be generated 2 ways:
# 1) vary W (meters) and Z (meters)
# 2) vary r (meters) and theta (degrees), and convert from polar to cartesian
#    r is distance from center of earth
#    theta is geocentric latitude (not geodetic)
#    polar-to-cartesian.py is good enough for this case, even though it's input is geocentric latitude.

$(INPUT_DIR)/ecef.2d.region-0.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py 0 50_000 1_000 0 90 30 | python3 $(SCRIPTS_DIR)/polar-to-cartesian.py > $@

$(INPUT_DIR)/ecef.2d.region-1.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py 0 7_000_000 100_000 0 90 30 | python3 $(SCRIPTS_DIR)/polar-to-cartesian.py > $@

$(INPUT_DIR)/ecef.2d.region-2.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py 6_300_000 6_500_000 1_000 0 90 30 | python3 $(SCRIPTS_DIR)/polar-to-cartesian.py > $@

$(INPUT_DIR)/ecef.2d.region-3.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py 6_350_000 6_400_000 100 0 90 30 | python3 $(SCRIPTS_DIR)/polar-to-cartesian.py > $@

$(INPUT_DIR)/ecef.2d.region-4.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py 0 100_000_000 1_000_000 0 90 30 | python3 $(SCRIPTS_DIR)/polar-to-cartesian.py > $@

$(INPUT_DIR)/ecef.2d.region-all.txt: $(INPUT_DIR)/ecef.2d.region-0.txt $(INPUT_DIR)/ecef.2d.region-1.txt $(INPUT_DIR)/ecef.2d.region-2.txt $(INPUT_DIR)/ecef.2d.region-3.txt $(INPUT_DIR)/ecef.2d.region-4.txt
	LC_ALL=C sort -u -- $^ > $@

$(INPUT_DIR)/ecef.2d.speed.txt: $(SCRIPTS_DIR)/create-speed-points.py
	python3 $(SCRIPTS_DIR)/create-speed-points.py > $@

# Geodetic points
# vary geodetic latitude (degrees) and ellipsoid height (meters)

$(INPUT_DIR)/geod.2d.region-0.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 0.0001 0 0 1 > $@

$(INPUT_DIR)/geod.2d.region-1.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 0.001 -100 1000 100 > $@

$(INPUT_DIR)/geod.2d.region-2.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 0.01 -10_000 100_000 1_000 > $@

$(INPUT_DIR)/geod.2d.region-3.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 0.1 -1_000_000 10_000_000 10_000 > $@

$(INPUT_DIR)/geod.2d.region-4.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 1 -5_000_000 500_000_000 100_000 > $@

$(INPUT_DIR)/geod.2d.region-all.txt: $(INPUT_DIR)/geod.2d.region-0.txt $(INPUT_DIR)/geod.2d.region-1.txt $(INPUT_DIR)/geod.2d.region-2.txt $(INPUT_DIR)/geod.2d.region-3.txt $(INPUT_DIR)/geod.2d.region-4.txt
	LC_ALL=C sort -u -- $^ > $@

$(INPUT_DIR)/geod.2d.neg-ht-1.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 5 -6_383_000 0 1_000 > $@

$(INPUT_DIR)/geod.2d.neg-ht-2.txt:
	python3 $(SCRIPTS_DIR)/Nd-arange.py -90 90 1 -6_383_000 0 10_000 > $@

# https://www.gnu.org/software/make/manual/html_node/Double_002dColon.html

plot-ecef:: $(INPUT_DIR)/ecef.2d.region-all.txt
	for F in $^; do python3 $(SCRIPTS_DIR)/plot-points.py -v --ell --evo --lim --km --dpi=$(DPI) < $$F; done

plot-ecef:: $(INPUT_DIR)/ecef.2d.speed.txt
	for F in $^; do python3 $(SCRIPTS_DIR)/plot-points.py -v --ell --evo       --km --dpi=$(DPI) < $$F; done

plot-geod:: $(INPUT_DIR)/geod.2d.region-0.txt $(INPUT_DIR)/geod.2d.region-1.txt $(INPUT_DIR)/geod.2d.region-2.txt $(INPUT_DIR)/geod.2d.region-3.txt $(INPUT_DIR)/geod.2d.region-4.txt
	for F in $^; do python3 $(SCRIPTS_DIR)/plot-points.py -v -g --ell --evo --lim --km --dpi=$(DPI) < $$F; done

plot-geod:: $(INPUT_DIR)/geod.2d.neg-ht-1.txt $(INPUT_DIR)/geod.2d.neg-ht-2.txt
	for F in $^; do python3 $(SCRIPTS_DIR)/plot-points.py -v -g --ell --evo       --km --dpi=$(DPI) < $$F; done

acc: $(BIN_ACC) input | $(OUTPUT_DIR)
	./$< -v -t -g -a < $(INPUT_DIR)/geod.2d.region-all.txt > $(OUTPUT_DIR)/$@.json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' \
		< $(OUTPUT_DIR)/$@.json | sponge $(OUTPUT_DIR)/$@.json

acc1: $(BIN_ACC) input | $(OUTPUT_DIR)
	@# NOTE: Only run this test with a few input points
	./$< -v -t -1 < $(INPUT_DIR)/ecef.2d.speed.txt > $(OUTPUT_DIR)/$@.json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' \
		< $(OUTPUT_DIR)/$@.json | sponge $(OUTPUT_DIR)/$@.json

speed: $(BIN_SPEED) input | $(OUTPUT_DIR)
	@# NOTE: The input data format must be ECEF, not Geodetic
	./$< \
		--benchmark_enable_random_interleaving=true \
		--benchmark_repetitions=$(BENCHMARK_REPS) \
		--benchmark_report_aggregates_only=true \
		--benchmark_out_format=json \
		--benchmark_out=$(OUTPUT_DIR)/$@.json \
		< $(INPUT_DIR)/ecef.2d.speed.txt

	@# Preserve the given order because --benchmark_enable_random_interleaving=true shuffles the order of the tests.
	jq '.benchmarks |= sort_by(.family_index)' \
		< $(OUTPUT_DIR)/$@.json | sponge $(OUTPUT_DIR)/$@.json

	@# Insert compile options
	jq --rawfile compile_opts $<.opts '. + {compile_opts: $$compile_opts}' \
		< $(OUTPUT_DIR)/$@.json | sponge $(OUTPUT_DIR)/$@.json

# Run both tests from the same make invocation, so they share $(OUTPUT_DIR), then combine and
# plot them.
# .WAIT keeps the tests from running in parallel under make -j.
full: $(BIN_ACC) $(BIN_SPEED) input .WAIT acc .WAIT speed
	@bash $(SCRIPTS_DIR)/process-results.bash $(OUTPUT_DIR)
	@python3 $(SCRIPTS_DIR)/plot-results.py -o $(OUTPUT_DIR)/acc-speed.png \
		$(OUTPUT_DIR)/acc-speed.filtered.csv

$(OUTPUT_DIR) $(INPUT_DIR)/.:
	@mkdir --verbose --parents -- $@

clean:
	@$(RM) --verbose -- $(DEPS) $(BINS) *.opts

clean-input:
	@$(RM) --verbose -- \
		$(ALL_INFILES_ECEF) $(ALL_INFILES_GEOD)
	@if [ -d $(INPUT_DIR) ]; then rmdir --verbose --ignore-fail-on-non-empty -- $(INPUT_DIR); fi

clean-all: clean clean-input

lint:
	-clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CXXFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: all input plot-ecef plot-geod acc acc1 speed full clean clean-input clean-all lint

# https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#index-removing-targets-on-failure
.DELETE_ON_ERROR:

-include $(DEPS)

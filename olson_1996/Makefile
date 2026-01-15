# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

export LC_ALL = C

SRCS = $(wildcard *.c)
DEPS = $(SRCS:.c=.d)
OBJS = $(SRCS:.c=.o)

BIN = ecef_to_geodetic

CPPFLAGS = -MMD -MP
CPPFLAGS += -D_GNU_SOURCE

CFLAGS = -std=c23
CFLAGS += -pipe -Wall -Wextra -Wpedantic -Wfatal-errors
CFLAGS += -O3 -flto=auto -march=native -fno-math-errno
CFLAGS += -Wno-old-style-definition

#LDFLAGS =

LDLIBS = -lm

$(BIN): $(OBJS)
	$(CC) $^ -o $@ $(LDLIBS)

clean:
	@$(RM) --verbose -- $(DEPS) $(OBJS) $(BIN)

lint:
	-clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: clean lint

-include $(DEPS)

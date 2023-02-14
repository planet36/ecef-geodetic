# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

SRCS = $(wildcard *.c)
DEPS = $(SRCS:.c=.d)
OBJS = $(SRCS:.c=.o)

BIN = ecef_to_geodetic

CPPFLAGS += -MMD -MP
CFLAGS += -Wno-old-style-definition
#LDFLAGS +=
LDLIBS += -lm

$(BIN): $(OBJS)
	$(CC) $^ -o $@ $(LDLIBS)

clean:
	@$(RM) --verbose -- $(DEPS) $(OBJS) $(BIN) *.so test *~ core core.[[:digit:]]*

lint:
	-clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: clean lint

-include $(DEPS)

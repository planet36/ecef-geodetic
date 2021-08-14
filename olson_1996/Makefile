# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

SRCS = $(wildcard *.c)
DEPS = $(SRCS:.c=.d)
OBJS = $(SRCS:.c=.o)

CPPFLAGS += -MMD -MP
LDFLAGS += -lm

ecef_to_geodetic: $(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)

clean:
	@$(RM) --verbose -- $(DEPS) $(OBJS) ecef_to_geodetic *.so test *~ core core.[[:digit:]]*

lint:
	clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CFLAGS) $(LDFLAGS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: clean lint

-include $(DEPS)

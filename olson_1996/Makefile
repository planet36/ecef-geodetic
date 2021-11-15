# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

SRCS = $(wildcard *.c)
DEPS = $(SRCS:.c=.d)
OBJS = $(SRCS:.c=.o)

CPPFLAGS += -MMD -MP
CFLAGS += -Wno-old-style-definition
LDLIBS += -lm

ecef_to_geodetic: $(OBJS)
	$(CC) $^ -o $@ $(LDLIBS)

clean:
	@$(RM) --verbose -- $(DEPS) $(OBJS) ecef_to_geodetic *.so test *~ core core.[[:digit:]]*

lint:
	clang-tidy --quiet $(SRCS) -- $(CPPFLAGS) $(CFLAGS) $(LDLIBS)

# https://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: clean lint

-include $(DEPS)

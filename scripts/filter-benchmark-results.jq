# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# With repetitions, use only the median.  Without them, every entry is a single run.
(any(.benchmarks[]; .aggregate_name == "median")) as $has_median |
.benchmarks[] |
select(if $has_median then .aggregate_name == "median" else true end) |
"\(.name | sub("/threads:[0-9]+(_median)?$"; "")),\(.cpu_time)"

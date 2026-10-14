# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

.func_names |
keys[] as $k |
"\(.[$k].info.display_name),\(.[$k].acc.mean_dist_err*1e9),\(.[$k].acc.max_dist_err*1e9)"

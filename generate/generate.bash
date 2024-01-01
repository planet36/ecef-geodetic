#!/usr/bin/bash
# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# shellcheck disable=SC2034

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname -- "${BASH_SOURCE[0]}")"

set -e

cd "$SCRIPT_DIR"

OUT_DIR="$(git rev-parse --show-toplevel)/include"

calc -d -m 4 -f generate-series-approx-coeff-aux-lat.cal > "$OUT_DIR"/aux-lat-conv.hpp
calc -d -m 4 -f generate-series-approx-coeff-tm.cal      > "$OUT_DIR"/utm-ups-const.hpp

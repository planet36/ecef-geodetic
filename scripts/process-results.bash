# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname -- "${BASH_SOURCE[0]}")"

export LC_ALL=C

function print_usage
{
    cat <<EOT 1>&2
Usage: bash $SCRIPT_NAME RESULTS_DIR

"RESULTS_DIR" is a folder that holds the products of the "acc" and "speed" tests, acc.json and
speed.json.  The CSV files are written to the same folder.
EOT
}

if (($# < 1))
then
    print_usage
    exit 1
fi

declare -r RESULTS_DIR="$1"

if [[ ! -d "$RESULTS_DIR" ]]
then
    printf 'Error: %q is not a folder\n' "$RESULTS_DIR" 1>&2
    print_usage
    exit 1
fi

declare -r INFILE_ACC="${RESULTS_DIR}/acc.json"
declare -r INFILE_SPEED="${RESULTS_DIR}/speed.json"

if [[ ! -f "$INFILE_ACC" ]]
then
    printf 'Error: %q does not exist\n' "$INFILE_ACC" 1>&2
    print_usage
    exit 1
fi

if [[ ! -f "$INFILE_SPEED" ]]
then
    printf 'Error: %q does not exist\n' "$INFILE_SPEED" 1>&2
    print_usage
    exit 1
fi

declare -r OUTFILE="${RESULTS_DIR}/acc-speed.csv"
declare -r OUTFILE_FILTERED="${RESULTS_DIR}/acc-speed.filtered.csv"

printf 'Name,Mean dist. error (nm),Max dist. error (nm),Time/call (ns)\n' > "$OUTFILE"

join -t ',' \
    <(jq -r --from-file "${SCRIPT_DIR}/filter-accuracy-results.jq" "$INFILE_ACC" | sort) \
    <(jq -r --from-file "${SCRIPT_DIR}/filter-benchmark-results.jq" "$INFILE_SPEED" | sort) >> \
    "$OUTFILE" || exit

# Select algorithms with good accuracy (mean dist err < 10 nm).
# Change precision of output to 3 decimal places.
awk -F ',' '$2 < 10 || NR == 1 {print $0}' "$OUTFILE" \
    | numfmt --header=1 --delimiter=, --field=2- --format='%0.3f' \
    > "$OUTFILE_FILTERED" || exit

printf 'Created files:\n%q\n%q\n' "$OUTFILE" "$OUTFILE_FILTERED"

printf 'To show the plot in a window, run:\npython3 %q %q\n' \
    "${SCRIPT_DIR}/plot-results.py" \
    "$OUTFILE_FILTERED"

# Use datamash to get stats of the accurate algorithms.
# Example:
# datamash --header-in --field-separator=',' q1 4 mean 4 median 4 q3 4 iqr 4 < "$OUTFILE_FILTERED"

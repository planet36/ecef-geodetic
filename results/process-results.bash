# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

SCRIPT_NAME="$(basename -- "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname -- "${BASH_SOURCE[0]}")"

function print_usage
{
    cat <<EOT 1>&2
Usage: bash $SCRIPT_NAME INFILE_ACC INFILE_SPEED

"INFILE_ACC" and "INFILE_SPEED" are the respective products of the "acc" and "speed" tests.
EOT
}

if (($# < 2))
then
    print_usage
    exit 1
fi

declare -r INFILE_ACC="$1"
declare -r INFILE_SPEED="$2"

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

DATETIME="$(date -u +'%Y%m%dT%H%M%S')"
readonly DATETIME

declare -r OUTFILE="${SCRIPT_DIR}/acc-speed.${DATETIME}.csv"
declare -r OUTFILE_FILTERED="${SCRIPT_DIR}/acc-speed.${DATETIME}.filtered.csv"

printf '%q\n' "$INFILE_ACC" > "$OUTFILE"
printf '%q\n' "$INFILE_SPEED" >> "$OUTFILE"
printf 'Name,Mean dist. error (nm),Max dist. error (nm),M conversions/sec\n' >> "$OUTFILE"

if grep -q -E '/threads:[0-9]+_median\>' "$INFILE_SPEED"
then
    # Speed test was run with more than 1 repetition.
    join -t ',' \
        <(jq -r '.func_names | keys[] as $k | "\(.[$k].info.display_name),\(.[$k].acc.mean_dist_err*1e9),\(.[$k].acc.max_dist_err*1e9)"' "$INFILE_ACC" | sort) \
        <(jq -r '.benchmarks[] | select(.name | endswith("_median")) | "\(.name),\(.items_per_second/1e6)"' "$INFILE_SPEED" | sed -E -e 's/\/threads:[0-9]+_median//' | sort) >> \
        "$OUTFILE" || exit
else
    # Speed test was run with 1 repetition.
    join -t ',' \
        <(jq -r '.func_names | keys[] as $k | "\(.[$k].info.display_name),\(.[$k].acc.mean_dist_err*1e9),\(.[$k].acc.max_dist_err*1e9)"' "$INFILE_ACC" | sort) \
        <(jq -r '.benchmarks[] | "\(.name),\(.items_per_second/1e6)"' "$INFILE_SPEED" | sed -E -e 's/\/threads:[0-9]+//' | sort) >> \
        "$OUTFILE" || exit
fi

# Select algorithms with good accuracy (mean dist err < 10 nm).
# Change precision of output to 3 decimal places.
awk --csv '$2 < 10 || NR <= 3 {print $0}' "$OUTFILE" \
    | numfmt --header=3 --delimiter=, --field=2- --format='%0.3f' \
    > "$OUTFILE_FILTERED" || exit

printf 'Created files:\n%q\n%q\n' "$OUTFILE" "$OUTFILE_FILTERED"

# Use datamash to get stats of the accurate algorithms.
# Example:
# datamash --header-in --field-separator=',' q1 4 mean 4 median 4 q3 4 iqr 4 < "$OUTFILE_FILTERED"

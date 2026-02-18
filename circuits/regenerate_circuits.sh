#!/bin/bash
# Regenerate Bristol Fashion circuits for the given Y_TYPE_BITS.
# Usage: ./regenerate_circuits.sh [Y_TYPE_BITS] [OUT_DIR]
#   Default: 64

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

Y=${1:-64}
OUT_DIR=${2:-.}

mkdir -p "$OUT_DIR" "$OUT_DIR/dummy_check" "$OUT_DIR/replace_if_dummy"

echo "Regenerating circuits for Y_TYPE_BITS=$Y ..."

python3 xy_if_xs_equal.py "$Y" > "$OUT_DIR/xy_if_xs_equal.txt"
python3 compare_swap.py "$Y" > "$OUT_DIR/compare_swap.txt"
python3 add_mul_64.py > "$OUT_DIR/add_mul_64.txt"
python3 add_less_64.py > "$OUT_DIR/add_less_64.txt"
python3 add_eq_64.py > "$OUT_DIR/add_eq_64.txt"

# dummy_check and replace_if_dummy only operate on x_type, no Y_TYPE_BITS dependency
for logN in $(seq 6 31); do
    python3 dummy_check.py "$logN" > "$OUT_DIR/dummy_check/$logN.txt"
    python3 replace_if_dummy.py "$logN" > "$OUT_DIR/replace_if_dummy/$logN.txt"
done

# cht_lookup does not depend on Y_TYPE_BITS
python3 cht_lookup.py > "$OUT_DIR/cht_lookup.txt"

echo "Done."

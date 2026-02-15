#!/bin/bash
# Regenerate Bristol Fashion circuits for the given Y_TYPE_BITS.
# Usage: ./regenerate_circuits.sh [Y_TYPE_BITS]
#   Default: 64

set -e
cd "$(dirname "$0")"

Y=${1:-64}

echo "Regenerating circuits for Y_TYPE_BITS=$Y ..."

python3 xy_if_xs_equal.py "$Y" > xy_if_xs_equal.txt
python3 compare_swap.py "$Y" > compare_swap.txt

# dummy_check and replace_if_dummy only operate on x_type, no Y_TYPE_BITS dependency
for logN in $(seq 6 31); do
    python3 dummy_check.py "$logN" > "dummy_check/$logN.txt"
    python3 replace_if_dummy.py "$logN" > "replace_if_dummy/$logN.txt"
done

# cht_lookup does not depend on Y_TYPE_BITS
python3 cht_lookup.py > cht_lookup.txt

echo "Done."

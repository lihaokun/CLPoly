#!/bin/bash
# Run all CLPoly automated tests
# Usage: bash test/run_all_tests.sh

set -e

cd "$(dirname "$0")/.."

TESTS=(
    test/test_variable
    test/test_monomial
    test/test_upolynomial
    test/test_polynomial_arith
    test/test_polynomial_div
    test/test_polynomial_gcd
    test/test_polynomial_resultant
    test/test_polynomial_misc
    test/test_realroot
    test/test_charset
    test/test_graph
    test/test_number
    test/test_random
    test/test_parse
    test/test_factorize_zp
    test/test_hensel
    test/test_recombine
    test/test_factorize
    test/test_multivar_helpers
    test/test_wang_lc
    test/test_multivar_hensel
    test/test_factorize_multivar
)

TOTAL_PASS=0
TOTAL_FAIL=0
FAILED_TESTS=()

echo "========================================"
echo "  CLPoly Test Suite"
echo "========================================"

# Build all tests
echo ""
echo "Building tests..."
for t in "${TESTS[@]}"; do
    echo "  Building $t..."
    if ! make "$t" 2>&1 | tail -3; then
        echo "  FAILED to build $t"
        FAILED_TESTS+=("$t (build failed)")
        TOTAL_FAIL=$((TOTAL_FAIL + 1))
        continue
    fi
done

echo ""
echo "Running tests..."
echo "----------------------------------------"

for t in "${TESTS[@]}"; do
    if [ ! -f "$t" ]; then
        echo "SKIP $t (not built)"
        continue
    fi
    echo ""
    echo ">>> $t"
    if ./"$t"; then
        TOTAL_PASS=$((TOTAL_PASS + 1))
    else
        TOTAL_FAIL=$((TOTAL_FAIL + 1))
        FAILED_TESTS+=("$t")
    fi
done

echo ""
echo "========================================"
echo "  Summary: $TOTAL_PASS test files passed, $TOTAL_FAIL failed"
if [ ${#FAILED_TESTS[@]} -gt 0 ]; then
    echo "  Failed:"
    for f in "${FAILED_TESTS[@]}"; do
        echo "    - $f"
    done
fi
echo "========================================"

exit $TOTAL_FAIL

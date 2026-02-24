#!/bin/bash
# Run all CLPoly automated tests
# Usage: bash test/run_all_tests.sh

set -e

cd "$(dirname "$0")/.."

BIN_DIR=_build/debug/bin

TESTS=(
    test_variable
    test_monomial
    test_upolynomial
    test_polynomial_arith
    test_polynomial_div
    test_polynomial_gcd
    test_polynomial_resultant
    test_polynomial_misc
    test_realroot
    test_charset
    test_graph
    test_number
    test_random
    test_parse
    test_factorize_zp
    test_hensel
    test_recombine
    test_factorize
    test_multivar_helpers
    test_wang_lc
    test_multivar_hensel
    test_factorize_multivar
    test_groebner
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
    if ! make "test/$t" 2>&1 | tail -3; then
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
    if [ ! -f "$BIN_DIR/$t" ]; then
        echo "SKIP $t (not built)"
        continue
    fi
    echo ""
    echo ">>> $t"
    if "$BIN_DIR/$t"; then
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

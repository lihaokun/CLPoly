#!/usr/bin/env bash
# 回归测试套件：每 Pass 完成后跑一次，确保无回归。
# 用法: bash tests/run_all_smoke.sh
set -e

cd "$(dirname "$0")/.."

echo "========================================="
echo "  cpp2lean v2 smoke + unit tests"
echo "========================================="

# 单元测试
echo ""
echo "--- Pass 1 单元测试 ---"
python3 tests/test_pass1_parse.py

echo ""
echo "--- Pass 2 单元测试 ---"
python3 tests/test_pass2_ref_elim.py

echo ""
echo "--- Pass 3 单元测试 ---"
python3 tests/test_pass3_lambda_lift.py

echo ""
echo "--- Pass 4 单元测试 ---"
python3 tests/test_pass4_iter.py

echo ""
echo "--- Pass 5 单元测试 ---"
python3 tests/test_pass5_resolve.py

echo ""
echo "--- MIR ir_types 单元测试 ---"
python3 tests/test_mir_types.py

# 烟测（全 65 函数）
echo ""
echo "--- Pass 1 烟测（65 函数）---"
python3 tests/smoke_pass1_full.py 2>&1 | tail -8

echo ""
echo "--- Pass 2 烟测（65 函数）---"
python3 tests/smoke_pass2_full.py 2>&1 | tail -8

echo ""
echo "--- Pass 3 烟测（65 函数, factorize x3）---"
python3 tests/smoke_pass3_full.py 2>&1 | tail -8

echo ""
echo "--- Pass 4 烟测（65 函数, factorize x3）---"
python3 tests/smoke_pass4_full.py 2>&1 | tail -8

echo ""
echo "--- Pass 5 烟测（65 函数, factorize x3）---"
python3 tests/smoke_pass5_full.py 2>&1 | tail -8

# 未来 Pass 会在此添加：
# echo ""
# echo "--- Pass 6 烟测 ---"
# python3 tests/smoke_pass6_full.py 2>&1 | tail -8

echo ""
echo "========================================="
echo "  ✓ 全部通过"
echo "========================================="

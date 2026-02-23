#!/usr/bin/env bash
# run_bench.sh — 构建并运行完整基准套件，输出到 benchmarks/YYYY-MM-DD-HHMMSS.txt
# 包含：crosscheck (FLINT/NTL) + bench_clpoly + bench_comparative
#
# 用法：
#   bash test/run_bench.sh              # 保存到 benchmarks/YYYY-MM-DD-HHMMSS.txt（默认）
#   bash test/run_bench.sh my_label     # 保存到 benchmarks/YYYY-MM-DD-HHMMSS-my_label.txt

set -euo pipefail
cd "$(dirname "$0")/.."

# ---------- 参数 ----------
DATETIME=$(date +%Y-%m-%d-%H%M%S)
LABEL="${1:-}"
if [[ -n "$LABEL" ]]; then
    OUTFILE="benchmarks/${DATETIME}-${LABEL}.txt"
else
    OUTFILE="benchmarks/${DATETIME}.txt"
fi

mkdir -p benchmarks

echo ">>> 构建 debug 库 + crosscheck 二进制..."
make lib/debug/libclpoly.a                    2>&1
make _build/debug/bin/test_crosscheck_flint   2>&1
make _build/debug/bin/test_crosscheck_ntl     2>&1

echo ">>> 构建 release 库 + benchmark 二进制..."
make lib/libclpoly.a                          2>&1
make _build/release/bin/bench_clpoly          2>&1
make _build/release/bin/bench_comparative     2>&1

echo ">>> 开始输出到 $OUTFILE"

# ---------- 文件头 ----------
{
echo "# CLPoly Benchmark — ${DATETIME}"
echo "# Format: fork-based (time + memory per benchmark)"
echo ""
echo ""

# ---------- Crosscheck: FLINT ----------
echo "=== Crosscheck: CLPoly vs FLINT ==="
_build/debug/bin/test_crosscheck_flint 2>&1
echo ""

# ---------- Crosscheck: NTL ----------
echo "=== Crosscheck: CLPoly vs NTL ==="
_build/debug/bin/test_crosscheck_ntl 2>&1
echo ""
echo ""

# ---------- bench_clpoly ----------
echo "========================================"
echo "  bench_clpoly (pure CLPoly)"
echo "========================================"
echo ""
_build/release/bin/bench_clpoly 2>&1
echo ""
echo ""

# ---------- bench_comparative ----------
echo "========================================"
echo "  bench_comparative (CLPoly vs FLINT/NTL)"
echo "========================================"
echo ""
_build/release/bin/bench_comparative 2>&1

} | tee "$OUTFILE"

echo ""
echo ">>> 已保存到 $OUTFILE"

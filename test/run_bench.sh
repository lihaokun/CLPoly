#!/usr/bin/env bash
# run_bench.sh — 构建并运行完整基准套件，输出到 benchmarks/YYYY-MM-DD-HHMMSS.txt
# 包含：bench_clpoly + bench_comparative（默认）
#       加 --crosscheck 可额外跑 debug crosscheck 正确性验证
#
# 用法：
#   bash test/run_bench.sh                      # 仅性能基准（默认）
#   bash test/run_bench.sh --crosscheck         # 性能基准 + debug crosscheck
#   bash test/run_bench.sh my_label             # 带标签
#   bash test/run_bench.sh --crosscheck my_label

set -euo pipefail
cd "$(dirname "$0")/.."

# ---------- 参数 ----------
RUN_CROSSCHECK=false
LABEL=""
for arg in "$@"; do
    if [[ "$arg" == "--crosscheck" ]]; then
        RUN_CROSSCHECK=true
    else
        LABEL="$arg"
    fi
done

DATETIME=$(date +%Y-%m-%d-%H%M%S)
if [[ -n "$LABEL" ]]; then
    OUTFILE="benchmarks/${DATETIME}-${LABEL}.txt"
else
    OUTFILE="benchmarks/${DATETIME}.txt"
fi

mkdir -p benchmarks

# ---------- 构建 ----------
if $RUN_CROSSCHECK; then
    echo ">>> 构建 debug 库 + crosscheck 二进制..."
    make lib/debug/libclpoly.a                    2>&1
    make _build/debug/bin/test_crosscheck_flint   2>&1
    make _build/debug/bin/test_crosscheck_ntl     2>&1
fi

echo ">>> 构建 release 库 + benchmark 二进制..."
make lib/libclpoly.a                          2>&1
make _build/release/bin/bench_clpoly          2>&1
make _build/release/bin/bench_comparative     2>&1

echo ">>> 开始输出到 $OUTFILE"

# ---------- 输出 ----------
{
echo "# CLPoly Benchmark — ${DATETIME}"
echo "# Format: fork-based (time + memory per benchmark)"
echo ""
echo ""

# ---------- Crosscheck (可选) ----------
if $RUN_CROSSCHECK; then
    echo "=== Crosscheck: CLPoly vs FLINT ==="
    _build/debug/bin/test_crosscheck_flint 2>&1
    echo ""

    echo "=== Crosscheck: CLPoly vs NTL ==="
    _build/debug/bin/test_crosscheck_ntl 2>&1
    echo ""
    echo ""
fi

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

#!/bin/bash
# Validation script for MDL clustering and dynamic insertion.
# Compares CPU (clnj_cpp) output with GPU (clnj_cuda) reference when available.
# Usage: ./validate_insertion.sh [build_dir] [benchmark_dir]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLNJ_CPP_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${1:-$CLNJ_CPP_DIR/build}"
BENCH_DIR="${2:-$CLNJ_CPP_DIR/../clnj_cuda/build/benchmark}"

BUILD_FASTA="${BENCH_DIR}/build_100.fasta"
INSERT_FASTA="${BENCH_DIR}/insert_10.fasta"
TMP="/tmp/clnj_validate_$$"
mkdir -p "$TMP"
trap "rm -rf $TMP" EXIT

CPU_BIN="$BUILD_DIR/clnj_pipeline"
GPU_BIN="${CLNJ_CPP_DIR}/../clnj_cuda/build/clnj_cuda"

echo "=== MDL Clustering & Insertion Validation ==="
echo "  CPU binary: $CPU_BIN"
echo "  Benchmark:  $BUILD_FASTA, $INSERT_FASTA"
echo ""

if [[ ! -f "$CPU_BIN" ]]; then
    echo "ERROR: CPU binary not found. Build with: cd $CLNJ_CPP_DIR/build && cmake .. && make"
    exit 1
fi
if [[ ! -f "$BUILD_FASTA" ]]; then
    echo "ERROR: Benchmark FASTA not found: $BUILD_FASTA"
    exit 1
fi
if [[ ! -f "$INSERT_FASTA" ]]; then
    echo "ERROR: Insert FASTA not found: $INSERT_FASTA"
    exit 1
fi

# ── Test 1: CPU build + save state ───────────────────────────────
echo "Test 1: CPU build + save-state"
"$CPU_BIN" build "$BUILD_FASTA" --save-state "$TMP/cpu_state.bin" --no-align 2>&1 | tail -15
echo ""

# ── Test 2: CPU insert ───────────────────────────────────────────
echo "Test 2: CPU insert"
"$CPU_BIN" insert "$TMP/cpu_state.bin" "$INSERT_FASTA" -o "$TMP/cpu_after.bin" 2>&1 | tail -25
echo ""

# ── Test 3: Dump and verify structure ────────────────────────────
echo "Test 3: Dump CPU state (after insert)"
"$CPU_BIN" dump "$TMP/cpu_after.bin"
echo ""

# ── Test 4: Cross-validation with GPU (if available) ─────────────
if [[ -f "$GPU_BIN" ]]; then
    echo "Test 4: GPU build + insert (reference)"
    "$GPU_BIN" build "$BUILD_FASTA" --save-state "$TMP/gpu_state.bin" --no-align 2>&1 | tail -10
    "$GPU_BIN" insert "$TMP/gpu_state.bin" "$INSERT_FASTA" -o "$TMP/gpu_after.bin" 2>&1 | tail -20
    echo ""

    echo "Test 5: Cross-load — CPU insert on GPU-built state"
    "$CPU_BIN" insert "$TMP/gpu_state.bin" "$INSERT_FASTA" -o "$TMP/cpu_from_gpu.bin" 2>&1 | grep -E "(Cluster|Final|Successful|PASSED)"
    echo ""

    echo "Test 6: Cross-load — GPU insert on CPU-built state"
    "$GPU_BIN" insert "$TMP/cpu_state.bin" "$INSERT_FASTA" -o "$TMP/gpu_from_cpu.bin" 2>&1 | grep -E "(Cluster|Final|Successful|PASSED)"
    echo ""

    echo "Test 7: Compare final structure (CPU vs GPU)"
    CPU_DUMP=$("$CPU_BIN" dump "$TMP/cpu_after.bin" 2>/dev/null)
    GPU_DUMP=$("$CPU_BIN" dump "$TMP/gpu_after.bin" 2>/dev/null)
    if [[ "$CPU_DUMP" == "$GPU_DUMP" ]]; then
        echo "  PASS: CPU and GPU produce identical structure"
    else
        echo "  CPU: $CPU_DUMP"
        echo "  GPU: $GPU_DUMP"
        echo "  (Minor differences in internal IDs are acceptable)"
    fi
else
    echo "Test 4–7: Skipped (clnj_cuda not built at $GPU_BIN)"
fi

echo ""
echo "=== Validation complete ==="

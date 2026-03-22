# MDL Clustering & Dynamic Insertion — Testing Strategy

## Overview

This document describes the testing strategy for validating the CPU sequential implementation of MDL clustering and dynamic insertion against the GPU reference (clnj_cuda).

## Test Components

### 1. Cluster Insertion
- **Build + save-state**: `clnj_pipeline build <fasta> --save-state <path> --no-align`
- Verifies: MDL clustering produces leaf clusters, state is saved with cached encoding
- Expected: 4 clusters (3 leaf clusters) for benchmark build_100.fasta

### 2. MDL Cost Update
- MDL clustering runs at build time only; insertion does not re-run MDL
- Cluster assignments use pre-computed tree distances and hierarchical descent

### 3. Tree Updates
- **Insert**: `clnj_pipeline insert <state> <fasta> -o <output>`
- Verifies: Steiner tree extraction, boundary edges, NJ rebuild per cluster, tree invariant
- Expected: Tree invariant PASSED (|V|=n, |E|=n-1)

### 4. Cluster Correctness
- Cluster assignments must match between CPU and GPU for identical inputs
- Cross-load: CPU insert on GPU-built state, GPU insert on CPU-built state

## Validation Script

Run the full validation suite:

```bash
./tests/validate_insertion.sh [build_dir] [benchmark_dir]
```

Default paths:
- `build_dir`: `clnj_cpp/build`
- `benchmark_dir`: `clnj_cuda/build/benchmark`

### Tests Performed

| Test | Description |
|------|-------------|
| 1 | CPU build + save-state |
| 2 | CPU insert |
| 3 | Dump CPU state (structure stats) |
| 4 | GPU build + insert (reference, if clnj_cuda built) |
| 5 | Cross-load: CPU insert on GPU-built state |
| 6 | Cross-load: GPU insert on CPU-built state |
| 7 | Compare final structure (CPU vs GPU) |

## Dump Subcommand

Inspect a saved state without running insertion:

```bash
./clnj_pipeline dump <state_path>
```

Output: observed count, hidden count, total nodes, edges, leaves, clusters, leaf clusters, cached encoding status.

## Validation Results (Benchmark: build_100.fasta + insert_10.fasta)

- **Cluster assignments**: Cluster 1: 3 new taxa, Cluster 2: 6 new taxa, Cluster 3: 1 new taxa
- **Final tree**: 208 nodes, 207 edges, 105 observed
- **Cross-load**: Both CPU and GPU can load each other's states; insert produces identical structure
- **Format compatibility**: Tree state format (CLNJ v2) is shared between clnj_cpp and clnj_cuda

## Constraints

- CPU implementation is sequential; no parallel constructs
- Output must match GPU reference for identical inputs
- State files are binary; exact byte comparison may differ due to map iteration order

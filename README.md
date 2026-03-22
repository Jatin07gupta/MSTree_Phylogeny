# CLNJ Pipeline (C++)

A high-performance C++17 implementation of the **Multi-Model Distance-based Phylogenetic Tree Reconstruction Pipeline** using Latent Tree Graphical Models (LTGM). This pipeline reconstructs phylogenetic trees with both observed (leaf) and hidden (ancestral) nodes from DNA/RNA sequence data.

The implementation is a faithful port of the original Python pipeline, producing numerically identical results on all core algorithms while leveraging C++ for performance.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Pipeline Stages](#pipeline-stages)
3. [Distance Models](#distance-models)
4. [Tree Reconstruction Algorithms](#tree-reconstruction-algorithms)
5. [Key Features](#key-features)
6. [Project Structure](#project-structure)
7. [Requirements](#requirements)
8. [Building](#building)
9. [Usage](#usage)
10. [CLI Reference](#cli-reference)
11. [Examples](#examples)
12. [Diagnostics](#diagnostics)
13. [Pipeline Constants](#pipeline-constants)

---

## Pipeline Overview

The pipeline takes a FASTA file (aligned or unaligned DNA/RNA sequences) as input and produces a phylogenetic tree with observed and hidden (ancestral) nodes. The workflow proceeds through 8 stages:

```
FASTA Input
    |
    v
[Stage 0] Alignment Detection (+ auto-MAFFT if needed)
    |
    v
[Stage 1] Load, Clean, Filter Sequences
    |
    v
[Stage 2] One-Hot Encoding + Alignment Diagnostics
    |                              |
    |                  [Stage 2c] IQ-TREE Model Selection (optional)
    |
    v
[Stage 3] Pairwise Distance Matrix (alignment-based or k-mer)
    |
    v
[Stage 4] Algorithm 2 -- Laminar Family (F_C) and Union Graph (G_U)
    |
    v
[Stage 5] MLVMST -- Initial Tree Topology
    |
    v
[Stage 6] CLNJ -- Cluster-based Neighbor-Joining (local tree refinement)
    |
    v
[Stage 7] Tree Analysis and Statistics
    |
    v
[Stage 8] Pipeline Summary
```

---

## Pipeline Stages

### Stage 0: Alignment Detection

Checks whether the input FASTA file contains already-aligned sequences (uniform length). If not, automatically invokes **MAFFT** with a size-dependent strategy:

| Dataset Size | MAFFT Strategy |
|---|---|
| n < 500 | `--auto` (MAFFT chooses best) |
| 500 <= n < 10,000 | `--retree 2` (FFT-NS-2) |
| n >= 10,000 | `--retree 1` (FFT-NS-1, fastest) |

For RNA data with `--rna-struct`, uses `--kimura 1 --xinsi` for structure-aware alignment.

Skipped entirely when using `--distance-method kmer`.

### Stage 1: Load, Clean, Filter

Parses the FASTA file and applies quality filters:

- **DIVIDER removal**: Skips separator records (names containing "DIVIDER").
- **Short sequence removal**: Sequences with fewer than `min_non_gap` valid (non-gap) nucleotides are removed.
- **Ambiguity filter**: Sequences exceeding `max_ambiguity` fraction of ambiguous bases (N, R, Y, etc.) are removed.
- **Length-percentile cutoff**: Sequences below the `length_percentile_cutoff`-th percentile of valid-site counts are removed.
- **Optional subsampling**: Random subset of N sequences (deterministic with seed=42).
- **U to T normalization**: RNA uracil bases are converted to thymine.

### Stage 2: One-Hot Encoding and Diagnostics

Encodes each nucleotide position into a 4-element vector: A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1]. Gaps and ambiguous characters are encoded as [0,0,0,0].

**Alignment Diagnostics** (Stage 2b) compute:
- Base composition frequencies (gA, gC, gG, gT)
- Transition/transversion ratio (Ts/Tv)
- Mean p-distance across sampled pairs
- Frequency uniformity test

These diagnostics drive the AUTO model selection heuristic.

### Stage 2c: IQ-TREE Model Selection (optional)

When `--model IQTREE` is specified, invokes `iqtree2` externally for likelihood-based model selection using BIC. The selected IQ-TREE model is mapped to the pipeline's supported models:

| IQ-TREE Model | Pipeline Model |
|---|---|
| JC, JC+... | JC69 |
| K2P, K80, K2P+... | K2P |
| HKY, TN, TN93, TNe, TIM, ... | TN93 |
| GTR, SYM, TPM, TVM, ... | TN93 |

Gamma rate heterogeneity parameters are extracted when present.

### Stage 3: Distance Matrix

Computes an n x n pairwise distance matrix using the selected model. Two distance computation paths:

**Alignment-based** (default): Uses one-hot encoded data with pairwise deletion for gap handling. A `min_shared_sites` threshold ensures reliable estimates.

**K-mer based** (`--distance-method kmer`): Alignment-free Mash distance using MinHash sketching. Does not require aligned sequences.

### Stage 4: Algorithm 2 (F_C and G_U)

Constructs the laminar family F_C and union graph G_U from the distance matrix using a Kruskal-like sweep with a Disjoint Set Union (DSU) data structure. Edges are sorted by distance and processed to build hierarchical cluster structure.

### Stage 5: MLVMST (Algorithm 3)

Builds the Minimum Leaves Vertex-Order Minimum Spanning Tree. Computes `delta_max` scores for each taxon and constructs the initial tree topology using a vertex-ordered Kruskal MST algorithm.

### Stage 6: CLNJ (Cluster-based Neighbor-Joining)

The core reconstruction step. Iterates over internal nodes of the MLVMST (sorted by degree, descending) and performs local tree refinement:

1. For each internal node, extract its neighborhood (center + adjacent nodes).
2. Run a local NJ variant (NJ, MFNJ, or BioNJ) on the neighborhood.
3. Splice the locally reconstructed subtree back into the global tree.
4. Validate tree invariant (|E| = |V| - 1) after each splice.

A **context-aware distance oracle** provides distances: observed-observed pairs use the original distance matrix; pairs involving hidden nodes use BFS tree-path distances for additivity guarantees.

Branch-length redistribution enforces non-negative edge weights: when delta_i < floor, set delta_i = floor and delta_j = d_ij - floor.

### Stage 7: Tree Analysis

Computes comprehensive statistics on the final tree:
- Node counts (observed, hidden, total)
- Edge weight statistics (min, max, mean, median, std)
- Edge type breakdown (obs-obs, obs-hidden, hidden-hidden)
- Degree distribution (leaves, max degree)
- Tree validity check

Runs on both the floored (min_branch_length = 1e-6) and raw (min_branch_length = 0) trees for comparison.

### Stage 8: Pipeline Summary

Prints a complete summary including dataset info, timings, distance matrix statistics, and tree characteristics.

---

## Distance Models

| Model | Description | When to Use |
|---|---|---|
| **JC69** | Jukes-Cantor (1969). Equal base frequencies, equal substitution rates. Simplest model. | Low divergence, uniform composition |
| **K2P** | Kimura 2-parameter (1980). Distinguishes transitions from transversions. | Moderate divergence, Ts/Tv bias |
| **TN93** | Tamura-Nei (1993). Unequal base frequencies, two transition rates, one transversion rate. | Unequal composition, high Ts/Tv ratio |
| **LOGDET** | Log-determinant (paralinear). Model-free, handles composition heterogeneity across lineages. | Compositional bias between taxa |
| **AUTO** | Heuristic selection using 5 alignment diagnostics: invariant sites, site-rate CV, gamma alpha, per-taxon composition, Ts/Tv ratio. | General purpose |
| **IQTREE** | External likelihood-based model selection via IQ-TREE2 (BIC criterion). | Best-fit model when IQ-TREE is available |

All alignment-based models support **gamma rate heterogeneity correction** via the `--gamma-alpha` parameter. The gamma transform is: `d = alpha * (x^(-1/alpha) - 1)`.

---

## Tree Reconstruction Algorithms

| Algorithm | Flag | Description |
|---|---|---|
| **Standard NJ** | `nj` | Classic Neighbor-Joining (Saitou & Nei 1987). Binary tree, Q-matrix based pair selection. |
| **MFNJ** | `mfnj` | MultiFurcating Neighbor-Joining (Fernandez et al. 2023). Detects Q-matrix ties and creates polytomies instead of arbitrary zero-length binary splits. Produces a unique tree independent of input order. **Default algorithm.** |
| **BioNJ** | `bionj` | Variance-minimizing NJ (Gascuel 1997). Uses optimal weighting parameter lambda* to minimize variance of distance estimates during agglomeration. |

---

## Key Features

- **Modular architecture**: Each pipeline stage is a separate, self-contained module.
- **Multi-model distance support**: 5 alignment-based models + alignment-free k-mer distance.
- **Automatic alignment**: MAFFT integration with size-dependent strategies.
- **IQ-TREE integration**: External likelihood-based model selection.
- **Three NJ variants**: Standard NJ, MFNJ (polytomy-aware), BioNJ (variance-minimizing).
- **Context-aware distance oracle**: Guarantees additive distances for hidden-node neighborhoods.
- **Branch-length redistribution**: Non-negative edge weights via redistribution (not clamping).
- **Diagnostic counters**: Track negative values and zero-length edges across all NJ variants.
- **Gamma rate heterogeneity**: Optional correction for among-site rate variation.
- **Pairwise deletion**: Gap handling with configurable minimum shared sites.
- **Deterministic subsampling**: Reproducible results with fixed seed (42).
- **Comprehensive tree analysis**: Edge statistics, degree distribution, validity checks.
- **Dual-run comparison**: Automatically compares floored vs. raw branch lengths.
- **Zero external dependencies at build time**: Eigen and CLI11 are fetched automatically by CMake.

---

## Project Structure

```
clnj_cpp/
├── CMakeLists.txt                  # Build configuration (CMake 3.14+)
├── README.md                       # This file
│
├── include/                        # Header files
│   ├── types.h                     # Constants, type aliases, structs, enums
│   ├── fasta_parser.h              # FASTA parsing and sequence cleaning
│   ├── one_hot.h                   # One-hot encoding of nucleotide sequences
│   ├── distance.h                  # Distance model computation and AUTO selection
│   ├── kmer_distance.h             # Alignment-free k-mer (Mash) distance
│   ├── alignment_utils.h           # Alignment detection and MAFFT wrapper
│   ├── iqtree_interface.h          # IQ-TREE2 model selection interface
│   ├── algorithm2.h                # Algorithm 2: F_C and G_U construction
│   ├── mlvmst.h                    # MLVMST (Algorithm 3): initial tree topology
│   ├── distance_oracle.h           # Context-aware distance queries (D_obs + BFS)
│   ├── nj.h                        # Standard Neighbor-Joining
│   ├── mfnj.h                      # MultiFurcating Neighbor-Joining
│   ├── bionj.h                     # BioNJ (variance-minimizing NJ)
│   ├── clnj.h                      # CLNJ dispatcher (iterative local refinement)
│   ├── tree_analysis.h             # Final tree statistics and reporting
│   └── pipeline.h                  # Pipeline orchestration and PipelineArgs
│
├── src/                            # Source files
│   ├── main.cpp                    # CLI entry point (CLI11 argument parsing)
│   ├── pipeline.cpp                # End-to-end pipeline orchestration (Stages 0-8)
│   ├── fasta_parser.cpp            # FASTA I/O, filtering, cleaning
│   ├── one_hot.cpp                 # Nucleotide -> one-hot lookup table encoding
│   ├── distance.cpp                # JC69, K2P, TN93, LogDet, AUTO, gamma correction
│   ├── kmer_distance.cpp           # Canonical k-mers, FNV-1a hash, MinHash, Mash
│   ├── alignment_utils.cpp         # is_aligned(), align_with_mafft()
│   ├── iqtree_interface.cpp        # iqtree2 invocation and .iqtree output parsing
│   ├── algorithm2.cpp              # DSU-based F_C/G_U construction
│   ├── mlvmst.cpp                  # delta_max computation, vertex-ordered MST
│   ├── distance_oracle.cpp         # BFS tree-path distance, DistanceOracle class
│   ├── nj.cpp                      # Standard NJ with diagnostics
│   ├── mfnj.cpp                    # MFNJ with tie detection and polytomies
│   ├── bionj.cpp                   # BioNJ with variance matrix and lambda*
│   ├── clnj.cpp                    # CLNJ iteration, splicing, tree invariant checks
│   └── tree_analysis.cpp           # Edge/node statistics, degree analysis
│
└── build/                          # Build output directory (created by CMake)
```

---

## Requirements

### Build Requirements

| Requirement | Version | Notes |
|---|---|---|
| **C++ compiler** | C++17 support | GCC 7+, Clang 5+, or MSVC 2017+ |
| **CMake** | 3.14+ | Build system |
| **Internet access** | (first build only) | To download Eigen and CLI11 via FetchContent |

The following libraries are **automatically downloaded** during the first CMake build:

| Library | Version | Purpose |
|---|---|---|
| **Eigen** | 3.4.0 | Dense matrix operations (distance matrices, one-hot encoding) |
| **CLI11** | 2.4.1 | Command-line argument parsing |

### Runtime Requirements (Optional External Tools)

These are only needed if you use the corresponding features:

| Tool | Required For | Installation |
|---|---|---|
| **MAFFT** | Automatic alignment of unaligned FASTA input | `sudo apt install mafft` or from [mafft.cbrc.jp](https://mafft.cbrc.jp/alignment/software/) |
| **IQ-TREE2** | `--model IQTREE` model selection | Download from [iqtree.org](http://www.iqtree.org/) or GitHub releases |

Both tools must be accessible on `PATH` (or symlinked into your environment's `bin/` directory).

### System Requirements

- **OS**: Linux (tested), macOS (should work), Windows (untested)
- **RAM**: Depends on dataset size. The n x n distance matrix requires O(n^2) memory.
- **Disk**: Minimal. Temporary files are created for MAFFT/IQ-TREE and cleaned up automatically.

---

## Building

```bash
# 1. Navigate to the project directory
cd clnj_cpp

# 2. Create and enter the build directory
mkdir -p build && cd build

# 3. Configure with CMake (downloads Eigen and CLI11 on first run)
cmake ..

# 4. Build (use all available cores)
cmake --build . -j$(nproc)
```

The executable `clnj_pipeline` will be created in the `build/` directory.

To rebuild after code changes:

```bash
cd build && cmake --build . -j$(nproc)
```

---

## Usage

### Basic Syntax

```bash
./clnj_pipeline <fasta_file> [min_non_gap] [subsample_n] [tree_algo] [options]
```

### Minimal Run

```bash
./clnj_pipeline /path/to/sequences.fasta
```

This uses all defaults: JC69 model, MFNJ algorithm, no subsampling, alignment-based distance.

### Typical Run with Options

```bash
./clnj_pipeline /path/to/sequences.fasta 100 0 mfnj \
    --model TN93 \
    --gamma-alpha 0.5 \
    --report-negatives \
    --trace-zero-edges
```

---

## CLI Reference

### Positional Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `fasta` | Yes | — | Path to FASTA file (aligned or unaligned) |
| 2 | `min_non_gap` | No | 100 | Minimum valid (non-gap) nucleotides per sequence |
| 3 | `subsample_n` | No | 0 | Subsample to N sequences; 0 = use all |
| 4 | `tree_algo` | No | `mfnj` | Tree algorithm: `nj`, `mfnj`, or `bionj` |

### Named Options

| Flag | Short | Default | Description |
|---|---|---|---|
| `--model` | `-m` | `JC69` | Distance model: `JC69`, `K2P`, `TN93`, `LOGDET`, `AUTO`, `IQTREE` |
| `--gamma-alpha` | `-g` | (none) | Gamma shape parameter for rate heterogeneity |
| `--no-align` | | false | Skip MAFFT alignment even if input is unaligned |
| `--rna-struct` | | false | Use MAFFT X-INS-i for RNA structure-aware alignment |
| `--min-shared-sites` | | 100 | Minimum shared valid sites for a valid pairwise distance |
| `--max-ambiguity` | | 0.5 | Maximum fraction of ambiguity codes per sequence |
| `--length-percentile-cutoff` | | 5.0 | Remove sequences below this percentile of valid-site count |
| `--distance-method` | | `alignment` | Distance method: `alignment` or `kmer` |
| `--kmer-size` | | 16 | K-mer size for alignment-free distance |
| `--sketch-size` | | 1000 | MinHash sketch size for k-mer distance |
| `--report-negatives` | | false | Track and report negative-value diagnostic counters |
| `--trace-zero-edges` | | false | Track and report zero-length edges during CLNJ |

---

## Examples

### 1. Basic run with default settings (JC69 + MFNJ)

```bash
./clnj_pipeline data/carnivora_cytb_aligned.fasta
```

### 2. Using K2P model with BioNJ

```bash
./clnj_pipeline data/sequences.fasta 100 0 bionj --model K2P
```

### 3. AUTO model selection with gamma correction

```bash
./clnj_pipeline data/sequences.fasta 100 0 mfnj --model AUTO --gamma-alpha 0.5
```

### 4. IQ-TREE model selection (requires iqtree2 on PATH)

```bash
./clnj_pipeline data/sequences.fasta 100 0 mfnj --model IQTREE
```

### 5. Alignment-free k-mer distance

```bash
./clnj_pipeline data/unaligned_sequences.fasta 50 0 mfnj \
    --distance-method kmer --kmer-size 16 --sketch-size 1000
```

### 6. Subsampling to 50 sequences with diagnostics

```bash
./clnj_pipeline data/large_dataset.fasta 100 50 mfnj \
    --model TN93 --report-negatives --trace-zero-edges
```

### 7. Unaligned input with automatic MAFFT alignment

```bash
./clnj_pipeline data/unaligned.fasta 100 0 mfnj --model JC69
```

### 8. RNA sequences with structure-aware alignment

```bash
./clnj_pipeline data/16s_rRNA.fasta 100 0 mfnj --rna-struct --model AUTO
```

---

## Diagnostics

Two optional diagnostic modes provide detailed insight into the numerical behavior of the tree reconstruction:

### Negative-Value Counters (`--report-negatives`)

Tracks how many times raw distance or branch-length values fall below -1e-12 during local NJ reconstruction. This is normal for non-additive distance matrices but useful for assessing data quality.

**Tracked categories:**

| Category | Algorithm | Meaning |
|---|---|---|
| `initial_pair` | All | Oracle distance for a 2-node neighborhood was negative |
| `oracle_dist` | All | Pairwise oracle distance in local D matrix was negative |
| `delta_i` | NJ, BioNJ | Raw branch length delta_i was negative before redistribution |
| `delta_j` | NJ, BioNJ | Raw branch length delta_j was negative before redistribution |
| `NJ_update_d_uk` | NJ | Distance update for new hidden node produced negative value |
| `BioNJ_update_d_uk` | BioNJ | Weighted distance update produced negative value |
| `final_edge` | All | Final edge between last 2 active nodes was negative |
| `mfnj_branch` | MFNJ | Branch length in MFNJ group (with complement) was negative |
| `mfnj_branch_final` | MFNJ | Branch length in final MFNJ group (no complement) was negative |
| `mfnj_update_d_uk` | MFNJ | Distance update for new MFNJ hidden node was negative |

### Zero-Edge Log (`--trace-zero-edges`)

Records every edge created with weight < 1e-9, providing a detailed log of near-zero branch lengths.

**Case types:**

| Case | Meaning |
|---|---|
| `initial_pair` | A 2-node neighborhood had near-zero oracle distance |
| `create_hidden` | A delta_i or delta_j was near-zero after redistribution (NJ/BioNJ) |
| `final_edge` | The final edge between the last 2 active nodes was near-zero |
| `final_star` | A branch in the 3-node star topology was near-zero (MFNJ) |
| `mfnj_polytomy` | A branch in an MFNJ polytomy group was near-zero |

Each entry records: the edge endpoints (a, b), the weight (w), raw values before clamping, the reason string, and the full node list of the local neighborhood.

---

## Pipeline Constants

All numerical constants match the Python implementation exactly:

| Constant | Value | Purpose |
|---|---|---|
| `MAX_DIST` | 10.0 | Maximum allowed pairwise distance (saturation cap) |
| `MIN_BRANCH` | 1e-6 | Minimum branch length floor (default for floored run) |
| `EPS` | 1e-12 | Floating-point comparison tolerance |
| `TIE_TOL` | 1e-9 | Q-matrix tie tolerance for MFNJ |
| `SEED` | 42 | Random seed for reproducible subsampling |
| `LOG_FLOOR` | 1e-300 | Floor to prevent log(0) in distance corrections |
| `FREQ_FLOOR` | 1e-10 | Floor for base frequency estimates |

---

## GitHub repository (`MSTree_Phylogeny`)

This project is intended to live at **`https://github.com/Jatin07gupta/MSTree_Phylogeny`**.

**First-time publish**

1. On GitHub: **New repository** → name **`MSTree_Phylogeny`** → leave it **empty** (no README/license).
2. From this directory:

   ```bash
   git init
   git add .
   git status   # confirm build/ is not listed
   git commit -m "Initial commit: MSTree phylogeny C++ pipeline"
   git branch -M main
   git remote add origin https://github.com/Jatin07gupta/MSTree_Phylogeny.git
   git push -u origin main
   ```

Use SSH instead of HTTPS if you prefer: `git@github.com:Jatin07gupta/MSTree_Phylogeny.git`.

---

## References

- Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406-425.
- Fernandez, M., Segura-Alabart, N., & Serratosa, F. (2023). The MultiFurcating Neighbor-Joining Algorithm for Reconstructing Polytomic Phylogenetic Trees. *Journal of Molecular Evolution*.
- Gascuel, O. (1997). BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Molecular Biology and Evolution*, 14(7), 685-695.
- Jukes, T. H., & Cantor, C. R. (1969). Evolution of Protein Molecules. *Mammalian Protein Metabolism*, 3, 21-132.
- Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions. *Journal of Molecular Evolution*, 16(2), 111-120.
- Tamura, K., & Nei, M. (1993). Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA. *Molecular Biology and Evolution*, 10(3), 512-526.
- Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 30(4), 772-780.
- Nguyen, L. T., et al. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Molecular Biology and Evolution*, 32(1), 268-274.

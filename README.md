# CLNJ / MSTree phylogeny pipeline (C++)

This repository implements a **complete distance-based phylogenetic reconstruction system** in **C++17**: from **FASTA input** through **multi-model pairwise distances**, **Algorithm 2** (laminar family + union graph), **MLVMST**, and **CLNJ** (local NJ / MFNJ / BioNJ with hidden nodes), to **persistent tree state** and **dynamic (online) insertion** of new taxa without rebuilding the full pipeline from scratch.

The design is aligned with the **`clnj_fresh`** Python stack (`run_snj_pipeline.py`, `distance.py`, `algorithm2_FC_GU_fixed.py`, `mlvmst.py`, `clnj_clean_implementation_no_negative.py`, `distance_oracle.py`): same **stage ordering**, **distance models**, and **numerical constants** (`types.h`), with C++ providing a **single binary** and better scalability.

---

## Table of contents

1. [What “the pipeline” includes](#what-the-pipeline-includes)
2. [End-to-end lifecycle](#end-to-end-lifecycle)
3. [Phase 1 — Batch reconstruction (Stages 0–8)](#phase-1--batch-reconstruction-stages-08)
4. [Phase 1 — Stateful preparation (Stage 9)](#phase-1--stateful-preparation-stage-9)
5. [Phase 2 — Dynamic (online) insertion](#phase-2--dynamic-online-insertion)
6. [CLI: `build`, `insert`, and `dump`](#cli-build-insert-and-dump)
7. [Persistent state: `TreeState` binary format](#persistent-state-treestate-binary-format)
8. [Module catalog (every source file)](#module-catalog-every-source-file)
9. [Implementation note: sequences, one-hot, and planned direct encoding](#implementation-note-sequences-one-hot-and-planned-direct-encoding)
10. [Experimental: alignment-free k-mer distance](#experimental-alignment-free-k-mer-distance)
11. [Distance models](#distance-models)
12. [Tree reconstruction algorithms](#tree-reconstruction-algorithms)
13. [Key features](#key-features)
14. [Requirements and building](#requirements-and-building)
15. [Examples](#examples)
16. [Diagnostics](#diagnostics)
17. [Pipeline constants](#pipeline-constants)
18. [GitHub (`MSTree_Phylogeny`)](#github-mstree_phylogeny)
19. [References](#references)

---

## What “the pipeline” includes

| Component | Role |
|-----------|------|
| **Batch reconstruction** | Stages **0–8**: alignment → clean FASTA → diagnostics & model → **D** → Algorithm 2 → MLVMST → CLNJ → analysis → summary. |
| **Stateful preparation (Stage 9)** | After a satisfactory tree, **MDL clustering** builds a **hierarchy of observed-taxa clusters** on the tree; **`TreeState`** serializes sequences, **D**, topology, model parameters, cluster tree, and **cached encodings** for fast updates. |
| **Dynamic insertion** | **`insert`** loads `TreeState`, assigns each new sequence to a cluster, **extends** **D**, **locally** re-runs NJ/MFNJ/BioNJ on the affected subgraph, **splices** back, updates cluster membership, and **saves** state again. |

Together, these define the **full pipeline** as a **living** system: initial tree + **incremental evolution**.  

If you only need a **one-off** tree and will **never** add sequences, you may run **`build` without `--save-state`**; you then use **only** Phase 1 Stages 0–8. That is a **subset** of the system, not a separate product.

---

## End-to-end lifecycle

```
  FASTA (DNA/RNA)
        |
        v
  +-------------------+
  | Stage 0: Align?   |  MAFFT if unaligned (skipped in k-mer mode only)
  +-------------------+
        |
        v
  +-------------------+
  | Stage 1: Clean    |  Filters, subsampling, U->T
  +-------------------+
        |
        v
  +-------------------+
  | Stage 2: Diag.    |  Base freqs, Ts/Tv, p-dist; optional IQ-TREE
  +-------------------+
        |
        v
  +-------------------+
  | Stage 3: D (nxn)  |  Substitution models + pairwise deletion
  +-------------------+
        |
        v
  +-------------------+
  | Stage 4: Alg. 2   |  F_C laminar family, G_U union graph
  +-------------------+
        |
        v
  +-------------------+
  | Stage 5: MLVMST   |  Initial scaffold (Kruskal/Boruvka backend)
  +-------------------+
        |
        v
  +-------------------+
  | Stage 6: CLNJ     |  Local NJ/MFNJ/BioNJ + distance oracle
  +-------------------+
        |
        v
  +-------------------+
  | Stage 7: Analysis |  Stats; floored vs raw branch comparison
  +-------------------+
        |
        v
  +-------------------+
  | Stage 8: Summary  |
  +-------------------+
        |
        |    build --save-state <path>
        v
  +-------------------+
  | Stage 9: MDL      |  Cluster tree on observed taxa + save TreeState
  +-------------------+
        |
        v
   .state file(s)
        |
        |    insert <state> <new.fasta>  (repeat as needed)
        v
  +-------------------+
  | Phase 2: Insert   |  Assign cluster -> extend D -> local rebuild -> save
  +-------------------+
```

---

## Phase 1 — Batch reconstruction (Stages 0–8)

### Stage 0 — Alignment detection and MAFFT

**Files:** `alignment_utils.{h,cpp}`, `fasta_parser.cpp` (for counting sequences before MAFFT).

- Tests whether all sequences in the FASTA have the **same length** (`is_aligned`).
- If **unaligned** and **`--no-align` is not set**, runs **MAFFT** with **size-dependent** options (small \(n\): `--auto`; medium: `--retree 2`; large: `--retree 1`).
- **`--rna-struct`**: structure-aware MAFFT (X-INS-i style).
- In **`--distance-method kmer`** only, alignment is **skipped** (experimental path; see [Experimental k-mer](#experimental-alignment-free-k-mer-distance)).

### Stage 1 — Load, clean, filter

**Files:** `fasta_parser.{h,cpp}`.

- **`parse_fasta`**: reads `FastaRecord` {name, sequence}.
- **`load_clean_fasta`**: removes DIVIDER-like records, enforces **`min_non_gap`**, **`max_ambiguity`**, **`length_percentile_cutoff`**; optional uniform length check (relaxed for k-mer).
- **Subsampling** (if `subsample_n > 0`): deterministic shuffle with **`SEED = 42`** in `pipeline.cpp`.
- **U → T** normalization for RNA.

**Output:** `CleanResult`: parallel vectors **`names`**, **`sequences`**, alignment length **`m`** (for alignment mode).

### Stage 2 — Alignment diagnostics and model resolution

**Files:** `one_hot.{h,cpp}`, `distance.{h,cpp}`, `iqtree_interface.{h,cpp}`.

- **One-hot encoding** of cleaned sequences (`one_hot_encode`) → tensor **(n, m, 4)** for fast column statistics.
- **`compute_pair_counts`**: empirical **gA, gC, gG, gT**, **Ts/Tv**, mean **p-distance**, composition flags — used by **AUTO** (`select_model`).
- **`--model IQTREE`**: runs **`run_iqtree_model_selection`** (subprocess to **IQ-TREE2**), parses `.iqtree` output, maps to **JC69 / K2P / TN93**, may set **gamma** from IQ-TREE.

### Stage 3 — Pairwise distance matrix **D**

**Files:** `distance.{h,cpp}`, `one_hot.h` (input type for alignment path).

- **`compute_distance_matrix`**: full **n × n** matrix under **JC69, K2P, TN93, LogDet, AUTO** with **pairwise deletion** and **`min_shared_sites`**.
- Optional **gamma** correction via **`jc69_correct` / `k2p_correct` / `tn93_correct`** and **`--gamma-alpha`**.
- Distances capped by **`MAX_DIST`**.

**Conceptual input:** aligned **character** sequences; **current** implementation uses **OneHotData** as the concrete representation (see [Implementation note](#implementation-note-sequences-one-hot-and-planned-direct-encoding)).

### Stage 4 — Algorithm 2 (**F_C**, **G_U**)

**Files:** `algorithm2.{h,cpp}`.

- **`construct_FC_and_GU(D)`**: Kruskal-like sweep over sorted pairwise edges with **DSU** to build:
  - **F_C**: laminar family of clusters.
  - **G_U**: edge set of the union graph.

### Stage 5 — MLVMST (Algorithm 3)

**Files:** `mlvmst.{h,cpp}`.

- **`build_mlvmst`**: computes **delta_max**, builds **vertex-ordered MST** adjacency on the Algorithm 2 structure.
- **`--mst-algorithm`**: `kruskal` or `boruvka`.

### Stage 6 — CLNJ

**Files:** `clnj.{h,cpp}`, `distance_oracle.{h,cpp}`, `nj.{h,cpp}`, `mfnj.{h,cpp}`, `bionj.{h,cpp}`.

- **`clnj_clean`**: iterates over internal nodes of the MLVMST (by degree), extracts neighborhoods, builds **`DistanceOracle`** (matrix distances for leaves; **BFS path sums** when hidden nodes involved), runs **`nj_local` / `mfnj_local` / `bionj_local`**, splices subtrees, **redistributes** negative branch proposals to keep **non-negative** edges.
- Flags **`--report-negatives`**, **`--trace-zero-edges`** fill **`ClnjResult`** diagnostics.

### Stage 7 — Tree analysis

**Files:** `tree_analysis.{h,cpp}`.

- **`analyze_tree`**: node/edge counts, weight summaries, obs/hidden/hidden-hidden edge classes, degree distribution, validity (**|E| = |V| − 1** for tree).
- Invoked for **floored** (`MIN_BRANCH`) and **raw** (`0`) CLNJ outputs.

### Stage 8 — Pipeline summary

**File:** `pipeline.cpp`.

- Prints timings, **D** statistics, tree summaries (console log suitable for experiments / supplementary logs).

---

## Phase 1 — Stateful preparation (Stage 9)

**Files:** `mdl_clustering.{h,cpp}`, `tree_state.{h,cpp}`, `pipeline.cpp` (orchestration).

When **`build` is run with `--save-state <path>`**, after Stages 6–7 the pipeline runs **Stage 9**:

1. **`precompute_tree_distances`**  
   For each **observed** leaf **i**, BFS over the **current CLNJ tree** with **edge_weights** to fill an **n_obs × n_obs** cache of **patristic** distances between observed taxa. This is the metric geometry used for MDL on the **observed** set.

2. **`mdl_cluster_tree(...)`**  
   Recursively partitions the tree into **clusters** of observed taxa using **Minimum Description Length (MDL)**-style scores (`description_length`, parameter **`mj`** from **`--mdl-mj`**, merge control **`--mdl-merge-threshold`**). Each **`Cluster`** stores:
   - **`observed_members`**, **`center_node`** (tree node acting as representative),
   - **`parent_cluster_id` / `child_cluster_ids`** for a **hierarchy**,
   - **`node_to_cluster`** maps each observed index → cluster id.

3. **`TreeState` assembly** (`pipeline.cpp`)  
   Copies **names**, **sequences**, **aln_len**, full **D**, **adjacency**, **edge_weights**, **hidden_info**, resolved **`DistModel`**, **gamma**, **base frequencies**, **`min_shared_sites`**, **`tree_algo`**, **`distance_method`**, k-mer metadata (for consistency), **`cluster_tree`**, and for alignment mode builds **`cached_nuc_idx`** / **`cached_valid`** (per-site A/C/G/T index + validity) so insertion avoids rebuilding full one-hot tensors.

4. **`save_tree_state`** (`tree_state.cpp`)  
   Writes a **versioned binary** (magic **`CLNJ`**, **`VERSION = 2`**) containing all of the above for **`load_tree_state`** on **`insert`**.

**Stage 9 is the bridge** between “tree as a printed result” and **Phase 2**. Without **`--save-state`**, **Phase 2 cannot run** because there is no cluster hierarchy or serialized **D**/topology.

---

## Phase 2 — Dynamic (online) insertion

**Files:** `online_insertion.{h,cpp}`, `distance.{h,cpp}`, `tree_state.{h,cpp}`, `pipeline.cpp` (`run_insertion_pipeline`).

**Entry point:** `clnj_pipeline insert <state> <new.fasta> [-o updated.state]`

### Purpose

Add **one or many** new sequences (same **alignment length** **`aln_len`** as stored; typically **pre-aligned** to the same reference) into the **existing** tree while:

- Preserving the **substitution model** and **gamma** used at build time.
- **Extending** **D** only where needed (to cluster centers, cluster members, and between new taxa in the same cluster batch).
- **Locally** re-optimizing topology with the **same** NJ/MFNJ/BioNJ **local** functions as CLNJ.

### Algorithm (high level)

Implemented in **`insert_batch`** (`online_insertion.cpp`):

1. **Load state**  
   **`load_tree_state`** → **`TreeState`** must contain a **non-empty** **`cluster_tree`**.

2. **Encode new sequences**  
   **`encode_sequence_flat`**: per-site **nucleotide index** (−1 if gap/ambiguous) and **valid** flag — same information as one-hot rows, compact.

3. **Phase A — Distance to cluster centers**  
   For each new sequence, compute distances to every **unique** **`center_node`** among clusters using either:
   - **`compute_distance_to_subset_cached`** (if **`cached_nuc_idx`** present — **fast path**), or  
   - **`compute_distance_to_subset`** after **`one_hot_encode`** of existing sequences (**cache miss**).

4. **Hierarchical cluster assignment**  
   **`find_cluster_hierarchical`**: walk the **cluster tree** from the root, comparing the new taxon’s distances to **child cluster centers** (with fallback leaf scan). Produces **`target_cluster_id`** per new taxon.

5. **Remap hidden node ids**  
   **`remap_hidden_nodes`**: shift hidden-node indices by **`batch_size`** so new observed ids **`old_n .. old_n+batch−1`** stay contiguous **0 .. new_n−1** for observed ordering conventions.

6. **Extend **D****  
   - Copy old **D** into **`new_D`** (top-left block).
   - For each new taxon vs **old members** of its assigned cluster: pairwise model distances (cached or one-hot path).
   - For **multiple** new taxa in the **same** cluster: fill **new–new** blocks using the same **JC69/K2P/TN93** formulas on **column-wise** counts (inline in `insert_batch`).

7. **Augment graph**  
   Append **isolated** observed nodes for each new taxon to **adjacency**; set **`n_observed = new_n`**; extend **cache** arrays if present.

8. **Per-cluster local rebuild**  
   For each cluster that received new taxa:
   - **`old_members`** = previous **observed_members**; **`all_members`** = union with new ids.
   - **`find_steiner_tree`**: prune the current tree to the minimal subtree spanning **old_members** (Steiner nodes may include **hidden** nodes).
   - **`find_boundary_edges`**: edges from Steiner nodes to the **rest** of the tree (attachment points).
   - **`remove_steiner_tree`**: cut that subgraph while keeping boundary metadata.
   - Build **`DistanceOracle`** on **updated** **D** and working **adjacency** / **edge_weights**.
   - **`local_fn`** (`nj_local` / `mfnj_local` / `bionj_local`) on **`member_list`** — same local reconstruction as CLNJ.
   - **`install_nj_and_reconnect`** (or direct edge install for degenerate cases): splice local tree back and reconnect **boundary** edges.
   - Update **`cluster.observed_members`** and **`node_to_cluster`** for new ids.

9. **Results**  
   Each taxon gets **`InsertionResult`**: **`new_obs_id`**, **`target_cluster_id`**, **`success`**, **`message`**.

10. **Save**  
    **`save_tree_state`** to **`--output-state`** or overwrite input path.

### Invariants and assumptions

- **Alignment length** of new sequences must match **`state.aln_len`** (padding/truncation behavior follows `encode_sequence_flat` and FASTA content).
- **Cluster tree** must have been built from the **same** tree that produced **D** and **adjacency** (consistent **Stage 9**).
- Insertion **does not** re-run Algorithm 2 / MLVMST globally; it is **local** to assigned clusters — appropriate for **incremental** updates when the MDL partition is meaningful.

---

## CLI: `build`, `insert`, and `dump`

**File:** `main.cpp` (CLI11).

| Command | Purpose |
|---------|---------|
| **`build`** | Full Stages 0–8; **Stage 9** if **`--save-state`** is set. |
| **`insert`** | Phase 2: load state, insert FASTA records, save updated state. |
| **`dump`** | Load state and print **summary stats** (validation / debugging). |

**`build` positional arguments:**  
`fasta` (required), `min_non_gap` (default 100), `subsample_n` (default 0), `tree_algo` (`nj` \| `mfnj` \| `bionj`, default `mfnj`).

**Important flags for the full pipeline:**

- **`--save-state <path>`** — required to enable **Phase 2** later.
- **`--mdl-mj`**, **`--mdl-merge-threshold`** — tune **MDL** cluster tree (Stage 9).

Run **`clnj_pipeline build --help`** and **`clnj_pipeline insert --help`** for the complete flag list.

---

## Persistent state: `TreeState` binary format

**File:** `tree_state.cpp`.

| Field | Meaning |
|-------|---------|
| Magic / **VERSION** | **`0x434C4E4A` ("CLNJ")**, version **2** — readers must match. |
| **n_observed**, **next_hidden_id**, **aln_len** | Dimensions and hidden-id cursor. |
| **model**, **gamma_alpha**, **gA..gT**, **min_shared_sites** | Distance model snapshot. |
| **tree_algo**, **distance_method**, **kmer_size**, **sketch_size** | Reproducibility / consistency. |
| **names**, **sequences** | Full alignment strings. |
| **D** | Dense **n × n** `double` matrix (row-major Eigen storage). |
| **adjacency**, **edge_weights**, **hidden_info** | CLNJ tree topology and branch lengths. |
| **cluster_tree** | **MDL** hierarchy: clusters, members, centers, parent/child links, **node_to_cluster**. |
| **cached_nuc_idx**, **cached_valid** | Flat **n × aln_len** arrays for fast insertion (alignment mode). |

**Security / portability:** binary layout is **native-endian**; treat state files as **opaque** blobs produced and consumed only by this tool (or compatible versions).

---

## Module catalog (every source file)

| Path | Responsibility |
|------|----------------|
| **`main.cpp`** | CLI11 app: subcommands **`build`**, **`insert`**, **`dump`**; parses **`PipelineArgs`** / **`InsertionArgs`**. |
| **`pipeline.h`** | **`PipelineArgs`**, **`InsertionArgs`**, **`run_pipeline`**, **`run_insertion_pipeline`**, **`dump_state`**. |
| **`pipeline.cpp`** | Orchestrates Stages **0–9**; builds **`TreeState`**; calls **`save_tree_state`**; **`run_insertion_pipeline`** wraps **`insert_batch`**. |
| **`types.h`** | Global constants (**`MAX_DIST`**, **`MIN_BRANCH`**, **`EPS`**, …), Eigen aliases, **Adjacency** / **EdgeWeights**, result structs (**`ClnjResult`**, **`TreeStats`**, …), **`DistModel`**, **`TreeAlgo`**, **`Cluster`**, **`ClusterTree`**, **`TreeState`**, **`InsertionResult`**. |
| **`fasta_parser.h` / `.cpp`** | **`parse_fasta`**, **`load_clean_fasta`** — I/O and filtering. |
| **`alignment_utils.h` / `.cpp`** | **`is_aligned`**, **`align_with_mafft`**. |
| **`one_hot.h` / `.cpp`** | **`OneHotData`**, **`one_hot_encode`** — (n,m,4) layout for batch distance + diagnostics. |
| **`distance.h` / `.cpp`** | Pair statistics, **AUTO**, **`compute_distance_matrix`**, model corrections, **`compute_distance_to_subset`** / **`_cached`**, **`encode_sequence_flat`** for insertion. |
| **`iqtree_interface.h` / `.cpp`** | External **IQ-TREE2** invocation and log parsing. |
| **`kmer_distance.h` / `.cpp`** | **Experimental** MinHash / Mash-style distances (`kmer_distance_matrix`). |
| **`algorithm2.h` / `.cpp`** | **`construct_FC_and_GU`**. |
| **`mlvmst.h` / `.cpp`** | **`build_mlvmst`**. |
| **`distance_oracle.h` / `.cpp`** | **`DistanceOracle`**, **`get_tree_distance`** — matrix vs tree-path queries for local NJ. |
| **`nj.h` / `.cpp`** | Standard NJ (global driver + **local** neighborhood rebuild). |
| **`mfnj.h` / `.cpp`** | MFNJ polytomy-aware local/global NJ. |
| **`bionj.h` / `.cpp`** | BioNJ local/global. |
| **`clnj.h` / `.cpp`** | **`clnj_clean`** — main CLNJ loop over MLVMST internal nodes. |
| **`mdl_clustering.h` / `.cpp`** | **`precompute_tree_distances`**, **`mdl_cluster_tree`** — Stage 9. |
| **`tree_state.h` / `.cpp`** | **`save_tree_state`**, **`load_tree_state`** — binary persistence. |
| **`online_insertion.h` / `.cpp`** | **`insert_taxon`**, **`insert_batch`** — Phase 2. |
| **`tree_analysis.h` / `.cpp`** | **`analyze_tree`** — Stage 7 reporting. |

---

## Implementation note: sequences, one-hot, and planned direct encoding

| Layer | Today | Planned refactor (your direction) |
|-------|-------|-----------------------------------|
| **Biological input** | Aligned strings in **`CleanResult`** | Unchanged. |
| **Stage 3 (batch)** | **`one_hot_encode` → `compute_distance_matrix`** | Same **formulas**, but iterate **columns** on **A/C/G/T** (+ masks) **without** allocating full **(n,m,4)** float tensor. |
| **Stage 2 diagnostics** | **`compute_pair_counts(OneHotData)`** | Equivalent stats from **character** scans. |
| **Insertion** | Already uses **`encode_sequence_flat`** + **`compute_distance_to_subset_cached`** | Aligns with “direct” per-site encoding; batch path may follow the same style. |

Algorithm 2, MLVMST, CLNJ, MDL, and insertion **logic** depend on **D** and **tree topology**, not on the internal representation used to fill **D**.

---

## Experimental: alignment-free k-mer distance

**Files:** `kmer_distance.{h,cpp}`, branches in **`pipeline.cpp`**.

- **`--distance-method kmer`**: skips alignment Stage 0 shortcut, uses **`kmer_distance_matrix`** (canonical k-mers, MinHash, Mash distance).
- Included for **experimentation** / comparison, **not** the primary methodology for the **alignment-based** LTGM + CLNJ story.
- **Dynamic insertion** is designed around **alignment-based** **`TreeState`** (cached nuc indices, **`aln_len`**); k-mer **build** + **insert** combinations should be validated carefully before use in production.

---

## Distance models

| Model | Description |
|-------|-------------|
| **JC69** | Jukes–Cantor |
| **K2P** | Kimura 2-parameter (Ts vs Tv) |
| **TN93** | Tamura–Nei + unequal base frequencies |
| **LOGDET** | Log-determinant / paralinear |
| **AUTO** | Heuristic from alignment diagnostics |
| **IQTREE** | External IQ-TREE2 → mapped to JC69/K2P/TN93 + optional gamma |

Gamma: **`--gamma-alpha`** applies the same transform as in the Python pipeline.

---

## Tree reconstruction algorithms

| CLI | Algorithm |
|-----|-----------|
| **`nj`** | Neighbor-joining |
| **`mfnj`** (default) | Multi-furcating NJ (tie-aware polytomies) |
| **`bionj`** | BioNJ (variance-style weighting) |

Used in **CLNJ** (Stage 6) and **Phase 2** local rebuilds (`state.tree_algo`).

---

## Key features

- **Full lifecycle**: batch **D → Algorithm 2 → MLVMST → CLNJ** + **MDL cluster tree** + **serialized state** + **online insertion**.
- **Multi-model** distances with **pairwise deletion** and **gamma** option.
- **MAFFT** integration; **IQ-TREE2** optional model selection.
- **Distance oracle** for consistent local distances with **hidden** nodes.
- **MDL**-guided **cluster assignment** for scalable incremental updates.
- **Deterministic subsampling** (seed 42).
- **Eigen** + **CLI11** via CMake **FetchContent**.

---

## Requirements and building

**Build:** C++17, CMake ≥ 3.14, first configure downloads **Eigen 3.4.0** and **CLI11 2.4.1**.

**Runtime tools (optional):** **MAFFT** (auto-align), **IQ-TREE2** (`--model IQTREE`).

```bash
cd clnj_cpp
mkdir -p build && cd build
cmake ..
cmake --build . -j"$(nproc)"
```

Executable: **`clnj_pipeline`**.

---

## Examples

**1. Full pipeline: build with state, then insert**

```bash
./clnj_pipeline build data/aligned.fasta 100 0 mfnj \
  --model TN93 \
  --save-state project.state \
  --mdl-mj 7.0 \
  --mdl-merge-threshold 100

./clnj_pipeline insert project.state new_taxa.fasta -o project_updated.state
# Further waves:
./clnj_pipeline insert project_updated.state more_taxa.fasta -o project_v3.state
```

**2. Batch-only tree (no insertion / no Stage 9)**

```bash
./clnj_pipeline build data/aligned.fasta --model AUTO
```

**3. Inspect a state file**

```bash
./clnj_pipeline dump project.state
```

---

## Diagnostics

- **`--report-negatives`**: aggregate negative intermediate values during local NJ (non-additive **D**).
- **`--trace-zero-edges`**: log near-zero edges during CLNJ (`initial_pair`, `mfnj_polytomy`, …).

---

## Pipeline constants

| Symbol | Value | Role |
|--------|-------|------|
| `MAX_DIST` | 10.0 | Distance saturation cap |
| `MIN_BRANCH` | 1e-6 | Default branch floor |
| `EPS` | 1e-12 | Floating-point tolerance |
| `TIE_TOL` | 1e-9 | MFNJ tie tolerance |
| `SEED` | 42 | Subsampling RNG |
| `LOG_FLOOR` | 1e-300 | Numerical floor in logs |
| `FREQ_FLOOR` | 1e-10 | Base frequency floor |

---

## GitHub (`MSTree_Phylogeny`)

**URL:** https://github.com/Jatin07gupta/MSTree_Phylogeny

**Push reminder:** commits are local until **`git push`** succeeds (PAT or SSH). To create the remote repo and push in one step (with a classic PAT):

```bash
export GITHUB_TOKEN=ghp_your_token_here
./scripts/create_github_repo_and_push.sh
```

See script comments in **`scripts/create_github_repo_and_push.sh`** for security notes.

---

## References

- Saitou & Nei (1987) — Neighbor-joining. *MBE*.
- Fernandez et al. (2023) — Multi-furcating NJ. *J Mol Evol*.
- Gascuel (1997) — BioNJ. *MBE*.
- Jukes & Cantor (1969); Kimura (1980); Tamura & Nei (1993) — Substitution models.
- Katoh & Standley (2013) — MAFFT. *MBE*.
- Nguyen et al. (2015) — IQ-TREE. *MBE*.

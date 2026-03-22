# MDL Clustering Dynamic Insertion — Implementation Plan

## Step 1: Codebase Analysis (clnj_cpp)

### What Exists
- **Static pipeline**: Stages 0–7 (alignment → load/clean → one-hot → diagnostics → distance → Algorithm 2 → MLVMST → CLNJ → analysis)
- **Core modules**: fasta_parser, one_hot, distance, algorithm2, mlvmst, distance_oracle, nj, mfnj, bionj, clnj, tree_analysis
- **types.h**: CleanResult, PairStats, Algo2Result, MlvmstResult, ClnjResult, TreeStats, Adjacency, EdgeWeights

### What Is Missing
1. **Cluster structures**: Cluster, ClusterTree, InsertionResult, TreeState
2. **mdl_clustering**: precompute_tree_distances, mdl_cluster_tree, TreeDistCache
3. **tree_state**: save_tree_state, load_tree_state (binary serialization)
4. **online_insertion**: insert_taxon, insert_batch (Steiner tree, boundary edges, NJ rebuild, hierarchical assignment)
5. **distance**: compute_distance_to_subset, compute_distance_to_subset_cached, encode_sequence_flat
6. **pipeline**: Stage 9 (MDL + save), run_insertion_pipeline, save_state_path, mdl args
7. **main**: insert subcommand

---

## Step 2: GPU Reference (clnj_cuda)

### Dynamic Insertion Flow (Algorithm 4)
1. **Phase 1**: Compute distances from each new seq to cluster centers only (CPU: compute_distance_to_subset / compute_distance_to_subset_cached)
2. **Phase 2**: find_cluster_hierarchical — descend cluster tree, pick nearest child by center distance
3. **Phase 3**: remap_hidden_nodes(state, old_n, batch_size)
4. **Phase 4**: Expand D; compute new→old and new→new distances per cluster (CPU path)
5. **Phase 5**: For each affected cluster: Steiner tree → boundary edges → remove Steiner → NJ on all members → install_nj_and_reconnect

### MDL Clustering (Static, used at build time)
- precompute_tree_distances: n BFS runs, n×n cache
- recursive_cluster: find center, DL_unsplit, find_best_split_node, split_at_node, DL_split; if DL_split < DL_unsplit recurse
- merge_small_clusters: merge leaf clusters < threshold if DL improves

### CPU-Only Translation
- Remove: gpu_dispatch.h, cuda_available(), gpu_compute_distance_matrix_from_sequences
- Always use: compute_distance_to_subset, compute_distance_to_subset_cached
- Keep: use_gpu parameter in API but ignore (or remove for CPU-only)

---

## Step 3: Implementation Plan

### 3.1 types.h
Add:
- Cluster, ClusterTree (MDL structures)
- InsertionResult
- TreeState (with cached_nuc_idx, cached_valid)

### 3.2 distance.h / distance.cpp
Add:
- jc69_correct, k2p_correct, tn93_correct (expose if static)
- compute_distance_to_subset(OneHotData, new_sequence, target_indices, ...)
- compute_distance_to_subset_cached(nuc_idx, valid, ..., new_nuc, new_valid, target_indices, ...)
- encode_sequence_flat(seq, m, out_nuc, out_valid)
- char_to_nuc (static helper)

### 3.3 mdl_clustering.h / mdl_clustering.cpp
New files. Copy from clnj_cuda, no GPU deps. Uses:
- distance_oracle.h: get_tree_distance
- types.h: Cluster, ClusterTree, TreeDistCache

### 3.4 tree_state.h / tree_state.cpp
New files. Copy from clnj_cuda. Pure binary I/O.

### 3.5 online_insertion.h / online_insertion.cpp
New files. Copy from clnj_cuda, strip GPU:
- Remove #include "gpu_dispatch.h"
- Remove gpu_active, cuda_available(), gpu_compute_distance_matrix_from_sequences
- Always use CPU path (using_cache ? compute_distance_to_subset_cached : compute_distance_to_subset)
- Remove use_gpu parameter from insert_taxon/insert_batch, or keep and ignore

### 3.6 pipeline.h / pipeline.cpp
- PipelineArgs: add save_state_path, mdl_mj, mdl_merge_threshold
- Add InsertionArgs, run_insertion_pipeline
- After Stage 7: if save_state_path non-empty → Stage 9 (mdl_cluster_tree, build TreeState, save_tree_state)
- Implement run_insertion_pipeline: load_tree_state, parse_fasta, insert_batch, save_tree_state

### 3.7 main.cpp
- Add subcommands: build, insert
- build: existing logic + --save-state, --mdl-mj, --mdl-merge-threshold
- insert: state, fasta, --output-state, --verbose (no --gpu)

### 3.8 CMakeLists.txt
Add: src/mdl_clustering.cpp, src/tree_state.cpp, src/online_insertion.cpp

---

## Step 4: Dynamic Insertion Trigger
- User runs `insert` subcommand with state path and new FASTA
- run_insertion_pipeline loads state, parses new FASTA, calls insert_batch

## Step 5: Cluster Evaluation Logic
- find_cluster_hierarchical: at each level, pick child with minimum distance from new seq to child center
- Leaf cluster → assignment

## Step 6: MDL Score Computation
- compute_DL: n*log(variance) + mj*log(n) + n*(log(2π)+1)
- compute_DL_split: sum over parts
- Split if DL_split < DL_unsplit

## Step 7: Structure Update Procedure
- remap_hidden_nodes before adding new taxa
- remove_steiner_tree, install_nj_and_reconnect
- cluster_ptr->observed_members = all_members
- node_to_cluster[new_id] = cid

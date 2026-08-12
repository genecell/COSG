"""
GPU-accelerated core computation for COSG.

Replicates the exact COSG scoring formula using CuPy + cuSPARSE:

    cos_sim(g, k) = cosine_similarity(expr_g, λ_k)
    e²(g, k) = cos_sim(g, k)²
    Σe²(g) = Σ_k e²(g, k)    [row-wise sum across all groups]
    
    For sparse input:
        score(g, k) = e²(g,k) / ((1-μ)·e²(g,k) + μ·Σe²(g)) × cos_sim(g, k)
        When μ=1:  score(g, k) = cos_sim(g, k)³ / Σe²(g)
    
    For dense input:
        Same formula applied only at nonzero positions.

GPU strategy:
    1. Transfer expression matrix and indicator matrix to GPU
    2. Compute cosine similarity via SpMM + L2 normalization (cuSPARSE)
    3. Apply the COSG penalty formula element-wise on GPU
    4. Return genexlambda as dense numpy array (genes × groups)

All intermediate results stay in GPU memory until final transfer.
"""

import numpy as np
import time
from scipy import sparse as sp_sparse

import cupy as cp
import cupyx.scipy.sparse as cp_sparse


def _cosine_similarity_gpu(X_gpu, Y_gpu):
    """
    GPU cosine similarity between rows of X and rows of Y.
    
    Equivalent to sklearn.metrics.pairwise.cosine_similarity(X, Y, dense_output=True).
    
    Parameters
    ----------
    X_gpu : cupy sparse CSR or dense, shape (n_samples_X, n_features)
    Y_gpu : cupy sparse CSR or dense, shape (n_samples_Y, n_features)
    
    Returns
    -------
    similarity : cupy dense array, shape (n_samples_X, n_samples_Y)
    """
    # Compute L2 norms for X
    if cp_sparse.issparse(X_gpu):
        Xsq = X_gpu.copy()
        Xsq.data = Xsq.data ** 2
        X_norms = cp.sqrt(cp.array(Xsq.sum(axis=1)).ravel())
        del Xsq
    else:
        X_norms = cp.sqrt(cp.sum(X_gpu ** 2, axis=1))
    
    # Compute L2 norms for Y
    if cp_sparse.issparse(Y_gpu):
        Ysq = Y_gpu.copy()
        Ysq.data = Ysq.data ** 2
        Y_norms = cp.sqrt(cp.array(Ysq.sum(axis=1)).ravel())
        del Ysq
    else:
        Y_norms = cp.sqrt(cp.sum(Y_gpu ** 2, axis=1))
    
    # Avoid division by zero
    X_norms = cp.maximum(X_norms, 1e-12)
    Y_norms = cp.maximum(Y_norms, 1e-12)
    
    # Dot product: X @ Y^T
    if cp_sparse.issparse(Y_gpu):
        Yt = Y_gpu.T.tocsr()
    else:
        Yt = Y_gpu.T
    
    dots = X_gpu @ Yt
    
    # Convert to dense for normalization
    if cp_sparse.issparse(dots):
        dots = dots.toarray()
    
    # Normalize: divide rows by X_norms, columns by Y_norms
    dots /= X_norms[:, cp.newaxis]
    dots /= Y_norms[cp.newaxis, :]
    
    return dots


def _scipy_to_cupy_csr(scipy_csr):
    """Convert scipy CSR matrix to CuPy CSR matrix."""
    scipy_csr = scipy_csr.tocsr()
    if scipy_csr.dtype != np.float32:
        scipy_csr = scipy_csr.astype(np.float32)
    return cp_sparse.csr_matrix(
        (cp.array(scipy_csr.data),
         cp.array(scipy_csr.indices),
         cp.array(scipy_csr.indptr)),
        shape=scipy_csr.shape
    )


def _apply_cosg_penalty_gpu(cosine_sim, mu):
    """
    Apply the COSG penalty formula on GPU.
    
    Matches the exact formula from cosg.py lines 381-415.
    
    Parameters
    ----------
    cosine_sim : cupy dense array, shape (n_genes, n_groups)
    mu : float
    
    Returns
    -------
    genexlambda : cupy dense array, shape (n_genes, n_groups)
    """
    # e_power2 = cos_sim²
    e_power2 = cp.multiply(cosine_sim, cosine_sim)
    
    # e_power2_sum = row-wise sum of cos_sim² (sum across groups for each gene)
    e_power2_sum = cp.sum(e_power2, axis=1)  # shape: (n_genes,)
    
    # Mask of nonzero positions
    pos_nonzero = cosine_sim != 0
    
    if mu == 1:
        # When mu=1: score = cos_sim² / Σcos_sim² × cos_sim = cos_sim³ / Σcos_sim²
        # Matches: genexlambda.data = genexlambda.data / np.repeat(e_power2_sum, ...)
        safe_sum = cp.maximum(e_power2_sum, 1e-30)
        genexlambda = e_power2 / safe_sum[:, cp.newaxis]
        genexlambda = cp.where(pos_nonzero, cp.multiply(genexlambda, cosine_sim), 0.0)
    else:
        # General mu: score = cos_sim² / ((1-mu)*cos_sim² + mu*Σcos_sim²) × cos_sim
        denom = (1 - mu) * e_power2 + mu * e_power2_sum[:, cp.newaxis]
        safe_denom = cp.where(pos_nonzero, denom, 1.0)
        genexlambda = cp.where(pos_nonzero, e_power2 / safe_denom, 0.0)
        genexlambda = cp.where(pos_nonzero, cp.multiply(genexlambda, cosine_sim), 0.0)
    
    return genexlambda


def gpu_core_cosg(cellxgene, cluster_mat, mu, is_sparse_input,
                  remove_lowly_expressed=False, expressed_pct=0.1,
                  expressed_min_num_cells_in_target_group=3, verbosity=0):
    """
    GPU-accelerated COSG core computation (non-batch mode).

    Drop-in replacement for the cosine_similarity + penalty formula block
    (lines ~366-415) in the original cosg.py.

    Parameters
    ----------
    cellxgene : scipy.sparse.csr_matrix or np.ndarray
        Expression matrix, shape (n_cells, n_genes).
    cluster_mat : scipy.sparse.csr_matrix
        Group indicator matrix, shape (n_groups, n_cells).
    mu : float
        Penalty parameter (non-negative).
    is_sparse_input : bool
        Whether cellxgene is sparse.
    remove_lowly_expressed : bool
        Whether to filter lowly expressed genes.
    expressed_pct : float
        Minimum fraction of cells expressing a gene.
    expressed_min_num_cells_in_target_group : int
        Minimum absolute number of cells expressing a gene.
    verbosity : int
        Verbosity level.

    Returns
    -------
    genexlambda : np.ndarray
        COSG score matrix, shape (n_genes, n_groups).
        Always returns dense ndarray for consistent downstream processing.
    """
    t0 = time.time()
    
    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]
    
    if verbosity > 0:
        print(f"[gpu_cosg] Transferring {n_cells:,} × {n_genes:,} matrix to GPU...")
    
    # ── Transfer to GPU ──────────────────────────────────────────────────
    t_xfer = time.time()
    
    if is_sparse_input:
        X_gpu = _scipy_to_cupy_csr(cellxgene)
    else:
        X_gpu = cp.array(np.asarray(cellxgene, dtype=np.float32))
    
    cm_gpu = _scipy_to_cupy_csr(cluster_mat)
    
    if verbosity > 0:
        print(f"[gpu_cosg] GPU transfer: {time.time() - t_xfer:.3f}s")
    
    # ── Cosine similarity: (genes × groups) ──────────────────────────────
    # We want cosine_similarity(X.T, cluster_mat) → (n_genes, n_groups)
    # But X.T.tocsr() causes OOM for large matrices.
    # Instead compute directly: dots = X.T @ cm.T = (cm @ X).T
    # cm_gpu is (n_groups, n_cells), X_gpu is (n_cells, n_genes)
    # cm_gpu @ X_gpu → (n_groups, n_genes) — efficient SpMM, no big transpose
    t_cos = time.time()

    # Convert small cluster_mat to dense for efficient SpMM (csrmm)
    # cm_gpu is (n_groups, n_cells) — small number of rows, safe to densify
    cm_dense = cm_gpu.toarray()  # (n_groups, n_cells)
    # X_gpu.T gives CSC view (no copy), then dense @ CSC works via BLAS
    # But we want: cm_dense @ X_gpu → (n_groups, n_genes)
    # For sparse X_gpu: use X_gpu.T (CSC) @ cm_dense.T → (n_genes, n_groups)
    if cp_sparse.issparse(X_gpu):
        # X_gpu is CSR (cells, genes). X_gpu.T is CSC (genes, cells).
        # We want (genes, groups) = X.T @ cm.T
        # cuSPARSE csrmm: CSR @ dense → dense
        # So: (cm_gpu @ X_gpu) where cm_gpu is now dense won't help with sparse X.
        # Instead: X_gpu.T.tocsr() is expensive. Use: result = (cm_dense @ X_gpu).T
        # CuPy supports dense @ sparse → dense via __matmul__
        dots = cm_dense @ X_gpu  # dense(n_groups,n_cells) @ sparse(n_cells,n_genes) → dense
    else:
        dots = cm_dense @ X_gpu
    dots = dots.T  # → (n_genes, n_groups), cheap for dense

    # Keep references for expressed_pct filtering later
    cm_dense_ref = cm_dense if remove_lowly_expressed else None
    if not remove_lowly_expressed:
        del cm_dense

    # L2 norms for X.T rows (= column norms of X)
    if cp_sparse.issparse(X_gpu):
        Xsq = X_gpu.copy()
        Xsq.data = Xsq.data ** 2
        gene_norms = cp.sqrt(cp.array(Xsq.sum(axis=0)).ravel())  # (n_genes,)
        del Xsq
    else:
        gene_norms = cp.sqrt(cp.sum(X_gpu ** 2, axis=0))  # (n_genes,)

    # L2 norms for cluster_mat rows
    cm_sq = cm_gpu.copy()
    cm_sq.data = cm_sq.data ** 2
    group_norms = cp.sqrt(cp.array(cm_sq.sum(axis=1)).ravel())  # (n_groups,)
    del cm_sq

    gene_norms = cp.maximum(gene_norms, 1e-12)
    group_norms = cp.maximum(group_norms, 1e-12)

    # Normalize: rows by gene norms, columns by group norms
    dots /= gene_norms[:, cp.newaxis]
    dots /= group_norms[cp.newaxis, :]

    cosine_sim = dots
    # Keep X_gpu reference for expressed_pct filtering
    X_gpu_ref = X_gpu if remove_lowly_expressed else None
    del dots
    if not remove_lowly_expressed:
        del X_gpu
    del cm_gpu
    
    if verbosity > 0:
        print(f"[gpu_cosg] Cosine similarity: {time.time() - t_cos:.3f}s")
    
    # ── COSG penalty formula ─────────────────────────────────────────────
    t_pen = time.time()
    
    genexlambda = _apply_cosg_penalty_gpu(cosine_sim, mu)

    del cosine_sim

    if verbosity > 0:
        print(f"[gpu_cosg] Penalty formula: {time.time() - t_pen:.3f}s")

    # ── Expressed-pct filtering (GPU) ────────────────────────────────────
    if remove_lowly_expressed:
        t_filt = time.time()

        if is_sparse_input:
            X_binary = X_gpu_ref.copy()
            X_binary.data = cp.ones_like(X_binary.data)
        else:
            X_binary = (X_gpu_ref != 0).astype(cp.float32)

        # nnz per gene per group: cm_dense_ref @ X_binary → (groups, genes)
        nnz_per_gene_per_group = (cm_dense_ref @ X_binary).T  # → (genes, groups)
        del X_binary

        # Group sizes and thresholds
        n_cells_per_group = cm_dense_ref.sum(axis=1)  # (groups,)
        thresholds = cp.maximum(
            n_cells_per_group * expressed_pct,
            expressed_min_num_cells_in_target_group
        )

        # Set lowly expressed to -1 (matching CPU behavior)
        lowly_expressed_mask = nnz_per_gene_per_group < thresholds[cp.newaxis, :]
        genexlambda[lowly_expressed_mask] = -1.0
        del nnz_per_gene_per_group, lowly_expressed_mask

        if verbosity > 0:
            print(f"[gpu_cosg] Expressed filtering: {time.time() - t_filt:.3f}s")

    del X_gpu_ref, cm_dense_ref

    # ── Transfer back to CPU ─────────────────────────────────────────────
    t_back = time.time()
    result = cp.asnumpy(genexlambda)

    del genexlambda
    cp.get_default_memory_pool().free_all_blocks()

    if verbosity > 0:
        print(f"[gpu_cosg] GPU → CPU transfer: {time.time() - t_back:.3f}s")
        print(f"[gpu_cosg] Total GPU computation: {time.time() - t0:.3f}s")

    return result


def gpu_core_cosg_batched(cellxgene, cluster_mat, mu, batch_info, unique_batches,
                          batch_cell_number_threshold=3, is_sparse_input=True,
                          verbosity=0):
    """
    GPU-accelerated COSG with batch correction.
    
    Computes cosine similarities separately per batch on GPU, then averages.
    
    Parameters
    ----------
    cellxgene : scipy sparse or np.ndarray, (n_cells, n_genes)
    cluster_mat : scipy sparse, (n_groups, n_cells)
    mu : float
    batch_info : np.ndarray of batch labels per cell
    unique_batches : np.ndarray of unique batch labels
    batch_cell_number_threshold : int
    is_sparse_input : bool
    verbosity : int
    
    Returns
    -------
    genexlambda : np.ndarray, (n_genes, n_groups)
    """
    t0 = time.time()
    
    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]
    
    if verbosity > 0:
        print(f"[gpu_cosg] Batch mode: {len(unique_batches)} batches, "
              f"{n_cells:,} cells × {n_genes:,} genes")
    
    # Transfer full expression matrix to GPU once
    if is_sparse_input:
        X_gpu = _scipy_to_cupy_csr(cellxgene)
    else:
        X_gpu = cp.array(np.asarray(cellxgene, dtype=np.float32))
    
    # Keep cluster_mat on CPU for slicing (scipy is better at fancy indexing)
    cm_csr = cluster_mat.tocsr()
    
    accum = cp.zeros((n_genes, n_groups), dtype=cp.float32)
    valid_counts = cp.zeros(n_groups, dtype=cp.float32)
    
    for batch in unique_batches:
        indices = np.flatnonzero(batch_info == batch)
        if indices.size == 0:
            continue
        
        # Cluster counts per group for this batch (CPU side, cheap)
        cluster_counts = np.array(cm_csr[:, indices].sum(axis=1)).flatten()
        valid_mask = (cluster_counts >= batch_cell_number_threshold).astype(np.float32)
        if not np.any(valid_mask):
            continue
        
        # Slice expression matrix for this batch
        if is_sparse_input:
            X_batch_cpu = cellxgene[indices]
            X_batch_gpu = _scipy_to_cupy_csr(X_batch_cpu)
        else:
            X_batch_gpu = cp.array(np.asarray(cellxgene[indices], dtype=np.float32))

        # Slice indicator matrix for this batch
        cm_batch = cm_csr[:, indices].tocsr()
        cm_batch_gpu = _scipy_to_cupy_csr(cm_batch)

        # Cosine similarity via memory-efficient path: cm @ X → (groups, genes)
        dots = cm_batch_gpu @ X_batch_gpu  # (n_groups, n_genes_batch)
        if cp_sparse.issparse(dots):
            dots = dots.toarray()
        dots = dots.T  # (n_genes, n_groups)

        # L2 norms
        if cp_sparse.issparse(X_batch_gpu):
            Xsq = X_batch_gpu.copy()
            Xsq.data = Xsq.data ** 2
            gene_norms = cp.sqrt(cp.array(Xsq.sum(axis=0)).ravel())
            del Xsq
        else:
            gene_norms = cp.sqrt(cp.sum(X_batch_gpu ** 2, axis=0))
        cm_sq = cm_batch_gpu.copy()
        cm_sq.data = cm_sq.data ** 2
        group_norms = cp.sqrt(cp.array(cm_sq.sum(axis=1)).ravel())
        del cm_sq
        gene_norms = cp.maximum(gene_norms, 1e-12)
        group_norms = cp.maximum(group_norms, 1e-12)
        dots /= gene_norms[:, cp.newaxis]
        dots /= group_norms[cp.newaxis, :]
        sim_batch = dots

        # Mask invalid clusters
        valid_mask_gpu = cp.array(valid_mask, dtype=cp.float32)
        sim_batch *= valid_mask_gpu[cp.newaxis, :]

        accum += sim_batch
        valid_counts += valid_mask_gpu

        del X_batch_gpu, cm_batch_gpu, sim_batch, dots
    
    # Average across batches
    safe_counts = cp.maximum(valid_counts, 1e-30)
    cosine_sim = accum / safe_counts[cp.newaxis, :]
    
    # Where valid_counts == 0, set to NaN (matching CPU behavior)
    zero_mask = valid_counts == 0
    if cp.any(zero_mask):
        cosine_sim[:, zero_mask] = cp.nan
    
    del accum, X_gpu
    
    # Apply COSG penalty formula
    # Replace NaN with 0 for computation (NaN groups get 0 scores)
    cosine_sim_clean = cp.nan_to_num(cosine_sim, nan=0.0)
    genexlambda = _apply_cosg_penalty_gpu(cosine_sim_clean, mu)
    
    result = cp.asnumpy(genexlambda)
    
    del cosine_sim, cosine_sim_clean, genexlambda
    cp.get_default_memory_pool().free_all_blocks()
    
    if verbosity > 0:
        print(f"[gpu_cosg] Batch total: {time.time() - t0:.3f}s")

    return result


def _auto_chunk_size(n_genes, is_sparse, cellxgene=None,
                     target_fraction=0.5, min_chunk=1000, max_chunk=100_000):
    """
    Auto-select chunk size based on available GPU memory.

    Strategy: Use at most `target_fraction` of free GPU memory for the chunk.
    The chunk memory is approximately:
      - Sparse: nnz_per_cell * chunk_size * 8 bytes (data + indices) + overhead
      - Dense: chunk_size * n_genes * 4 bytes
    Plus the transposed version (similar size).

    Parameters
    ----------
    n_genes : int
    is_sparse : bool
    cellxgene : sparse matrix or None (used to estimate density)
    target_fraction : float
        Fraction of free GPU memory to use for the chunk (default 0.5)
    min_chunk : int
        Minimum chunk size
    max_chunk : int
        Maximum chunk size

    Returns
    -------
    chunk_size : int
    """
    mem_free, mem_total = cp.cuda.Device().mem_info
    available = mem_free * target_fraction

    if is_sparse and cellxgene is not None:
        # Estimate average nnz per cell
        avg_nnz_per_cell = cellxgene.nnz / cellxgene.shape[0]
        # Each nonzero needs: 4 bytes (data) + 4 bytes (index) in CSR
        # We need both original and transposed -> x2
        # Plus indicator matrix overhead (negligible)
        # Plus intermediate results during SpMM
        bytes_per_cell = avg_nnz_per_cell * 8 * 2.5  # 2.5x safety factor
    else:
        # Dense: chunk x genes x 4 bytes, x2 for transpose, x1.5 safety
        bytes_per_cell = n_genes * 4 * 3.0

    chunk_size = int(available / bytes_per_cell)
    chunk_size = max(min_chunk, min(chunk_size, max_chunk))

    return chunk_size


def _should_use_chunked(n_cells, n_genes, is_sparse, cellxgene=None):
    """
    Decide whether to use chunked streaming or monolithic transfer.

    For small datasets, monolithic is faster (no per-chunk overhead).
    For large datasets that don't fit in GPU memory, chunked is necessary.

    Heuristic: Use chunked if estimated GPU memory > 70% of free memory.
    The monolithic path is faster when the data fits in GPU memory.
    """
    mem_free, _ = cp.cuda.Device().mem_info

    # Estimate full-matrix GPU memory
    if is_sparse and cellxgene is not None:
        # Sparse: data + indices + indptr + working memory for SpMM
        est_bytes = (cellxgene.data.nbytes + cellxgene.indices.nbytes +
                     cellxgene.indptr.nbytes) * 2.5
    else:
        est_bytes = n_cells * n_genes * 4 * 2.5

    # Use chunked only when matrix won't fit comfortably in GPU memory
    if est_bytes > mem_free * 0.7:
        return True
    return False


def gpu_core_cosg_chunked(cellxgene, cluster_mat, mu, is_sparse_input,
                           remove_lowly_expressed=True, expressed_pct=0.1,
                           expressed_min_num_cells_in_target_group=3,
                           chunk_size=None, verbosity=0):
    """
    GPU-accelerated COSG with chunked cell-axis streaming.

    Instead of transferring the entire expression matrix to GPU at once,
    processes cells in chunks and accumulates partial results. This:
      - Reduces GPU memory from O(n_cells * n_genes) to O(chunk_size * n_genes)
      - Can improve speed by avoiding giant CuPy CSR construction overhead
      - Enables processing of arbitrarily large datasets on any GPU

    Mathematical basis:
        X.T @ L = Sum_chunk (X_chunk.T @ L_chunk)  [dot products decompose]
        ||gene_g||^2 = Sum_chunk ||gene_g in chunk||^2  [norms decompose]
        nnz(g,k) = Sum_chunk nnz(g,k in chunk)  [expressed pct decomposes]

    Parameters
    ----------
    cellxgene : scipy.sparse.csr_matrix or np.ndarray
        Expression matrix, shape (n_cells, n_genes).
    cluster_mat : scipy.sparse.csr_matrix
        Group indicator matrix, shape (n_groups, n_cells).
    mu : float
        COSG penalty parameter.
    is_sparse_input : bool
        Whether cellxgene is sparse.
    remove_lowly_expressed : bool
        Whether to filter genes with low expression percentage.
    expressed_pct : float
        Minimum fraction of cells in a group that must express a gene.
    expressed_min_num_cells_in_target_group : int
        Minimum number of cells expressing a gene in target group.
    chunk_size : int or None
        Number of cells per chunk. None = auto-select based on GPU memory.
    verbosity : int
        0 = silent, 1 = progress, 2 = detailed timing.

    Returns
    -------
    genexlambda : np.ndarray, shape (n_genes, n_groups)
        COSG penalty scores. Dense numpy array.
    """
    t0 = time.time()

    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]

    # -- Auto-select chunk size --
    if chunk_size is None:
        chunk_size = _auto_chunk_size(n_genes, is_sparse_input, cellxgene)
    chunk_size = min(chunk_size, n_cells)  # don't exceed dataset
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    if verbosity > 0:
        print(f"[gpu_cosg_chunked] {n_cells:,} cells x {n_genes:,} genes, "
              f"{n_groups} groups")
        print(f"[gpu_cosg_chunked] Chunk size: {chunk_size:,} cells, "
              f"{n_chunks} chunks")

    # -- Prepare cell-to-group mapping (CPU, tiny) --
    # cluster_mat is (n_groups x n_cells), row k has 1s for cells in group k
    # Extract from cluster_mat: for each cell, which group index?
    cm_csc = cluster_mat.tocsc()
    cell_group_indices = np.array(cm_csc.indices, dtype=np.int32)

    # Group sizes (for group norms)
    group_sizes_np = np.array(cluster_mat.sum(axis=1)).ravel().astype(np.float32)
    group_sizes_gpu = cp.array(group_sizes_np)
    group_norms_gpu = cp.sqrt(group_sizes_gpu)  # (n_groups,)

    # -- Initialize accumulators on GPU --
    dot_accum = cp.zeros((n_genes, n_groups), dtype=cp.float32)
    gene_sq_accum = cp.zeros(n_genes, dtype=cp.float32)

    if remove_lowly_expressed and expressed_pct > 0:
        nnz_accum = cp.zeros((n_genes, n_groups), dtype=cp.float32)
    else:
        nnz_accum = None

    # -- Ensure CSR format on CPU --
    if is_sparse_input:
        cellxgene_csr = cellxgene.tocsr()
        if cellxgene_csr.dtype != np.float32:
            cellxgene_csr = cellxgene_csr.astype(np.float32)
    else:
        cellxgene_np = np.asarray(cellxgene, dtype=np.float32)

    t_transfer_total = 0.0
    t_compute_total = 0.0

    # -- Stream chunks --
    # Strategy: transfer X_chunk as (cells x genes) CSR directly (no CPU transpose),
    # then compute Lam_dense.T @ X_chunk = (groups x cells) @ (cells x genes) -> (groups x genes)
    # This avoids the expensive CPU-side .T.tocsr() which dominates runtime.
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        actual_chunk = end - start

        if verbosity > 1:
            print(f"  Chunk {chunk_idx+1}/{n_chunks}: cells [{start:,}:{end:,}]",
                  end="", flush=True)

        # -- Slice chunk on CPU and transfer to GPU --
        t_xfer = time.time()

        if is_sparse_input:
            X_chunk_cpu = cellxgene_csr[start:end]  # (chunk x genes), CSR — O(1) slice
            # Transfer directly to GPU without transposing
            X_chunk_gpu = cp_sparse.csr_matrix(
                (cp.array(X_chunk_cpu.data),
                 cp.array(X_chunk_cpu.indices),
                 cp.array(X_chunk_cpu.indptr)),
                shape=X_chunk_cpu.shape
            )
        else:
            X_chunk_np = cellxgene_np[start:end]  # (chunk x genes)
            X_chunk_gpu = cp.array(X_chunk_np)  # keep as (chunk x genes), dense

        # Build dense indicator for this chunk: (groups x chunk)
        # Lam_dense_T[k, i] = 1 if cell i belongs to group k
        chunk_groups = cell_group_indices[start:end]
        Lam_dense_T = cp.zeros((n_groups, actual_chunk), dtype=cp.float32)
        Lam_dense_T[cp.array(chunk_groups), cp.arange(actual_chunk)] = 1.0

        cp.cuda.Stream.null.synchronize()
        t_transfer_total += time.time() - t_xfer

        # -- Compute partial results on GPU --
        t_comp = time.time()

        # 1. Dot products: (groups x chunk) @ (chunk x genes) -> (groups x genes)
        # dense @ sparse is supported by CuPy and avoids spgemm OOM
        D_chunk = Lam_dense_T @ X_chunk_gpu  # (groups x genes)
        if cp_sparse.issparse(D_chunk):
            D_chunk = D_chunk.toarray()
        dot_accum += D_chunk.T  # -> (genes x groups)

        # 2. Squared gene norms: sum of X^2 per gene across chunk cells
        if is_sparse_input:
            Xsq = X_chunk_gpu.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum += cp.array(Xsq.sum(axis=0)).ravel()  # sum over cells (axis=0)
            del Xsq
        else:
            gene_sq_accum += cp.sum(X_chunk_gpu ** 2, axis=0)  # (genes,)

        # 3. Nonzero counts for expressed_pct filtering
        if nnz_accum is not None:
            if is_sparse_input:
                X_bin = X_chunk_gpu.copy()
                X_bin.data = cp.ones_like(X_bin.data)
                nnz_chunk = Lam_dense_T @ X_bin  # (groups x genes)
                if cp_sparse.issparse(nnz_chunk):
                    nnz_chunk = nnz_chunk.toarray()
                nnz_accum += nnz_chunk.T  # -> (genes x groups)
                del X_bin, nnz_chunk
            else:
                X_bin = (X_chunk_gpu > 0).astype(cp.float32)
                nnz_accum += (Lam_dense_T @ X_bin).T
                del X_bin

        cp.cuda.Stream.null.synchronize()
        t_compute_total += time.time() - t_comp

        # -- Free chunk memory --
        del X_chunk_gpu, D_chunk, Lam_dense_T
        cp.get_default_memory_pool().free_all_blocks()

        if verbosity > 1:
            print(f" ({time.time() - t_xfer:.3f}s)")

    if verbosity > 0:
        print(f"[gpu_cosg_chunked] Transfer total: {t_transfer_total:.3f}s")
        print(f"[gpu_cosg_chunked] Compute total:  {t_compute_total:.3f}s")

    # -- Compute cosine similarity from accumulated sums --
    t_post = time.time()

    gene_norms = cp.sqrt(gene_sq_accum)
    gene_norms = cp.maximum(gene_norms, 1e-12)
    group_norms_safe = cp.maximum(group_norms_gpu, 1e-12)

    cosine_sim = dot_accum / (gene_norms[:, None] * group_norms_safe[None, :])

    del dot_accum, gene_sq_accum

    # -- Apply COSG penalty formula --
    genexlambda = _apply_cosg_penalty_gpu(cosine_sim, mu)

    del cosine_sim

    # -- Apply expressed_pct filter --
    if nnz_accum is not None:
        for k in range(n_groups):
            n_cells_k = group_sizes_gpu[k]
            if n_cells_k > 0:
                threshold = max(
                    float(n_cells_k * expressed_pct),
                    float(expressed_min_num_cells_in_target_group)
                )
                genexlambda[nnz_accum[:, k] < threshold, k] = -cp.inf
        del nnz_accum

    if verbosity > 0:
        print(f"[gpu_cosg_chunked] Post-processing: {time.time() - t_post:.3f}s")

    # -- Transfer result to CPU --
    result = cp.asnumpy(genexlambda)

    del genexlambda
    cp.get_default_memory_pool().free_all_blocks()

    if verbosity > 0:
        print(f"[gpu_cosg_chunked] TOTAL: {time.time() - t0:.3f}s")

    return result


def gpu_core_cosg_full(
    cellxgene,
    cluster_mat,
    mu,
    n_genes_user,
    var_names,
    remove_lowly_expressed=False,
    expressed_pct=0.1,
    expressed_min_num_cells_in_target_group=3,
    is_sparse_input=True,
    verbosity=0,
):
    """
    Full GPU COSG pipeline: cosine sim + penalty + filtering + top-k.

    Keeps everything on GPU including expressed_pct filtering and top-k
    extraction, bypassing the CPU per-group loop entirely.

    Parameters
    ----------
    cellxgene : scipy.sparse.csr_matrix or np.ndarray, (n_cells, n_genes)
    cluster_mat : scipy.sparse.csr_matrix, (n_groups, n_cells)
    mu : float
    n_genes_user : int
    var_names : np.ndarray of str, shape (n_genes,)
    remove_lowly_expressed : bool
    expressed_pct : float
    expressed_min_num_cells_in_target_group : int
    is_sparse_input : bool
    verbosity : int

    Returns
    -------
    top_gene_names : np.ndarray, shape (n_genes_user, n_groups), dtype object
    top_gene_scores : np.ndarray, shape (n_genes_user, n_groups), dtype float32
    top_gene_indices : np.ndarray, shape (n_genes_user, n_groups), dtype int
    """
    t0 = time.time()

    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]

    if verbosity > 0:
        print(f"[gpu_cosg] Transferring {n_cells:,} × {n_genes:,} matrix to GPU...")

    # ── Transfer to GPU ──────────────────────────────────────────────────
    t_xfer = time.time()

    if is_sparse_input:
        X_gpu = _scipy_to_cupy_csr(cellxgene)
    else:
        X_gpu = cp.array(np.asarray(cellxgene, dtype=np.float32))

    cm_gpu = _scipy_to_cupy_csr(cluster_mat)

    if verbosity > 0:
        print(f"[gpu_cosg] Transfer to GPU:      {time.time() - t_xfer:.3f}s")

    # ── Cosine similarity: (genes × groups) ──────────────────────────────
    t_cos = time.time()

    cm_dense = cm_gpu.toarray()  # (n_groups, n_cells)
    if cp_sparse.issparse(X_gpu):
        dots = cm_dense @ X_gpu
    else:
        dots = cm_dense @ X_gpu
    dots = dots.T  # → (n_genes, n_groups)

    # L2 norms for X.T rows (= column norms of X)
    if cp_sparse.issparse(X_gpu):
        Xsq = X_gpu.copy()
        Xsq.data = Xsq.data ** 2
        gene_norms = cp.sqrt(cp.array(Xsq.sum(axis=0)).ravel())
        del Xsq
    else:
        gene_norms = cp.sqrt(cp.sum(X_gpu ** 2, axis=0))

    # L2 norms for cluster_mat rows
    cm_sq = cm_gpu.copy()
    cm_sq.data = cm_sq.data ** 2
    group_norms = cp.sqrt(cp.array(cm_sq.sum(axis=1)).ravel())
    del cm_sq

    gene_norms = cp.maximum(gene_norms, 1e-12)
    group_norms = cp.maximum(group_norms, 1e-12)

    dots /= gene_norms[:, cp.newaxis]
    dots /= group_norms[cp.newaxis, :]

    cosine_sim = dots
    del dots

    if verbosity > 0:
        print(f"[gpu_cosg] Cosine similarity:    {time.time() - t_cos:.3f}s")

    # ── COSG penalty formula ─────────────────────────────────────────────
    t_pen = time.time()

    genexlambda_gpu = _apply_cosg_penalty_gpu(cosine_sim, mu)
    del cosine_sim

    if verbosity > 0:
        print(f"[gpu_cosg] Penalty formula:      {time.time() - t_pen:.3f}s")

    # ── Expressed-pct filtering (GPU) ────────────────────────────────────
    t_filt = time.time()

    if remove_lowly_expressed and is_sparse_input:
        # Build binary matrix: 1 where nonzero
        X_binary_gpu = X_gpu.copy()
        X_binary_gpu.data = cp.ones_like(X_binary_gpu.data)

        # nnz per gene per group: cm_dense @ X_binary → (n_groups, n_genes)
        nnz_per_gene_per_group = cm_dense @ X_binary_gpu  # dense @ sparse → dense
        nnz_per_gene_per_group = nnz_per_gene_per_group.T  # → (n_genes, n_groups)
        del X_binary_gpu

        # n_cells per group
        n_cells_per_group = cp.array(cm_gpu.sum(axis=1)).ravel()  # (n_groups,)

        # Threshold: max(n_cells_per_group * expressed_pct, min_cells)
        thresholds = cp.maximum(
            n_cells_per_group * expressed_pct,
            expressed_min_num_cells_in_target_group
        )

        # Mask lowly expressed
        lowly_expressed_mask = nnz_per_gene_per_group < thresholds[cp.newaxis, :]
        genexlambda_gpu[lowly_expressed_mask] = -1.0
        del nnz_per_gene_per_group, lowly_expressed_mask
    elif remove_lowly_expressed and not is_sparse_input:
        # Dense path: count nonzeros per column per group
        X_binary = (X_gpu != 0).astype(cp.float32)
        nnz_per_gene_per_group = (cm_dense @ X_binary).T  # (n_genes, n_groups)
        del X_binary

        n_cells_per_group = cp.array(cm_gpu.sum(axis=1)).ravel()
        thresholds = cp.maximum(
            n_cells_per_group * expressed_pct,
            expressed_min_num_cells_in_target_group
        )
        lowly_expressed_mask = nnz_per_gene_per_group < thresholds[cp.newaxis, :]
        genexlambda_gpu[lowly_expressed_mask] = -1.0
        del nnz_per_gene_per_group, lowly_expressed_mask

    del X_gpu, cm_gpu, cm_dense

    if verbosity > 0:
        print(f"[gpu_cosg] Expressed filtering:  {time.time() - t_filt:.3f}s")

    # ── Top-k extraction (GPU) ───────────────────────────────────────────
    t_topk = time.time()

    # Clamp n_genes_user to available genes
    actual_k = min(n_genes_user, n_genes)

    # argpartition for top-k per column
    top_k_idx = cp.argpartition(genexlambda_gpu, -actual_k, axis=0)[-actual_k:]

    # Gather scores at top-k positions
    top_k_scores = cp.take_along_axis(genexlambda_gpu, top_k_idx, axis=0)

    # Sort within top-k (descending)
    sort_order = cp.argsort(-top_k_scores, axis=0)
    top_k_idx_sorted = cp.take_along_axis(top_k_idx, sort_order, axis=0)
    top_k_scores_sorted = cp.take_along_axis(top_k_scores, sort_order, axis=0)

    del genexlambda_gpu, top_k_idx, top_k_scores, sort_order

    if verbosity > 0:
        print(f"[gpu_cosg] Top-k extraction:     {time.time() - t_topk:.3f}s")

    # ── Transfer only small top-k results to CPU ─────────────────────────
    t_back = time.time()

    top_k_idx_cpu = cp.asnumpy(top_k_idx_sorted)        # (n_genes_user, n_groups)
    top_k_scores_cpu = cp.asnumpy(top_k_scores_sorted)   # (n_genes_user, n_groups)

    del top_k_idx_sorted, top_k_scores_sorted
    cp.get_default_memory_pool().free_all_blocks()

    # Map indices to gene names
    top_gene_names = var_names[top_k_idx_cpu]  # (n_genes_user, n_groups)

    if verbosity > 0:
        print(f"[gpu_cosg] Transfer to CPU:      {time.time() - t_back:.3f}s")
        print(f"[gpu_cosg] Total GPU pipeline:   {time.time() - t0:.3f}s")

    return top_gene_names, top_k_scores_cpu.astype(np.float32), top_k_idx_cpu

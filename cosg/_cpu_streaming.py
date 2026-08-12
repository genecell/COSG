"""
CPU streaming core computation for COSG.

Chunked cell-axis accumulation using numpy/scipy. Processes cells in blocks,
accumulating partial results to avoid holding the full matrix transpose or
dense cosine similarity matrix in memory.

Peak memory: O(chunk_size * n_genes) + O(n_genes * n_groups) accumulators
instead of O(n_cells * n_genes) for the full matrix.
"""

import numpy as np
import time
from scipy import sparse as sp_sparse


def cpu_core_cosg_chunked(cellxgene, cluster_mat, mu, is_sparse_input,
                          remove_lowly_expressed=True, expressed_pct=0.1,
                          expressed_min_num_cells_in_target_group=3,
                          chunk_size=None, verbosity=0):
    """
    CPU COSG with chunked cell-axis streaming for low-memory environments.

    Drop-in replacement for the CPU cosine_similarity + penalty block in cosg.py.
    Same interface as gpu_core_cosg_chunked() but uses numpy/scipy.

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
        Absolute minimum nonzero count threshold.
    chunk_size : int or None
        Number of cells per chunk. None = auto-select based on available RAM.
    verbosity : int
        0 = silent, 1 = progress, 2 = per-chunk timing.

    Returns
    -------
    genexlambda : np.ndarray, shape (n_genes, n_groups)
        COSG penalty scores. Dense numpy array.
    """
    t0 = time.time()

    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]

    # ── Auto-select chunk size ──────────────────────────────────────────
    if chunk_size is None:
        chunk_size = _auto_chunk_size_cpu(n_cells, n_genes, is_sparse_input,
                                          cellxgene)
    chunk_size = min(chunk_size, n_cells)
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    if verbosity > 0:
        print(f"[cpu_cosg_chunked] {n_cells:,} cells x {n_genes:,} genes, "
              f"{n_groups} groups")
        print(f"[cpu_cosg_chunked] Chunk size: {chunk_size:,} cells, "
              f"{n_chunks} chunks")

    # ── Prepare cell-to-group mapping ───────────────────────────────────
    # cluster_mat is (n_groups x n_cells) with exactly one 1 per column
    cm_csc = cluster_mat.tocsc()
    cell_group_indices = np.array(cm_csc.indices, dtype=np.int32)

    # Group sizes
    group_sizes = np.array(cluster_mat.sum(axis=1)).ravel().astype(np.float64)
    group_norms = np.sqrt(group_sizes)

    # ── Initialize accumulators ─────────────────────────────────────────
    # Use float64 for accumulation precision (matches original CPU path)
    dot_accum = np.zeros((n_genes, n_groups), dtype=np.float64)
    gene_sq_accum = np.zeros(n_genes, dtype=np.float64)

    if remove_lowly_expressed and expressed_pct > 0:
        nnz_accum = np.zeros((n_genes, n_groups), dtype=np.float64)
    else:
        nnz_accum = None

    # ── Ensure CSR format ───────────────────────────────────────────────
    if is_sparse_input:
        cellxgene_csr = cellxgene.tocsr()
        if cellxgene_csr.dtype != np.float32 and cellxgene_csr.dtype != np.float64:
            cellxgene_csr = cellxgene_csr.astype(np.float64)
    else:
        cellxgene_np = np.asarray(cellxgene, dtype=np.float64)

    t_compute = 0.0

    # ── Stream chunks ───────────────────────────────────────────────────
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        actual_chunk = end - start

        if verbosity > 1:
            t_chunk = time.time()
            print(f"  Chunk {chunk_idx+1}/{n_chunks}: cells [{start:,}:{end:,}]",
                  end="", flush=True)

        t_comp = time.time()

        # ── Slice chunk (CSR row slice is O(nnz_in_range)) ──────────────
        if is_sparse_input:
            X_chunk = cellxgene_csr[start:end]  # (chunk x genes), CSR
        else:
            X_chunk = cellxgene_np[start:end]    # (chunk x genes), dense

        # ── Build indicator matrix for this chunk ───────────────────────
        chunk_groups = cell_group_indices[start:end]
        Lam_chunk = sp_sparse.csr_matrix(
            (np.ones(actual_chunk, dtype=np.float64),
             (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
            shape=(actual_chunk, n_groups)
        )

        # ── Accumulate dot products ─────────────────────────────────────
        # (genes x chunk) @ (chunk x groups) -> (genes x groups)
        if is_sparse_input:
            D_chunk = X_chunk.T @ Lam_chunk  # scipy sparse @ sparse
            if sp_sparse.issparse(D_chunk):
                D_chunk = D_chunk.toarray()
        else:
            D_chunk = X_chunk.T @ Lam_chunk.toarray()
        dot_accum += D_chunk

        # ── Accumulate squared gene norms ───────────────────────────────
        if is_sparse_input:
            Xsq = X_chunk.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum += np.asarray(Xsq.sum(axis=0)).ravel()
            del Xsq
        else:
            gene_sq_accum += np.sum(X_chunk ** 2, axis=0)

        # ── Accumulate nonzero counts ───────────────────────────────────
        if nnz_accum is not None:
            if is_sparse_input:
                X_bin = X_chunk.copy()
                X_bin.data = np.ones_like(X_bin.data)
                nnz_chunk = X_bin.T @ Lam_chunk
                if sp_sparse.issparse(nnz_chunk):
                    nnz_chunk = nnz_chunk.toarray()
                nnz_accum += nnz_chunk
                del X_bin, nnz_chunk
            else:
                X_bin = (X_chunk > 0).astype(np.float64)
                nnz_accum += X_bin.T @ Lam_chunk.toarray()
                del X_bin

        t_compute += time.time() - t_comp

        del X_chunk, Lam_chunk, D_chunk

        if verbosity > 1:
            print(f" ({time.time() - t_chunk:.3f}s)")

    if verbosity > 0:
        print(f"[cpu_cosg_chunked] Chunk compute total: {t_compute:.3f}s")

    # ── Compute cosine similarity from accumulated sums ─────────────────
    t_post = time.time()

    gene_norms = np.sqrt(gene_sq_accum)
    gene_norms = np.maximum(gene_norms, 1e-12)
    group_norms_safe = np.maximum(group_norms, 1e-12)

    cosine_sim = dot_accum / (gene_norms[:, None] * group_norms_safe[None, :])

    del dot_accum, gene_sq_accum

    # ── Apply COSG penalty formula ──────────────────────────────────────
    genexlambda = _apply_cosg_penalty_cpu(cosine_sim, mu)

    del cosine_sim

    # ── Apply expressed_pct filter ──────────────────────────────────────
    if nnz_accum is not None:
        for k in range(n_groups):
            n_cells_k = group_sizes[k]
            if n_cells_k > 0:
                threshold = max(
                    n_cells_k * expressed_pct,
                    expressed_min_num_cells_in_target_group
                )
                genexlambda[nnz_accum[:, k] < threshold, k] = -1
        del nnz_accum

    if verbosity > 0:
        print(f"[cpu_cosg_chunked] Post-processing: {time.time() - t_post:.3f}s")
        print(f"[cpu_cosg_chunked] TOTAL: {time.time() - t0:.3f}s")

    return genexlambda


def cpu_cosine_sim_chunked(cellxgene, cluster_mat, is_sparse_input,
                           chunk_size=None, verbosity=0):
    """
    Compute cosine similarity via chunked cell-axis streaming (CPU).

    Returns the cosine_sim matrix without applying the penalty formula.
    Used for batch mode where cosine similarities are averaged across batches
    before applying the penalty.

    Parameters
    ----------
    cellxgene : scipy.sparse or np.ndarray, shape (n_cells, n_genes)
    cluster_mat : scipy.sparse, shape (n_groups, n_cells)
    is_sparse_input : bool
    chunk_size : int or None
    verbosity : int

    Returns
    -------
    cosine_sim : np.ndarray, shape (n_genes, n_groups), float64
    """
    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]

    if chunk_size is None:
        chunk_size = _auto_chunk_size_cpu(n_cells, n_genes, is_sparse_input,
                                          cellxgene)
    chunk_size = min(chunk_size, n_cells)
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    cm_csc = cluster_mat.tocsc()
    cell_group_indices = np.array(cm_csc.indices, dtype=np.int32)

    group_sizes = np.array(cluster_mat.sum(axis=1)).ravel().astype(np.float64)
    group_norms = np.sqrt(group_sizes)

    dot_accum = np.zeros((n_genes, n_groups), dtype=np.float64)
    gene_sq_accum = np.zeros(n_genes, dtype=np.float64)

    if is_sparse_input:
        cellxgene_csr = cellxgene.tocsr()
        if cellxgene_csr.dtype != np.float32 and cellxgene_csr.dtype != np.float64:
            cellxgene_csr = cellxgene_csr.astype(np.float64)
    else:
        cellxgene_np = np.asarray(cellxgene, dtype=np.float64)

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        actual_chunk = end - start

        if is_sparse_input:
            X_chunk = cellxgene_csr[start:end]
        else:
            X_chunk = cellxgene_np[start:end]

        chunk_groups = cell_group_indices[start:end]
        Lam_chunk = sp_sparse.csr_matrix(
            (np.ones(actual_chunk, dtype=np.float64),
             (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
            shape=(actual_chunk, n_groups)
        )

        if is_sparse_input:
            D_chunk = X_chunk.T @ Lam_chunk
            if sp_sparse.issparse(D_chunk):
                D_chunk = D_chunk.toarray()
        else:
            D_chunk = X_chunk.T @ Lam_chunk.toarray()
        dot_accum += D_chunk

        if is_sparse_input:
            Xsq = X_chunk.copy()
            Xsq.data = Xsq.data ** 2
            gene_sq_accum += np.asarray(Xsq.sum(axis=0)).ravel()
            del Xsq
        else:
            gene_sq_accum += np.sum(X_chunk ** 2, axis=0)

        del X_chunk, Lam_chunk, D_chunk

    gene_norms = np.sqrt(gene_sq_accum)
    gene_norms = np.maximum(gene_norms, 1e-12)
    group_norms_safe = np.maximum(group_norms, 1e-12)

    cosine_sim = dot_accum / (gene_norms[:, None] * group_norms_safe[None, :])
    return cosine_sim


def cpu_nnz_counts_chunked(cellxgene, cluster_mat, is_sparse_input,
                           chunk_size=None, verbosity=0):
    """
    Compute per-gene per-group nonzero counts via chunked streaming.

    Returns nnz_counts[g, k] = number of cells in group k where gene g is nonzero.
    Used for expressed_pct filtering in batch mode where nnz is not accumulated
    during cosine similarity computation.

    Parameters
    ----------
    cellxgene : sparse or dense, shape (n_cells, n_genes)
    cluster_mat : sparse, shape (n_groups, n_cells)
    is_sparse_input : bool
    chunk_size : int or None
    verbosity : int

    Returns
    -------
    nnz_counts : np.ndarray, shape (n_genes, n_groups), float64
    """
    n_cells, n_genes = cellxgene.shape
    n_groups = cluster_mat.shape[0]

    if chunk_size is None:
        chunk_size = _auto_chunk_size_cpu(n_cells, n_genes, is_sparse_input,
                                          cellxgene)
    chunk_size = min(chunk_size, n_cells)
    n_chunks = (n_cells + chunk_size - 1) // chunk_size

    cm_csc = cluster_mat.tocsc()
    cell_group_indices = np.array(cm_csc.indices, dtype=np.int32)

    nnz_accum = np.zeros((n_genes, n_groups), dtype=np.float64)

    if is_sparse_input:
        cellxgene_csr = cellxgene.tocsr()

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_cells)
        actual_chunk = end - start

        if is_sparse_input:
            X_chunk = cellxgene_csr[start:end]
            X_bin = X_chunk.copy()
            X_bin.data = np.ones_like(X_bin.data)
        else:
            X_bin = (cellxgene[start:end] > 0).astype(np.float64)

        chunk_groups = cell_group_indices[start:end]
        Lam_chunk = sp_sparse.csr_matrix(
            (np.ones(actual_chunk, dtype=np.float64),
             (np.arange(actual_chunk, dtype=np.int32), chunk_groups)),
            shape=(actual_chunk, n_groups)
        )

        nnz_chunk = X_bin.T @ Lam_chunk
        if sp_sparse.issparse(nnz_chunk):
            nnz_chunk = nnz_chunk.toarray()
        nnz_accum += nnz_chunk

        del X_bin, Lam_chunk, nnz_chunk

    return nnz_accum


def _apply_cosg_penalty_cpu(cosine_sim, mu):
    """
    Apply the COSG penalty formula on CPU.

    Replicates the exact logic from cosg.py lines ~381-445.

    For sparse input:
        e^2 = cosine_sim^2
        sum_e^2 = row-wise sum of e^2
        if mu == 1: score = e^2 / sum_e^2 * cosine_sim = cosine_sim^3 / sum_e^2
        else: score = e^2 / ((1-mu)*e^2 + mu*sum_e^2) * cosine_sim

    For dense input: same formula applied only at nonzero positions.
    """
    e_power2 = np.multiply(cosine_sim, cosine_sim)
    e_power2_sum = np.sum(e_power2, axis=1)  # (n_genes,)

    # Mask of nonzero positions
    pos_nonzero = cosine_sim != 0

    if mu == 1:
        safe_sum = np.maximum(e_power2_sum, 1e-30)
        genexlambda = e_power2 / safe_sum[:, None]
        genexlambda = np.where(pos_nonzero, np.multiply(genexlambda, cosine_sim), 0.0)
    else:
        denom = (1 - mu) * e_power2 + mu * e_power2_sum[:, None]
        safe_denom = np.where(pos_nonzero, denom, 1.0)
        genexlambda = np.where(pos_nonzero, e_power2 / safe_denom, 0.0)
        genexlambda = np.where(pos_nonzero, np.multiply(genexlambda, cosine_sim), 0.0)

    return genexlambda


def _auto_chunk_size_cpu(n_cells, n_genes, is_sparse, cellxgene=None,
                         target_ram_mb=2000, min_chunk=1000, max_chunk=100_000):
    """
    Auto-select chunk size for CPU streaming.

    Strategy: Keep per-chunk RAM usage under target_ram_mb.

    Per-chunk memory:
      - Sparse: nnz_per_cell * chunk_size * 12 bytes (data + indices + overhead)
                + transpose + indicator + intermediates
      - Dense: chunk_size * n_genes * 8 bytes * ~3 (chunk + transpose + intermediates)

    Parameters
    ----------
    n_cells : int
    n_genes : int
    is_sparse : bool
    cellxgene : sparse matrix or None
    target_ram_mb : int
        Target per-chunk RAM in MB (default 2000 = 2 GB).
    min_chunk : int
    max_chunk : int

    Returns
    -------
    chunk_size : int
    """
    target_bytes = target_ram_mb * 1024 * 1024

    if is_sparse and cellxgene is not None:
        avg_nnz_per_cell = cellxgene.nnz / cellxgene.shape[0]
        # data (8 bytes float64) + indices (4 bytes) + overhead
        # x ~3 for chunk + transpose + binary copy
        bytes_per_cell = avg_nnz_per_cell * 12 * 3
    else:
        # Dense: chunk x genes x 8 bytes x 3 (chunk + transpose + intermediates)
        bytes_per_cell = n_genes * 8 * 3

    chunk_size = int(target_bytes / max(bytes_per_cell, 1))
    chunk_size = max(min_chunk, min(chunk_size, max_chunk))

    return chunk_size

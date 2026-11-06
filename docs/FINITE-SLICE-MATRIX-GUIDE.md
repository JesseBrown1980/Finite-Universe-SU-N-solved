# Finite universe, one spatial slice at a time: exact matrix laws

This guide covers the **2×2×2** periodic spatial graph with 24 *oriented* links and 8 vertices; its time evolution is a four-dimensional finite transfer problem. The gauge groups Z₃ and non-Abelian S₃ are finite subgroups of SU(3), **not** continuous SU(3). The number of finite time slices can be chosen without making the spatial state space infinite.

## Law A — remove vertex gauge redundancy exactly

Choose one spanning tree of 7 edges. For every oriented tree edge parent→child, recursively set the local gauge choice so its link becomes the identity. The remaining E−V+1=24−8+1=**17** cycle variables are complete; the root gauge freedom acts by **one common conjugation** of all 17.

Thus the physical orbit count is the number of simultaneous conjugacy classes of 17-tuples:

    dim(H_gauge) = (1/|G|) sum_{g in G} |C_G(g)|^17

The centralizer sizes are 3,3,3 for Z₃; for S₃ they are 6 once, 2 three times, and 3 twice. Hence:

    Z3: 3^17 = 129140163
    S3: (6^17 + 3*2^17 + 2*3^17)/6 = 2821153019713

`orbit_reduce.py` computes the formula from multiplication tables, checks it against full conjugation-orbit enumeration for 0–4 cycles, and tests deterministic full-cube gauge changes. For S₃, quotienting by `|G|^8` instead of accounting for the root stabilizer gives an incorrect count; the checker rejects it.

These are *exact dimensions*, not an exponentially small representation. Large dimensions remain.

## Law B — eliminate temporal links before forming matrices

For a given before/after spatial boundary, there are 8 temporal link variables. Because the spatial side length is 2, the 24 temporal plaquettes group into 12 **pair factors**. Exact variable elimination sums over the 8 links with maximum intermediate scope 4. No `|G|^8` tuple enumeration is necessary once the pair factors have been built.

Use the independent bipartite contraction and brute-force checks to prove the local contraction. Cache a factor by (group multiplication table hash, plaquette weight hash, ordered boundary holonomies, edge and orientation). Apply spatial gauge fixing **before** cache-key creation, while preserving orbit multiplicities. A code-only or hash-only match without verified semantic invariance is not a valid cache hit.

## Law C — jump only a proven closed evolution

For a fixed transfer map `T` on a **proven invariant state representation**, compute `T^n` by repeated squaring/cubing and verify against sequential updates at small n. The tested four-state catalogue is **not** invariant, so its powers are diagnostic: they do not give the full finite-universe spectral or partition results. An exact shift of the time index must also preserve the complete boundary state and its gauge constraints.

No theorem says that the exact tensor-network bond dimension stays small. Contracting a fixed boundary tensor can still grow rapidly in space or time. If approximating, record a certified, additive error bound for observables and eigenvalues. **Freeze and play** reuses independently authenticated full boundary states; the internal SHA chain alone does not prove origin.

## Law D — build the closed gauge-invariant transfer operator

Next implementation target: construct and test an operator on character/irrep-decomposed link variables with Gauss-law projectors at each vertex. For Z₃, links have three one-dimensional irreducible characters, and Gauss constraints are trit-flux conservation at each vertex. For S₃, the irreducible representations have dimensions 1, 1 and 2; coupling their indices at the vertices produces gauge-invariant spin-network intertwiners. Matrix-vector products can be represented as contractions of local transfer tensors, without materializing a `dim(H_gauge)` square matrix.

Acceptance gates:

- Independence: compare the character / intertwiner contraction with the original exact finite-group plaquette formula at small lattices and boundary samples.
- Closure: prove the transfer map preserves Gauss-law sectors and the chosen compressed representation; never use a diagnostic boundary catalogue as a closed sector.
- Exactness: compare full observable output and hash receipts across direct updates and matrix jumps on all tractable instances.
- Scale: record bond dimensions and peak intermediate factor size, not just wall time; stop if memory grows unbounded.
- Spectral proof: only call an eigenvalue a finite gauge mass if it belongs to the physical transfer sector and is certified with an independently checked residual or inertia bound.

## Work status

The complete transfer *entry* is implemented and independently verified. The gauge-orbit count is implemented and independently checked. **No full 2×2×2 transfer eigensolve or full gauge-state propagation has been certified.** This is the honest boundary of the current finite-universe matrix result. Independent Opus review remains required before declaring a new closed-sector theorem.

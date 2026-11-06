# Full gauge-space finite 4D slice transfer: exact tensor solution to state storage

**2026-11-06, prediction-first.** This package is a continuation of the merged
[`2×2×2` exact-slice test](../finite-4d-slices-2x2x2-20261106/README.md)
and the [finite-slice matrix guide](../../docs/FINITE-SLICE-MATRIX-GUIDE.md).
It addresses the exponentially large **representation and exact matrix-free application** problem, not the unsolved numerical interacting mass spectrum. The prediction was committed at `77454a7e2700216cdcc8efb548a175d79b17fa2c` before source. The **six independently computed source/result digests** were pinned in the separate Git commit `508ba40d8734414fba4de186e3efff64dbfea5c9` before the implementation commit.

## Measured complete-space results

| Quantity | Z₃ (center subgroup) | S₃ (non-Abelian subgroup) |
|---|---:|---:|
| Raw 24-link configurations | 282,429,536,481 | 4,738,381,338,321,616,896 |
| Physical Gauss-gauge orbits | 129,140,163 | 2,821,153,019,713 |
| Number of stored entries of full state | **0** | **0** |
| Exact (one step) transfer of uniform input | `M(U)·3⁸·6²⁴` | `M(U)·6⁸·12²⁴` |
| Independently checked full transfer entries | 3 plus gauge transforms | 3 plus gauge transforms |
| Full electric-only first nonvacuum eigenvalue ratio | **1/4** | **1/4** |
| Interacting transfer gap | **Strictly positive; tiny rigorous bound** | **Strictly positive; tiny rigorous bound** |

Here `M(U)` is the product of the **24** spatial Wilson plaquette weights.
Each transfer includes **24** temporal plaquettes, so **all six 4D plane orientations** have equal, positive integer weights, with the same finite group on every link. The two directions of a periodic two-site edge are distinct.

For Z₃ only, an exact duality reduces the full magnetic partition to **59,049** admissible closed 2-form flux assignments (`rank ∂₂ = 14`, `dim ker ∂₂ = 10`):

| Full-space result | Exact value |
|---|---:|
| Sum over dual closed flux weights | **99,186,762** |
| Sum over all 3²⁴ raw spatial link states of `M(U)` | **28,013,271,216,711,264,522** |
| Sum over the physical 3¹⁷ gauge configurations of `M(U)` | **12,808,994,612,122,206** |
| Complete one-open-time-slab partition, raw-link counting | **870,891,141,386,018,892,529,920,237,038,951,498,514,432** |
| Normalized spatial plaquette character (`Re Tr Uₚ/3`) | **27985349/33062254** |
| Two selected plaquette characters | **37034332/49593381** |

This is a **complete finite full-space observable**, not a small projected matrix. The open-slab partition has free independent spatial fields on the two boundary slices; it is not `Tr(T^Nₜ)` on a time-periodic lattice. The Z₃ plaquette expectation is the spatial marginal of this open-slab measure.

The 59,049 dual fluxes and Wilson insertions are enumerated in `gauge_tensor.py`. The independent `primal_z3_check.c` does **not use Fourier characters**: it enumerates all **4,782,969** distinct plaquette-flux vectors in the image of the spatial link-boundary map. It reproduces the full partition and both Wilson numerator sums **exactly** with C11 integer arithmetic.

## The exact one-slice transfer law

Write `G=Z₃` or `S₃`, `q=|G|`, and let `U_e, V_e` be the 24 spatial link variables before and after one slice. With eight temporal links `g_x`, define

```
M(U) = ∏_{24 spatial plaquettes p} w(hol_p(U))
K(U,V) = Σ_{g ∈ G⁸} ∏_{e=(x,y)} w(U_e g_y V_e⁻¹ g_x⁻¹)
T(U,V) = M(U) K(U,V).
```

Let `C` be the product of **24 individual q×q link convolution matrices** `C_e[a,b]=w(a b⁻¹)`, and `P` be the *full, exact eight-vertex gauge average*. As an operator,

```
T = M · q⁸ P C = M · q⁸ P C P
```

when acting on gauge-invariant inputs (the second equality is for the gauge-invariant sector). Centrality of `w` makes `C` commute with `P`, while gauge invariance makes `M` commute with `P`. **The full gauge-invariant space is therefore exactly closed**, without a physical-state table, a truncated basis, or a bounded tensor rank assumption.

For a factorized input `f(V)=∏_e f_e(V_e)`, the 24 sums over `V_e` become local link factors. The eight temporal variables are eliminated by a bipartite cube contraction. This computes any selected output component `(Tf)(U)` with **no `q²⁴` intermediate**. The independent verifier compares exact output components to the previously merged four-dimensional plaquette kernel and tests independent gauge transforms on each boundary.

For the uniform input `f(V)=1` the algebra closes completely:

```
(T1)(U) = M(U) · q⁸ · (Σ_g w(g))²⁴.
```

This is an exact compressed wavefunction over **every one of the trillions of physical S₃ states**. Applying another `T` gives an exact finite tensor circuit but may increase its contraction width. We **do not** claim a constant-size exact representation after arbitrary time evolution.

## Exact spectrum facts, without pretending the full interacting spectrum is computed

The local link-convolution eigenvalues are:

- `Z3`: 6, 3, 3;
- `S3`: 12 (trivial), 6 (sign), 3 (two-dimensional standard irrep, multiplicity four).

All are positive, and a nontrivial gauge-invariant network on this graph uses **at least two edges**. The length-two cycle with nontrivial Z₃ flux or the S₃ sign representation realizes the leading nonvacuum eigenvalue ratio `1/4`. Thus the **exact full physical electric-only** transfer gap is `ln 4` for *both* finite groups. The independent checker exhaustively reproduces the needed eigenmode on the two-link gauge-invariant graph (`9` or `36` configurations).

For the **interacting isotropic** model, every temporal kernel entry is between `q⁸` and `q⁸ w_max²⁴`. Perron normalization followed by a common Doeblin minorization yields an explicit *global* transfer ratio bound,

```
0 < lambda_1/lambda_0 <= (w_max²⁴ - 1) / w_max²⁴,
```

where `lambda_1` is the largest nonvacuum eigenvalue on the **complete physical gauge-invariant space**. Therefore the strictly positive finite transfer gap is at least `log(w_max²⁴/(w_max²⁴-1))`, which is extremely small:

- `Z3`: ratio ≤ `281474976710655/281474976710656`;
- `S3`: ratio ≤ `59604644775390624/59604644775390625`.

The similarity with the positive symmetric operator `M^(1/2) K M^(1/2)` and positivity of the local convolution spectra establish that all eigenvalues are positive. This bound is **not** a computed physical glueball mass or an estimate of the actual gap. It is a rigorous finite positive-gap certificate for the specified discrete model, independent of the number of physical states.

## Verification and provenance

From the repository root:

```sh
p=results/finite-4d-gauge-tensor-20261106
python3 "$p/gauge_tensor.py" --output /tmp/gauge-tensor-result.hbp
cmp "$p/RESULT.hbp" /tmp/gauge-tensor-result.hbp
python3 "$p/verify_gauge_tensor.py" > /tmp/gauge-tensor-review.hbp
cmp "$p/REVIEW.hbp" /tmp/gauge-tensor-review.hbp
cc -std=c11 -O2 -Wall -Wextra -Werror "$p/primal_z3_check.c" -o /tmp/z3-primal
/tmp/z3-primal > /tmp/z3-primal.hbp
cmp "$p/PRIMAL.hbp" /tmp/z3-primal.hbp
(cd "$p" && sha256sum --strict -c SHA256SUMS)
```

The CI workflow checks the **earlier Git-pinned anchor commit** `508ba40d8734414fba4de186e3efff64dbfea5c9` separately from the current branch files. The anchor is a content-addressed Git record, **not a cryptographic signature of who authored the files**. The original predictions and old frozen receipts are left unchanged.

## Unsolved problems

The interacting **numeric** 2×2×2 `Z3` and `S3` mass eigenvalues, static-charge potentials and multi-slice time-periodic partition functions remain OPEN. For `S3`, a full interacting partition value is not yet evaluated. For arbitrary evolved input, exact factor networks can acquire very high width: no assertion is made that a polynomial-time contraction exists. Opus independent review is required before this result is considered accepted. Full continuous SU(3) and physical QCD are **not** claimed.

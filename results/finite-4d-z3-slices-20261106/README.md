# Finite four-dimensional universe: exact one-slice transfer (2026-11-06)

**Status: implementation and local exact checks PASS; independent Opus review OPEN.**

This is an intentionally small, fully discrete **Euclidean** 4D lattice gauge model with link group **Z3**, embedded as the center of SU(3). It is **not** full continuous SU(3), not the 27-element Heisenberg subgroup H3-27, not dynamical QCD, and not a physical 4D glueball-mass determination. The spatial periodic lattice is **2 × 1 × 1** vertices and the temporal circle is **Nt = 6** slices for its principal partition calculation (with Nt = 2, 3 and 9 controls). The length-one spatial directions are degenerate. No infinite-volume or continuum limit is required or claimed: the requested object is a finite matrix system.

## Frozen prediction first

The independent [PREDICTION.hbp](PREDICTION.hbp) was committed at `451e0249c74fda758877f1b52cbf8da83e072428` **before** the solver. It specifies exact matrix identities, the Polyakov static-charge ratios, the scalar symmetry-sector polynomial, negative controls, and excluded physical claims. All receipt lines have `json=0`. The predicted scalar quantity is the **transfer-channel eigenvalue gap**, not the mass of full SU(3) glueballs.

## The four-dimensional discrete Wilson action

Let each oriented link carry `a ∈ {0,1,2}` (the SU(3)-center matrix `ω^a I₃` with `ω³=1`). For **all six** plane orientations `(xy,xz,xt,yz,yt,zt)` use the same nonnegative plaquette weight `W(a_p)=4` for `a_p=0`, else `1`; this equals the standard SU(3) Wilson *relative weight* restricted to the Z3 center at `β=(4/3)ln(2)`. Normalization common to every plaquette is immaterial to the potential and normalized correlator. The action is **not** a unitary real-time evolution.

Each spatial slice has two vertices and six link variables. Its full discrete configuration space has `3⁶=729` elements; local Z3 gauge transformations identify triples, leaving **243** orbit states. Fix `x₀=0`, retain `h=x₀+x₁`, `y₀,y₁,z₀,z₁ ∈ Z3`. The temporal gauge sum uses `d=g₁−g₀` and all six time-like plaquettes; every slice also includes the six spatial plaquettes. When passing to the gauge-orbit transfer, the omitted constant gauge-volume contribution gives `Z_full=9^Nt Tr(T^Nt)` for the unnormalized full-link partition. Such constants cancel normalized ratios.

The exact symmetric transfer entry is

`T(C,C') = 2^(M(C)+M(C')) * Σ_(d=0,1,2) 4^(number of zero-flux temporal plaquettes(C,C',d))`,

where `M` counts the zero-flux *spatial* plaquettes. Its 243×243 array is not needed for the fast computation because it **exactly factors**:

`T = 16 * K_x ⊗ B_y ⊗ B_z`, where `K_x` is 3×3 and each `B` is 9×9.

The engine applies `T` to a 243-entry integer state **one time slice at a time**, writing an SHA-256-linked row for each state from `t=0` through `t=81`. Five power-jump checks reproduce the individually propagated state at `t=3,6,9,27,81`; a stored checkpoint at `t=6` also reproduces `t=81` via a 75-slice jump. These are *exact* integer equalities, not rounded scientific observables.

## Results from the finite transfer matrix

The *external fundamental charge-pair* Polyakov-line correlator has exactly the ratio

`<P(0) P(1)*>_Nt = (2*2^Nt+1)/(4^Nt+2)`.

| Finite Nt | Exact Polyakov pair ratio | Effective potential `-ln(ratio)/Nt` (illustrative decimal) |
|---:|---:|---:|
| 2 | 1/2 | 0.346574 |
| 3 | 17/66 | 0.452147 |
| 6 | 43/1366 | 0.576407 |
| 9 | 1025/262146 | 0.616023 |

The finite transfer eigenvalue ratio in the external-charge channel is exactly `18/36 = 1/2`, so the *finite-model spectral static-pair potential* is `ln(2)` in lattice-time units. This statement concerns the fixed finite spatial lattice and a comparison of charge sectors, not an infinite spatial system.

The invariant local-spatial-plaquette scalar sector is controlled by the exact quadratic `λ² − 315λ + 5184`. Its two relevant eigenvalues are `(315 ± √78489)/2`. The corresponding transfer-channel scalar gap is exactly `ln[(315+√78489)/(315−√78489)]`, approximately **2.838030** in lattice units. It is a **Z3 finite-slice scalar gap**, *not* a 4D SU(3) glueball mass. Six finite-`Nt` connected spatial-plaquette correlators are reported as exact reduced rational numbers in `RESULT.hbp` and checked independently; short-time effective-mass ratios need not equal this eigenvalue gap because finite temporal wraparound and other finite sectors contribute.

## Independent checks and integrity

The implementation's **59,049** 243×243 entries were compared against raw 4D plaquette construction. Direct dense propagation and factored propagation agree for six steps; the independent `verify_slices.py` constructs the *unfactored* transfer separately and replays **all 81 transitions**, matching all archived state digests and the saved t=6 checkpoint. It independently recomputes six scalar correlators. It refuses altered leaves and altered internally valid chains when checked against the Git-pinned expected root, and rejects changed orientation, plaquette, and charge definitions. `SHA256SUMS`, `SU3-SLICES.hbi`, `ROOT.hbp` and the separate commit SHA act as the evidence and external anchor. Cryptographic hashing alone does not prove that the underlying scientific model is correct.

**Measured implementation speed (same container, three trials):** 81 complete dense-matrix slice updates took a median **523,938,769 ns**, versus **1,916,661 ns** for the factored 81-step matrix-power jump, a measured **273×** speedup *excluding* the **246,108,636 ns** dense-matrix construction. The values were exactly equal. This is an algorithm benchmark for this 243-state matrix, not a speedup claim for full continuous SU(3), a larger finite gauge group, or nondegenerate spatial lattices. See `BENCHMARK.hbp` for individual timings.

The original prediction, the preliminary floating-point 4D Monte Carlo experiment, the continuous SU(3) leading-order benchmark, and these exact Z3 finite-slice results are **separate experiments**. This implementation does not modify any frozen blind run.

## Reproduce

From this directory with Python 3.10 or later, without third-party dependencies:

```sh
python3 -I slice_engine.py
python3 -I verify_slices.py --root 72a1eec4275eb719ce70142c09270c266a765b96726e4c309b04218442e5ab8f
sha256sum --check SHA256SUMS
```

The root in the command is fixed independently of the leaves by the reviewed Git commit containing this README/workflow; a future changed run requires a new pinned root and review. Negative controls are recorded in `NEGATIVE-CONTROLS.hbp` and exercised by both implementations.

**Next:** independently review the gauge-volume normalization and charged Polyakov insertion, then scale to a nondegenerate spatial geometry and a non-Abelian finite subgroup of SU(3) under a separately frozen protocol. That next stage may need tensor/network boundary compression; factorization here is unusually favorable.

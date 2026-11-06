# Exact finite-universe 4D slices on 2×2×2 space

**Status (2026-11-06):** exact entry kernel, four-state *diagnostic projection*, and exact gauge-orbit counts **PASS locally**. The full gauge-invariant transfer spectrum is **OPEN**. Prediction frozen in commit `5c183683ebdab5ef74521d1f6a199792a302d48a` before the engine. This work is for independent Opus review; it is not a claim that a four-state projection solves the universe.

The spatial slice is a periodic 2×2×2 cube (8 vertices, 24 *oriented* links and 24 spatial plaquettes). Between slices there are 8 temporal links and 24 temporal plaquettes. Equal plaquette weights are used for all six spacetime orientation pairs. Spatial lengths equal two; the two oppositely directed lattice links between adjacent sites remain **distinct**, as required on a periodic 2-site direction.

The two exact finite gauge groups are **Z₃**, embedded in SU(3) as its center, and the **non-Abelian S₃**, embedded by sign-twisted 3×3 permutation matrices of determinant one. They are finite *subgroups*, not the full continuous SU(3) group. The Boltzmann weights are positive integers: Z₃ `w(e)=4, w(nonidentity)=1`, S₃ `w(e)=5, w(transposition)=1, w(3-cycle)=2`. No floating point or Monte Carlo occurs in this package.

## Transfer rule

At each step take one complete spatial-link configuration `U`, another `V`, and sum the 8 temporal links `g` exactly:

    M(U) = product_{24 spatial plaquettes p} w(hol_p(U))
    K(U,V) = sum_{g in G^8} product_{24 temporal plaquettes p} w(hol_p(U,V,g))
    T(U,V) = M(U) K(U,V)

The asymmetric `T` has a symmetric positive-weight kernel `K` satisfying `T(U,V) / M(U) = T(V,U) / M(V)`, and gauge transformations on U and V leave its entries unchanged. Temporal-variable elimination processes **one time boundary at a time** and reduces the contraction to 12 local pair factors with a maximum intermediate scope of 4 temporal link variables on the tested cube. This is exact, not a truncation.

A four-configuration boundary catalogue is used **only to test 16 kernel entries and 27 jumps**. This catalogue is **not an invariant subspace** and its projected trace and eigenvalues are **not** the finite gauge theory's partition function or mass spectrum.

## Verified results

- 16 catalogue entries for each group, each reconstructed independently through a parity-bipartite contraction.
- Two direct temporal-link sums at each group: 3⁸ = 6,561 and 6⁸ = 1,679,616 terms, matching variable elimination exactly.
- 96 direct four-dimensional plaquette comparisons per group; 2 independent gauge-transformation comparisons per group.
- 27 slice updates and 3, 6, 9, and 27-step matrix jumps with equal exact integer states, 28 chain leaves per group.
- External roots pinned in `ROOT.hbp` and `verify_slices.py`; forged leaf + self-consistent internal chain rejected against the root.
- Six targeted negative controls, including noncommutative multiplication, temporal orientation, plaquette omission, broken class function, broken gauge transform and forged root.
- Measured example: S₃ cold-to-cold temporal factor contraction 8.26 ms vs direct enumeration 1.13 s (local run; hardware-dependent). See `BENCHMARK.hbp` for both tested pairs and both groups.
- Independent exact gauge-orbit count on the full 24-link slice: **Z₃: 129,140,163 states; S₃: 2,821,153,019,713 states.** Direct small-cycle class enumerations and full-cube gauge samples agree. See `orbit_reduce.py` and [`docs/FINITE-SLICE-MATRIX-GUIDE.md`](../../docs/FINITE-SLICE-MATRIX-GUIDE.md).

The S₃ quotient count is still too large for dense eigensolvers. Nothing in this package certifies a static potential, lightest glueball or full-state mass gap. Gauge-invariant spin-network / representation-sector contractions and exact boundary-tensor evolution are the next work items.

## Reproduce

From the repository root:

```sh
p=results/finite-4d-slices-2x2x2-20261106
python3 "$p/finite_slice.py" --out /tmp/finite-slices-verify
for f in RESULT.hbp MATRIX.hbp ROOT.hbp SLICE-CHAIN.hbp; do cmp "$p/$f" "/tmp/finite-slices-verify/$f"; done
python3 "$p/verify_slices.py"
python3 "$p/orbit_reduce.py"
(cd "$p" && sha256sum --check SHA256SUMS)
```

`PREDICTION.hbp` was committed before implementation. `RESULT.hbp`, `MATRIX.hbp`, `SLICE-CHAIN.hbp`, `ROOT.hbp`, `BENCHMARK.hbp`, `ORBIT.hbp`, `SU3-SLICES.hbi` and `SHA256SUMS` carry the reproducibility evidence. `verify_slices.py` is independently written and does not import the engine. Hashes authenticate bytes against the committed root; they do not by themselves establish the authorship of an unpinned record.

## What is still open

1. Build a **closed** gauge-invariant, representation-adapted transfer algebra rather than projecting onto a four-state catalogue.
2. Apply exact tensor-network contractions to matrix-vector products, avoiding enumeration of 129 million or 2.8 trillion orbit states.
3. Prove stable truncation error bounds before attempting approximate low-rank / MPO compression; validate any exact sector against the original plaquette kernel.
4. Extract transfer eigenvalue intervals and static-charge ratios **only after** spectral and charge sectors are proved physical and complete for the intended observable.

No infinite-volume or continuum limit is required for this finite-universe research objective.

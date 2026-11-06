# Matrix laws for the enormous gauge-invariant slice space

**2026-11-06 — finite, discrete, one spatial slice per transfer.**

This extension to `docs/FINITE-SLICE-MATRIX-GUIDE.md` specifies how the complete
2×2×2 `Z3` and non-Abelian `S3` physical state spaces are handled *without*
forming a dense gauge-orbit transfer matrix.

## Law E — gauge projection is an exact tensor operator

Each spatial slice contains 8 vertices and 24 oriented links. The exact
space of local-gauge-invariant functions has 129,140,163 states (`Z3`) or
2,821,153,019,713 states (`S3`). Those counts cannot be mistaken for an
invariant four-state catalogue.

For any class weight `w` and any spatial link configuration `U`, define
`C_e[a,b]=w(a b^-1)`, `C=⊗_{e=1}^{24} C_e`, and `P` as the average over the
`|G|^8` independent vertex gauge transformations. Then

    K = |G|^8 P C = |G|^8 P C P

on gauge-invariant inputs, and the full isotropic finite transfer is

    T = M K,
    M(U)=∏_{24 spatial faces p} w(hol_p(U)).

Because `w` is a class function, `C` commutes with `P`. As `M` is gauge
invariant, `M` commutes with `P`. Hence `T` **preserves the entire physical
Gauss-law subspace** exactly. No sampled catalogue or variational truncation
is needed for this closure theorem. This theorem alone does not determine
the eigenvalues of the interacting `T`.

## Law F — eliminate the next slice before it exists in memory

On a product tensor input `f(V)=∏_e f_e(V_e)`, first sum the input link `v_e`
inside its own convolution factor:

    L_e(g_x,g_y;U_e)=Σ_{v∈G} w(U_e g_y v^-1 g_x^-1) f_e(v).

There remain only 8 temporal vertex variables and 24 factors. On the 2×2×2
cube, pair the two directed edges between each neighboring vertex pair,
then eliminate its four odd-parity vertices for each of the `q^4`
assignments to even-parity vertices. The calculation is **exact**, and its
largest independent input tensor is a length-`q` vector per link.

For a general tensor-network input, apply each local convolution, insert
`P` and multiply by `M`, retaining the expression as a new tensor circuit.
The number of local circuit tensors added per slice is bounded, but the
width of the tensor network needed to compute an observable may grow with
the number of slices. **Do not confuse succinct circuit storage with a
polynomial-time evaluation theorem.**

## Law G — closed flux duality for the entire Z3 magnetic partition

On the periodic 2×2×2 spatial graph, `24` faces map to `24` oriented links
through the mod-3 boundary matrix `∂₂`. Its rank is `14`, so the closed
2-form flux space `ker ∂₂` has dimension `10` and contains `3^10=59049`
elements. The magnetic plaquette weight `w(z)=4` at the identity and `1`
otherwise has exact three-character coefficients `(2,1,1)`.

Summing every spatial link variable analytically projects onto closed
plaquette fluxes. The exact partition is

    Σ_{U∈Z3^24} M(U) = 3^24 Σ_{σ∈ker ∂₂} ∏_p c(σ_p).

The dual sum is `99186762`. A separate C11 primal calculation enumerates
`3^14=4782969` distinct plaquette-flux images, weights each by the exact
number `3^10` of preimages, and gets the same integer. It also reproduces
two Wilson-character numerator sums, independent of the dual algorithm.

This is a genuine all-states calculation, including the full gauge quotient:

    (Σ_U M(U))/3^7 = 12808994612122206.

## Law H — a strictly positive finite interacting gap, not a measured mass

Local link-convolution Fourier eigenvalues are `(6,3,3)` for `Z3` and
`(12,6,3,3,3,3)` for `S3`. Therefore `K` is positive definite on the
physical gauge-invariant subspace. The transfer `T=M K` is similar to the
strictly positive symmetric operator `M^(1/2) K M^(1/2)`.

Furthermore `q^8 ≤ K(U,V) ≤ q^8 wmax^24` for every pair of full spatial
boundaries. Let `r` be the positive Perron eigenvector. The stochastic
Perron transform has entries

    Q(U,V) = K(U,V) r(V) / Σ_W K(U,W) r(W).

All its rows dominate the common probability vector proportional to `r(V)`
with weight at least `wmax^-24`. Consequently every nontrivial transfer
eigenvalue obeys

    0 < λ/λ0 ≤ 1 - wmax^-24.

The same proof passes to orbit coordinates with appropriate orbit-size
weights. It certifies a positive finite gap, but its bound is extremely
conservative and is not the interacting mass or a glueball-channel result.

## Acceptance and remaining work

The exact full-space representations, one-slice action on product tensors,
full `Z3` magnetic partition and Wilson correlators, electric-only spectrum,
and finite positive-gap *inequality* have their own distinct verdicts.

Still OPEN: a practical bounded-width contraction for arbitrary long times;
the interacting numeric eigenvalues and charged sectors on the full 2×2×2
space; complete `S3` interacting partition; and an independent Opus review.
A power jump `T^n` is allowed only after its complete physical state or
contractible tensor representation is retained; power on a nonclosed
four-state catalogue is not accepted as a spectrum calculation.

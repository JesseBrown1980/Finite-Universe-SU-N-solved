# zq8, cloud seat: Q8 lattice gauge theory, first SU(2)-subgroup milestone (Finite-Universe-SU-N-solved).
Pre-registration, frozen before any Q8 code is written or run.

Track A of the SU(2) programme (Q8 -> 2T -> 2O -> 2I). These are finite-subgroup tests, compared with each
other and later with the representation-basis (j_max) calculation of track B; they are not a convergence
sequence to continuous SU(2). No S3 result is inherited. Rust 1.81, clippy -D warnings, integer arithmetic
only, json=0, every computation a job of the HTTP-0 backend (hebc_ym_http0 serve), HBI -> HBP -> SHA -> SH ->
HASH with a separate seal. Claim boundary: exact integer-arithmetic certificates on fixed finite lattices;
the BHM/TQO volume-uniform argument is a separate proof obligation; production HOLD unchanged.

## Group, characters, Hamiltonian
- Q8 is generated from the 2x2 Gaussian-integer matrices i = diag(i, -i), j = [[0, 1], [-1, 0]]; elements by
  closure, multiplication table from exact matrix products. Fundamental character chi_f = trace (values
  2, -2, 0, all integers). The four 1-dimensional characters come from Q8/{+-1}.
- Electric operator, per link: E_link = sum over g in S of (1 - L_g), S = {+-i, +-j, +-k} (the six elements
  at rotation angle pi/2; three conjugacy classes, equal weight 1, closed under inverses).
  Eigenvalue on irrep R: 6 - sum_{g in S} chi_R(g) / d_R = 0 (trivial), 6 (fundamental), 8 (each nontrivial
  1-dimensional irrep). This ordering is not the SU(2) Casimir j(j+1); Q8 coefficients are not read as SU(2)
  approximations.
- H/J_E = E + r M, M = sum over plaquettes (2 - chi_f(U_p)) (values 0, 2, 4). Observable W(C) = chi_f(U_C),
  reported also as W/2 (the (1/N) Re Tr normalisation).
- Geometry: the open 2 x 4 ladder of zwilson (3 plaquettes; loops 1x1, 1x2, 1x3 on columns 0..A);
  couplings r = 1/27, 1/9, 1/3, 1.

## Acceptance checks, in order (each an exact PASS)
- G1 group and characters: closure, associativity, identity, inverses; 8 elements, 5 classes; row and column
  orthogonality of the character table; sum of d_R^2 = 8.
- G2 electric spectrum: the 8 x 8 integer matrix of E_link on L^2(Q8) has characteristic polynomial
  x (x - 6)^4 (x - 8)^3 exactly, and the character formula gives 0, 6, 8.
- G3 Gauss law (independent of G4): for every orbit representative, every vertex v and every h in Q8, the
  locally gauge-transformed raw configuration (tree links included) gauge-fixes into the same orbit, M and
  every W(C) are unchanged, and the electric move targets of the transformed configuration are the same
  orbits with the same multiplicities.
- G4 every matrix written is exactly W-self-adjoint (w_i Z_ij = w_j Z_ji).
- R1 regression of the generic finite-group lane against the sealed S3 results: with the S3 group it
  reproduces the zwilson series c_{A,k} exactly, and one S3 Wilson matrix byte for byte against `wilsonmat`.

## Exact series
<W(C_A)> = sum_k c_{A,k} r^k to order 5, Rayleigh-Schrodinger in H0 = E, V = M with every
<psi_j|W|psi_{k-j}> term and the normalisation, exact rationals:
- method A: gauge-orbit basis under global conjugation with Gauss weights, reduced resolvent by exact
  Gauss-Jordan elimination;
- method B: full gauge-fixed configuration basis (8^3 states, unit weights), reduced resolvent by the
  spectral polynomial over every possible electric eigenvalue (all sums over the 10 links of 0, 6, 8), with
  an exact residual check E x = b.

## Predictions
- W1 series: methods A and B identical for every k <= 5; c_{A,k} = 0 for k < A; and, derived before any code,
  (a) c_{A,A} = 2^-A x (S3 value): c_{1,1} = 1/12, c_{2,2} = 7/1728, c_{3,3} = 11/55296 (at leading order only
      fundamental links with trivial fusions occur, the group factors depend only on d = 2 and orthogonality,
      and every energy denominator doubles: fundamental link energy 6 against 3);
  (b) c_{A,A+1} = 0 for A = 1, 2, 3 (centre symmetry: multiplying the free links by central signs flips every
      plaquette, commutes with E and gauge transformations, and maps H(r) to H(-r) + 4 r N_p, so <W(C_A)>
      contains only r^A, r^{A+2}, ...);
  (c) vacuum energy E_1 = 6, E_2 = -1/8 (N_p s_vac, s_vac = -1/24), E_3 = E_5 = 0.
- W2 certificates: concavity bound (E0(lambda) - E0(0))/lambda <= <W> <= (E0(0) - E0(-lambda))/lambda with
  JE = 1350000000000000, JB = r JE, lambda = 10^-8 (LAM = 13500000); E0(0) from the second-order prediction
  with doubling, E0(+-lambda) from [E0(0) bracket -/+ 2 lambda]; every energy bisected to width <= 1 T unit;
  final brackets re-executed by verify. Every <W(C_A)> enclosed with width <= 2 * 10^-9 and lower end > 0.
- W3 leading order: at r = 1/27, <W(C_A)> within [0.8, 1.25] x c_{A,A} r^A, A = 1, 2, 3.
- W4 convergence: |<W>/(c_{A,A} r^A) - 1| smaller at r = 1/27 than at r = 1/9, every A.
- W5 area law: <W(C_{A+1})> < <W(C_A)> at every r; at r = 1/27 the ratio within [0.8, 1.25] x
  (c_{A+1,A+1}/c_{A,A}) r.
- W6 even corrections at finite coupling: for the 1x1 loop, D(r) = <W>/(r/12) - 1 satisfies
  D(1/9)/D(1/27) in [6, 12] or [54, 108] (an even leading power, r^2 or r^4; the S3 package measured 3.02,
  a linear correction).
Measured, not predicted: c_{A,k} for k >= A + 2, values at r = 1/3 and 1.
Fail-closed: any failed job, missing shard, failed verify, uncertain count or uncertified bracket makes the
affected check INCONCLUSIVE.

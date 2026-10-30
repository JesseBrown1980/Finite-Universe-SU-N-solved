# zt2, cloud seat: binary tetrahedral group 2T, second SU(2)-subgroup milestone (Finite-Universe-SU-N-solved).
Pre-registration, frozen before any 2T code is written or run.

Track A (Q8 -> 2T -> 2O -> 2I): finite-subgroup tests, compared with each other and later with the
representation-basis (j_max) calculation of track B. No S3 or Q8 result is inherited; the Q8 package and its
registered W4 FAIL stay as sealed. Rust 1.81, clippy -D warnings, integer arithmetic only, json=0, every
computation a job of the HTTP-0 server, HBI -> HBP -> SHA -> SH -> HASH with a separate seal.

## Group and exact arithmetic
- 2T = the 24 Hurwitz unit quaternions {+-1, +-i, +-j, +-k, (+-1 +-i +-j +-k)/2}, stored as integer doubled
  coordinates (a, b, c, d), q = (a + b i + c j + d k)/2; multiplication table from exact Hamilton products
  (division by 2 checked exact); closure checked.
- Fundamental character chi_f(q) = 2 Re(q) = a, in {2, -2, 0, 1, -1}.
- Character table in Z[w], w^2 + w + 1 = 0, as exact integer pairs x + y w (conjugation x + y w -> (x - y) - y w):
  trivial; chi_1'(g) = w^k(g), chi_1''(g) = w^(2 k(g)), k(g) in {0, 1, 2} the coset of g in 2T/Q8 relative to
  a = (1 + i + j + k)/2; chi_3(q) = 4 Re(q)^2 - 1; chi_f; chi_2' = chi_f chi_1'; chi_2'' = chi_f chi_1''.
  FGroup characters become Z[w]-valued; S3 and Q8 keep integer values (y = 0).
- Electric operator, per link: E_link = sum over g in S of (1 - L_g), S = the eight elements (1 +- i +- j +- k)/2
  (rotation angle 2 pi / 3; two conjugacy classes of four, each the inverse of the other).
  Eigenvalue E_R = 8 - sum_{g in S} chi_R(g)/d_R: trivial 0, fundamental 4, 3-dimensional 8, 2' and 2'' 10,
  1' and 1'' 12.
- H/J_E = E + r M, M = sum_p (2 - chi_f(U_p)) (values 0, 1, 2, 3, 4). W = chi_f(U_C), also reported as W/2.

## Geometries
- S2: open 2 x 3 strip (2 plaquettes; loops 1x1, 1x2): 576 configurations.
- S3: open 2 x 4 ladder (3 plaquettes; loops 1x1, 1x2, 1x3): 13824 configurations.

## Acceptance checks (each an exact PASS)
- G1: axioms; 7 classes; 7 irreps; row and column orthogonality in Z[w]; sum of d^2 = 24; 1-dimensional
  characters multiplicative.
- G2: the 24 x 24 integer matrix of E_link on L^2(2T) has characteristic polynomial
  x (x - 4)^4 (x - 8)^9 (x - 10)^8 (x - 12)^2 exactly (derived independently from the regular representation
  by the reviewing seat and from the characters here).
- G3: Gauss law on S2 for every orbit representative, every vertex and every element (same orbit, M and W
  invariant, electric moves equivariant).
- G4: every matrix written exactly W-self-adjoint.
- R1 regression after the Z[w] change: the lane reproduces the sealed zwilson (S3) and zq8 (Q8) series exactly,
  and one S3 and one Q8 Wilson matrix byte for byte (against `wilsonmat`, and against the sealed zq8 shard).

## Exact series, to order 5 (every <psi_j|W|psi_{k-j}> term and the normalisation)
- S2: method A gauge-orbit basis with exact Gauss-Jordan elimination; method B configuration basis with the
  spectral polynomial over every possible nonzero electric eigenvalue and exact residual checks (independent
  algorithms and bases).
- S3: method A' gauge-orbit basis with the weighted spectral polynomial (the algorithm validated against
  elimination on S2) and exact residual checks; method B configuration basis as on S2 (two bases).

## Certified finite-coupling values (S2 only; the S3 orbit matrices are too large for sturmrs in replay time)
Concavity bound (E0(lambda) - E0(0))/lambda <= <W> <= (E0(0) - E0(-lambda))/lambda, JE = 1350000000000000
(27 x 5 * 10^13), JB = r JE, lambda = 10^-8 (LAM = 13500000), r = 1/27, 1/9, 1/3, 1; E0(0) from the
second-order prediction 2 r N_p - N_p r^2 / 16 with doubling, E0(+-lambda) from [E0(0) bracket -/+ 2 lambda];
each energy bisected to width <= 1 T unit (7.4 * 10^-19 J_E), final brackets re-executed by verify.
At r = 1 the largest threshold, one doubling included, is 8.6 * 10^18, inside the 64-bit spec field.

## Predictions
- W1 series (S2 and S3): the two methods identical for every k <= 5; c_{A,k} = 0 for k < A; derived before any
  code,
  (a) c_{A,A} = (3/4)^A x S3: c_{1,1} = 1/8, c_{2,2} = 7/768, c_{3,3} = 11/16384 (normalised W/2: 1/16, 7/1536,
      11/32768), the same on S2 and S3 (fundamental link energy 4 against S3's 3; d = 2; 2 x 2 = 1 + 3);
  (b) c_{A,A+1} = 0 (the centre -1 flips every plaquette; H(r) -> H(-r) + 4 r N_p);
  (c) vacuum energy E_1 = 2 N_p, E_2 = -N_p/16, E_3 = E_5 = 0.
- W2 certificates (S2): every <W> enclosed with width <= 2.5 * 10^-9 and lower end > 0.
- W3 leading order (S2): at r = 1/27, <W(C_A)> within [0.8, 1.25] x c_{A,A} r^A, A = 1, 2.
- W4 convergence with a resolution rule (the Q8 lesson): for each loop let delta = |c_{A,A+2}/c_{A,A}| (1/9)^2,
  the series-predicted relative deviation at r = 1/9 from the in-package exact series. The loop is resolved if
  the certified ratio intervals <W>/(c_{A,A} r^A) at r = 1/27 and r = 1/9 both have width below delta; for a
  resolved loop |ratio - 1| must be smaller at r = 1/27 than at r = 1/9. W4 = PASS if at least one loop is
  resolved and every resolved loop passes; an unresolved loop is INCONCLUSIVE by resolution, not FAIL.
- W5 area law (S2): <W(1x2)> < <W(1x1)> at every r; at r = 1/27 the ratio within [0.8, 1.25] x
  (c_{2,2}/c_{1,1}) r.
- W6 even correction (S2, 1x1): D(1/9)/D(1/27) in [6, 12] or [54, 108], D(r) = <W>/(r/8) - 1.
Measured, not predicted: c_{A,k} for k >= A + 2; S3 1x3 coefficients beyond leading order; values at r = 1/3, 1.
Fail-closed: any failed job, missing shard, failed verify, uncertain count or uncertified bracket makes the
affected check INCONCLUSIVE.

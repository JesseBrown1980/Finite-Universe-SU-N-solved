# z2o, cloud seat: binary octahedral group 2O, third SU(2)-subgroup milestone, open 2 x 3 strip
# (Finite-Universe-SU-N-solved). Pre-registration, frozen before any 2O code or the sturmq2 counter exists.

Track A (Q8 -> 2T -> 2O -> 2I). No earlier result is inherited; Q8's registered W4 FAIL and every seal stay as
they are. Rust 1.81, clippy -D warnings, exact integer and quadratic-field arithmetic only (no floats),
json=0, every computation a job of the HTTP-0 server, HBI -> HBP -> SHA -> SH -> HASH with a separate seal.
The 2 x 4 ladder is a later package, gated on a solver-scaling milestone (4888 orbit states; candidate
spectral degree 582 under this electric set).

## Group and exact arithmetic
- 2O = the 48 unit quaternions: the 24 of 2T and the 24 elements (+-u +-v)/sqrt2 for distinct u, v in
  {1, i, j, k}. Coordinates stored exactly as (x + y sqrt2)/2 with integers x, y; multiplication table from
  exact Hamilton products in Z[sqrt2] (division by 2 checked exact); closure checked.
- Fundamental character chi_f(q) = 2 Re(q) in Z[sqrt2], values 2, -2, 0, +-1, +-sqrt2.
- Characters in Z[sqrt2] (8 irreps): trivial; sign (+1 on 2T, -1 on 2O \ 2T); chi_2 (2 on Q8, -1 on 2T \ Q8,
  0 on 2O \ 2T); chi_3 = 4 Re(q)^2 - 1; chi_3' = chi_3 sign; chi_f; chi_f' = chi_f sign; chi_4 = chi_f chi_2.
- Electric operator, per link: E_link = sum over g in S of (1 - L_g), S = {(1 +- i)/sqrt2, (1 +- j)/sqrt2,
  (1 +- k)/sqrt2} (the six elements nearest the identity, rotation angle pi/2, one conjugacy class, closed
  under inverses). The matrix of E is an integer matrix; its eigenvalues E_R = 6 - sum_{g in S} chi_R(g)/d_R
  are 0 (trivial), 6 - 3 sqrt2 (fundamental, the lowest excitation), 4 (3), 6 (2 and 4), 8 (3'),
  6 + 3 sqrt2 (fundamental'), 12 (sign).
- H/J_E = E + r M, M = sum_p (2 - chi_f(U_p)), diagonal with entries in Z[sqrt2]. W = chi_f(U_C), also W/2.
- Geometry: the open 2 x 3 strip (2 plaquettes; loops 1x1, 1x2): 2304 configurations, 136 gauge orbits.

## Acceptance checks (each an exact PASS)
- G1: axioms; 8 classes; 8 irreps; row and column orthogonality in Z[sqrt2]; sum of d^2 = 48; sign and chi_2
  characters of homomorphisms (sign multiplicative).
- G2: the 48 x 48 integer electric matrix has characteristic polynomial
  x (x - 4)^9 (x - 6)^20 (x - 8)^9 (x - 12) (x^2 - 12 x + 18)^4 exactly (independently factorised by the
  reviewing seat), and the character formula gives the spectrum above.
- G3: Gauss law on the 2 x 3 strip for every orbit representative, vertex and element.
- G4: every matrix exactly W-self-adjoint (Z[sqrt2] diagonal, integer electric part).
- R1: with quadratic-field arithmetic, the lane reproduces the sealed S3 (zwilson), Q8 (zq8) and 2T (zt2)
  series exactly, and the S3, Q8 and 2T integer matrices byte for byte.
- S1 sturmq2 regression: the new counter (exact LDL^T over Q(sqrt2)) gives, on the regenerated sealed integer
  matrices, exactly the sturmrs negative counts at the final bracket thresholds of every E0(0) certificate in
  zwilson, zq8 and zt2 (count 0 at the lower end, count >= 1 at the upper end).

## Methods
- Exact series to order 5 over Q(sqrt2) (vectors as pairs p + q sqrt2 with rational p, q): method A gauge-orbit
  elimination (the reduced resolvent of the integer E is rational); method B configuration basis with the
  conjugation-closed candidate spectral polynomial (integer coefficients), each solve checked exactly
  (E x = b, <1, x>_W = 0, <1, b>_W = 0).
- sturmq2: inertia of the exactly symmetric B(T) = W (1000 (JE E + JB M + LAM W_C) - T W) over Q(sqrt2) by
  LDL^T with symmetric pivoting on a nonzero diagonal entry, an explicit 2 x 2 pivot block when every remaining
  diagonal entry is zero, and the exact sign test of a + b sqrt2 (comparing a^2 with 2 b^2 when the signs
  differ). No intervals and no floats.
- Certified values: concavity bound (E0(lambda) - E0(0))/lambda <= <W> <= (E0(0) - E0(-lambda))/lambda with
  JE = 810000000000000 (27 x 3 * 10^13), JB = r JE, lambda = 4 * 10^-9 (LAM = 3240000), r = 1/27, 1/9, 1/3, 1;
  E0(0) bracketed from 2 r N_p - N_p r^2 (2 + sqrt2)/24 with sqrt2 replaced by 99/70 for the starting point
  only, doubling as before; E0(+-lambda) from [E0(0) bracket -/+ 2 sqrt2 lambda] rounded outward (|W| <= 2);
  every energy bisected to width <= 1 T unit; final brackets re-executed by verify.

## Predictions
- W1 series: methods A and B identical for every k <= 5; c_{A,k} = 0 for k < A; derived before any code from the
  leading law of the sealed S3, Q8 and 2T series (c_{1,1} = 1/(2 E_f), c_{2,2} = 7/(48 E_f^2)) with
  E_f = 6 - 3 sqrt2:
  (a) c_{1,1} = (2 + sqrt2)/12, c_{2,2} = (21 + 14 sqrt2)/864 (normalised W/2: halves);
  (b) c_{A,A+1} = 0 (centre -1);
  (c) vacuum energy E_1 = 4, E_2 = -(2 + sqrt2)/12, E_3 = E_5 = 0.
- W2 certificates: every <W(C_A)> enclosed with width <= 3 * 10^-9 and lower end > 0.
- W3 leading order: at r = 1/27, <W(C_A)> within [0.8, 1.25] x c_{A,A} r^A (exact comparison in Q(sqrt2)).
- W4 convergence with the zt2 resolution rule (delta = |c_{A,A+2}/c_{A,A}| (1/9)^2 from the in-package series;
  unresolved loops INCONCLUSIVE by resolution; PASS needs at least one resolved loop and every resolved loop
  converging).
- W5 area law: <W(1x2)> < <W(1x1)> at every r; at r = 1/27 the ratio within [0.8, 1.25] x (c_{2,2}/c_{1,1}) r.
- W6 even correction: D(1/9)/D(1/27) in [6, 12] or [54, 108], D(r) = <W(1x1)>/(c_{1,1} r) - 1.
Measured, not predicted: c_{A,k} for k >= A + 2; values at r = 1/3 and 1.
Fail-closed: any failed job, missing shard, failed verify, uncertain count or uncertified bracket makes the
affected check INCONCLUSIVE.

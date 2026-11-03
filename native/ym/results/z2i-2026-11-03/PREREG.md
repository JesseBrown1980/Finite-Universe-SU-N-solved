# z2i, acer seat: binary icosahedral group 2I, fourth SU(2)-subgroup milestone, open 2 x 3 strip
# (Finite-Universe-SU-N-solved). Pre-registration, frozen before any 2I code.

Track A (Q8 -> 2T -> 2O -> 2I), its last finite subgroup. No earlier result is inherited; Q8's registered W4 FAIL, the qutrit
network's N4 FAIL and every seal stay as they are. Rust 1.81, clippy -D warnings, exact integer and quadratic-field arithmetic only
(no floats), json=0, every computation a job of the HTTP-0 server, HBI -> HBP -> SHA -> SH -> HASH with a separate seal.
The 2 x 4 ladder is a later package, gated on the same solver-scaling milestone as 2O's: 29288 orbit states for 2I, against 4888 for 2O.

## Group and exact arithmetic
- phi = (1 + sqrt5)/2, phi^2 = 1 + phi, phi^-1 = phi - 1. Z[phi] and Q(phi) use the existing Field with p = 1, q = 1.
- 2I = the 120 unit icosians: the 24 of 2T, then the 96 even permutations of the coordinates of (+-phi, +-1, +-phi^-1, 0)/2.
  Coordinates stored doubled in Z[phi]: q = (X + Y i + Z j + U k)/2. Multiplication table from exact Hamilton products (division by 2
  checked exact); closure checked.
- Fundamental character chi_f(q) = 2 Re(q) = X in Z[phi], with values 2, -2, 0, +-1, +-phi, +-phi^-1.
- Classes, 9, by chi_f: 2 (1 element), -2 (1), 0 (30), 1 (20), -1 (20), phi (12), -phi (12), phi^-1 (12), -phi^-1 (12).
- Characters in Z[phi], 9 irreps, with t = chi_f and sigma the conjugation phi -> 1 - phi on values:
  chi_1 = 1; chi_2 = t; chi_3 = t^2 - 1; chi_4 = t^3 - 2t; chi_5 = t^4 - 3t^2 + 1; chi_6 = t^5 - 4t^3 + 3t;
  chi_2' = sigma(chi_2); chi_3' = sigma(chi_3); chi_4' = chi_2 chi_2'.
- Electric operator, per link: E_link = sum over g in S of (1 - L_g), S = the 12 elements with chi_f = phi (rotation angle 2 pi/5, the
  elements nearest the identity, one conjugacy class, closed under inverses). The matrix of E is an integer matrix. Its eigenvalues
  E_R = 12 - 12 chi_R(g)/d_R, g in S, are 0 (chi_1), 9 - 3 sqrt5 = 12 - 6 phi = 6/phi^2 (chi_2, the lowest excitation, 2.2917960675),
  10 - 2 sqrt5 (chi_3), 9 (chi_4), 12 (chi_5), 14 (chi_6), 15 (chi_4'), 10 + 2 sqrt5 (chi_3'), 9 + 3 sqrt5 (chi_2').
- H/J_E = E + r M, M = sum_p (2 - chi_f(U_p)), diagonal with entries in Z[phi]. W = chi_f(U_C), also W/2.
- Geometry: the open 2 x 3 strip (2 plaquettes; loops 1x1, 1x2): 14400 configurations, 296 gauge orbits.
  Derived: orbits = (1/|G|) sum_g |C(g)|^2 = (2 x 14400 + 30 x 16 + 40 x 36 + 48 x 100)/120 = 296. The same formula gives the recorded
  11, 28, 76 and 136 for S3, Q8, 2T and 2O, and with cubes gives 49, 176 and 1312 on 2 x 4, the 4888 recorded for 2O in its PREREG, and 29288.

## Acceptance checks (each an exact PASS)
- G1: axioms; 9 classes with the sizes above; 9 irreps; row and column orthogonality in Z[phi]; sum of d^2 = 120.
- G2: the 120 x 120 integer electric matrix has characteristic polynomial
  x (x - 9)^16 (x - 12)^25 (x - 14)^36 (x - 15)^16 (x^2 - 18 x + 36)^4 (x^2 - 20 x + 80)^9 exactly, and the character formula gives the
  spectrum above. Derived: degree 120, trace 1440, trace of the square 18720, each equal to its direct count from E_link.
- G3: Gauss law on the 2 x 3 strip for every orbit representative, vertex and element: 296 x 6 x 120 = 213120 local gauge transformations.
- G4: every matrix exactly W-self-adjoint (Z[phi] diagonal, integer electric part).
- R1: with the phi lane added, the sealed S3 (zwilson), Q8 (zq8), 2T (zt2) and 2O (z2o) series reproduce exactly, and their matrices byte for
  byte. The sealed packages' own replays stay in CI on the exact head.
- S1: interval sturmq2 reproduces the sturmrs counts at every final E0(0) bracket threshold of zwilson, zq8 and zt2 (0 at the lower end, >= 1
  at the upper end), and exact sturmq2 does on zwilson and zt2. The exact Q8 counts and z2o's X1 stay with their sealed replays.
- S2: the exact sign test of a + b phi (comparing (2a + b)^2 with 5 b^2 when the signs of 2a + b and b differ) agrees with its rational
  enclosure on every pair with |a|, |b| <= 64 and on consecutive Fibonacci pairs up to 2^4096, where a + b phi is smallest in size.

## Methods
- Exact series to order 5 over Q(phi): method A gauge-orbit elimination (the reduced resolvent of the integer E is rational); method B
  configuration basis with the conjugation-closed candidate spectral polynomial (integer coefficients), each solve checked exactly.
- sturmq2: exact mode as in z2o, with the sign test of S2. Interval mode as in z2o addendum 1, with phi in [L, L + 1]/10^80,
  L = floor((10^80 + isqrt(5 * 10^160))/2), checked exactly before counting (L^2 - 10^80 L - 10^160 <= 0 < (L + 1)^2 - 10^80 (L + 1) - 10^160) and
  recorded. A count is accepted only with status=COMPLETE, anything else is INCONCLUSIVE.
- Certified values: concavity bound (E0(lambda) - E0(0))/lambda <= <W> <= (E0(0) - E0(-lambda))/lambda with JE = 810000000000000,
  JB = r JE, lambda = 4 * 10^-9 (LAM = 3240000), r = 1/27, 1/9, 1/3, 1; E0(0) bracketed from 2 r N_p - N_p r^2 phi^2/24 with phi replaced by
  987/610 for the starting point only, doubling as before; E0(+-lambda) from [E0(0) bracket -/+ 2 phi lambda] rounded outward (|W| <= 2); every
  energy bisected to width <= 1 T unit; final brackets re-executed by verify.

## Gates, before any 2I count, series job or official run
Timing runs only, at thresholds outside every registered window. Each outcome goes in a sealed addendum before the run it governs:
- X1 gate. 2O took 80.5 and 81.8 min for its two exact counts at 136 states. Scaling that to 296 states projects tens of hours per count
  (CONJECTURE; the gate measures it). If one exact 2I count at the final thresholds projects beyond a 6-hour hosted job, X1 runs in exact and interval
  mode at coarser thresholds named in the addendum, and that scope is stated beside every X1 result.
- Method B gate. 2O's candidate spectral degree was 282 at 2304 configurations. 2I has 14400. If order 5 by method B projects beyond a 6-hour
  hosted job, method B stops at the highest order that fits, named in the addendum, and W1 states it. Method A still runs to order 5.

## Predictions
- X1: on the 2I matrix at r = 1/27 with lambda = 0, exact sturmq2 gives the same counts as interval mode at the thresholds the gate sets.
- W1 series: methods A and B identical at every order both reach; c_{A,k} = 0 for k < A; from the leading law of the four sealed series
  (c_{1,1} = 1/(2 E_f), c_{2,2} = 7/(48 E_f^2)) with E_f = 6/phi^2:
  (a) c_{1,1} = phi^2/12 = (1 + phi)/12 = (3 + sqrt5)/24 = 0.2181694991; c_{2,2} = 7 phi^4/1728 = (14 + 21 phi)/1728 = (49 + 21 sqrt5)/3456 = 0.0277654593
      (normalised W/2: halves);
  (b) c_{A,A+1} = 0 (centre -1);
  (c) vacuum energy E_1 = 4, E_2 = -phi^2/12, E_3 = E_5 = 0; and c_{1,1} = -E_2, c_{1,3} = -2 E_4.
- W2 certificates: every <W(C_A)> enclosed with width <= 3 * 10^-9 and lower end > 0.
- W3 leading order: at r = 1/27, <W(C_A)> within [0.8, 1.25] x c_{A,A} r^A (exact comparison in Q(phi)).
- W4 convergence with the zt2 resolution rule (delta = |c_{A,A+2}/c_{A,A}| (1/9)^2 from the in-package series; unresolved loops INCONCLUSIVE by
  resolution; PASS needs at least one resolved loop and every resolved loop converging).
- W5 area law: <W(1x2)> < <W(1x1)> at every r; at r = 1/27 the ratio within [0.8, 1.25] x (c_{2,2}/c_{1,1}) r.
- W6 even correction: D(1/9)/D(1/27) in [6, 12] or [54, 108], D(r) = <W(1x1)>/(c_{1,1} r) - 1.
Measured, not predicted: c_{A,k} for k >= A + 2; values at r = 1/3 and 1; the candidate spectral degree; the ratios E_R/E_f against the SU(2) Casimir
ratios j(j + 1)/(3/4).
Fail-closed: any failed job, missing shard, failed verify, uncertain count or uncertified bracket makes the affected check INCONCLUSIVE.
Any prediction that misses is recorded as a FAIL and kept.

## Scope
This package is a statement about one finite group on one open strip of 2 plaquettes. Continuous SU(2), 3+1 dimensions and continuum scaling are later steps of
the sequence in the README.

# zspec, cloud seat: exact characteristic polynomials of the finite-group strip Hamiltonians
# (Finite-Universe-SU-N-solved). Pre-registration, frozen before any zspec code.

A third route to the strong-coupling series, and exact algebraic values at finite coupling: the characteristic
polynomial det(e I - H) of the gauge-orbit Hamiltonians H/J_E = E + r M + lambda W_C of the sealed packages,
computed exactly. Registered after the review of van Vliet's tameness thesis (the eigenvalue branches of a
finite lattice are semialgebraic), with the reviewing seat's refinements: the lambda-derivative taken before
lambda = 0, the factor carrying E0 identified instead of reading the full discriminant, Kato's bound in the
weighted space, and total degree and coefficient heights recorded. The sealed packages zwilson (S3), zq8 (Q8),
zt2 (2T) and z2o (2O) are immutable regression anchors: every result registered there, zq8's W4 FAIL included,
stays as registered. Rust 1.81, clippy -D warnings, exact integer, modular and quadratic-field arithmetic only
(no floats), json=0, every computation a job of the HTTP-0 server, HBI -> HBP -> SHA -> SH -> HASH with a
separate seal.

## Systems
E, M, W_C and the orbit weights from fgroup::basis, the builder of the sealed packages (E integer; M and W_C
diagonal, in Z or Z[sqrt2]).

| system | group | open strip | d (orbit states) | loops C | sealed reference |
| --- | --- | --- | --- | --- | --- |
| S3/2x3 | S3 | 2 x 3 | 11 | 1x1, 1x2 | none |
| S3/2x4 | S3 | 2 x 4 | 49 | 1x1, 1x2, 1x3 | zwilson: series to order 4, 12 certified values |
| Q8/2x3 | Q8 | 2 x 3 | 28 | 1x1, 1x2 | none |
| Q8/2x4 | Q8 | 2 x 4 | 176 | 1x1, 1x2, 1x3 | zq8: series to order 5, 12 certified values |
| 2T/2x3 | 2T | 2 x 3 | 76 | 1x1, 1x2 | zt2: 2 x 3 series to order 5, 8 certified values |
| 2O/2x3 | 2O | 2 x 3 | 136 | 1x1, 1x2 | z2o: series to order 5, 8 certified values |

## Objects
- P0(e, r) = det(e I - E - r M), in Z[e, r] (Z[sqrt2][e, r] for 2O): monic of degree d in e, total degree d.
- P1_C(e, r) = the coefficient of lambda in det(e I - E - r M - lambda W_C), computed with lambda^2 = 0 (dual
  numbers), never by setting lambda = 0 first.
- The ground branch E0(r): the root of P0 with E0(0) = 0 (the electric vacuum, simple at r = 0); by
  Feynman-Hellmann <W_C>(r) = -P1_C(E0, r) / (d/de) P0(E0, r).
- mu(e, r): the minimal polynomial over Q(sqrt2)(r) of the electric vacuum v (the constant function on gauge
  orbits, the kernel of E) under E + r M. v and the ground state are positive in the configuration basis, so
  E0(r) is a root of mu, and mu divides P0.

## Algorithm (registered before any code)
- Primes: the largest primes below 2^62 with p = 1 mod 8, in decreasing order. For 2O the two embeddings
  sqrt2 -> +s, -s (s by Tonelli-Shanks) give x and y of x + y sqrt2 mod p.
- Coefficient bound B: the largest coefficient of the majorant prod_k (e + sum_l |E_kl| + r m_k + lambda w_k),
  with m_k, w_k integer bounds |x| + 2|y| of both embeddings of M_kk, W_kk (permanent bound of the Leibniz
  expansion). The number of primes is the least k with prod p_i > 2 B + 1; CRT with symmetric residues. No
  stopping on agreement.
- Per prime and r0 = 0, 1, ..., d: det(e I - E - r0 M - lambda W_C) over GF(p)[lambda]/(lambda^2) by reduction
  to upper Hessenberg form (similarity transforms; pivots with a nonzero lambda-free part; a prime where a
  column has nonzero entries but no such pivot is skipped and recorded) and the Hessenberg recurrence for the
  characteristic polynomial; Lagrange interpolation in r on 0..d.
- Independent evaluation: Gaussian elimination over Q(sqrt2) with reduced rationals (as sturmq2 exact), giving
  det A and d/dlambda det(A - lambda W_C) at lambda = 0 = -sum_k (W_C)_kk adj(A)_kk for A = e0 I - E - r0 M.
- mu: Krylov ranks of [v, Hv, H^2 v, ...] modulo a prime at fixed integer r0 bound deg mu from below; the
  candidate mu is reconstructed modularly and is accepted only when mu(E + r M, r) v = 0 holds exactly at 2m + 1
  integer values of r and mu divides P0 exactly.

## Checks (each an exact PASS)
- K1 exactness: P0 monic of degree d with total degree d and the coefficient of e^(d-1) equal to
  -Tr(E + r M); at the test points (e0, r0) = (-1, 0), (-1, 1), (-1/2, 1/3), (-5/3, 2) (e0 < 0 and H >= 0, so
  A is invertible) P0 and every P1_C equal the values of the independent elimination.
- K2 series: P0(0, 0) = 0 and (d/de) P0(0, 0) != 0 exactly; Newton iteration in Q(sqrt2)[[r]] modulo r^6 gives
  E0(r) = sum_k E_k r^k and <W_C>(r) = sum_k c_{C,k} r^k. Every E_k and c_{C,k} identical to the sealed lines
  (methods A/A' and B): S3/2x4 through order 4 (zwilson), Q8/2x4 (zq8), 2T/2x3 (zt2) and 2O/2x3 (z2o) through
  order 5.
- K3 values: for each of the 40 sealed certificates, the root of P0(e, r0) in the sealed E0(0) bracket
  [T_lo, T_hi] / (1000 J_E) is the only one there (Descartes' rule on the Moebius image, exact when the count
  is 0 or 1), refined by exact sign evaluation to width <= 10^-60; <W_C>(r0) is enclosed by interval division
  of the ranges of P1_C and (d/de) P0 on the final bracket (each monotone there, certified by Descartes on its
  derivative, so its range is spanned by its endpoint values; the denominator range excludes 0), with width
  <= 10^-40; every enclosure lies inside the sealed [W_lo, W_hi].
- K4 W4 resolved: with the K3 enclosures, for the ten loops of the four sealed systems, ratio(r) =
  <W_C>(r) / (c_{A,A} r^A) at r = 1/27 and 1/9. Resolved when one of the two strict inequalities between the
  enclosures of |ratio(1/27) - 1| and |ratio(1/9) - 1| holds; converging when it is
  max |ratio(1/27) - 1| < min |ratio(1/9) - 1|. A new measurement; zq8's W4 FAIL stays as registered.
- K5 ground factor: m = deg_e mu certified from both sides (Krylov rank at r0 = 7, 11, 13 modulo the first
  three primes; the exact identities above). Irreducibility over Q(sqrt2)(r) certified when mu(e, r0) modulo p
  (sqrt2 -> s) is squarefree and irreducible of degree m by distinct-degree factorisation, for some r0 in
  1..8 and one of the first 8 primes; otherwise INCONCLUSIVE. Recorded: m, deg_r mu, total degree, maximum
  coefficient bits, and the conditions identifying the branch (mu(0, 0) = 0, (d/de) mu(0, 0) != 0).
- K6 Kato radius: g = the smallest positive root of P0(e, 0) (isolated exactly), c = (max_k M_kk + min_k M_kk)/2,
  a = max_k |M_kk - c| (the norm of M - c I in the weighted space, M diagonal), R_K = g / (2 a): the ground
  eigenvalue and every <W_C> are analytic for |r| < R_K.

## Predictions
- K1, K2: exact identity throughout. On the two 2 x 3 systems without sealed series, from the laws note
  (docs/FINITE_GROUP_LAWS.md): c_{1,1} = 1/(2 E_f) and c_{2,2} = 7/(48 E_f^2), that is 1/6 and 7/432 (S3),
  1/12 and 7/1728 (Q8); E_1 = 4; for Q8 c_{A,A+1} = 0 and E_3 = E_5 = 0. On all four 2 x 3 systems the strip
  identity <W(1x1)> = 2 - (1/2) dE0/dr holds order by order: c_{1,k} = -(k+1) E_{k+1} / 2 for k <= 4.
- K3: all 40 exact values inside their sealed intervals.
- K4: all ten loops resolved and converging, Q8's 1x3 loop included.
- K5: mu irreducible over Q(sqrt2)(r) and m < d (the strip reflections fix v) for all six systems.
- K6: g = 4 E_f with E_f = 3, 6, 4, 6 - 3 sqrt2 for S3, Q8, 2T, 2O (one fundamental loop around a plaquette),
  and R_K > 1/3 for all six systems: R_K = 2 (S3/2x3), 4/3 (S3/2x4), 3 (Q8/2x3), 2 (Q8/2x4), 2 (2T/2x3),
  3 - (3/2) sqrt2 (2O/2x3).
Measured, not predicted: B and the coefficient heights of P0 and P1_C; the series coefficients beyond those
predicted (order 5 for S3/2x4 included); m, deg_r mu, total degree and coefficient bits; the 40 values to 40
digits. In a later package: the branch points of mu (its discriminant in e) and the exact radius; 2T on the
2 x 4 ladder (1312 states); 2I.
Fail-closed: any failed job, missing shard, failed verify, unmet prime bound, Descartes count other than 1,
enclosure wider than registered, failed exact identity or missing sealed line makes the affected check
INCONCLUSIVE. Registered FAILs are kept.

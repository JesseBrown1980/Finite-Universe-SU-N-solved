# Finite-group laws and the sizes of the next packages (2026-10-31)

Exact relations shared by the sealed finite-group packages (zwilson S3, zq8 Q8, zt2 2T), the predictions they
give for the 2O 2×4 ladder and for 2I, and the exact counts that size those packages. Recorded before any 2O
2×4 or 2I code.

## 1. Leading Wilson-loop coefficients follow the fundamental link energy

E_f is the electric energy of the fundamental irrep on one link, E_R = |S| − Σ_{g∈S} χ_R(g)/d_R, for the
registered electric set S of each group. W = χ_f(U_C), the unnormalised fundamental loop.

| Group | S | E_f | c_{1,1} | c_{2,2} | c_{3,3} |
| --- | --- | --- | --- | --- | --- |
| S3 | the 3 transpositions | 3 | 1/6 | 7/432 | 11/6912 |
| Q8 | ±i, ±j, ±k | 6 | 1/12 | 7/1728 | 11/55296 |
| 2T | the 8 elements (1 ± i ± j ± k)/2 | 4 | 1/8 | 7/768 | 11/16384 |

All nine sealed coefficients are exactly

    c_{1,1} = 1/(2 E_f),   c_{2,2} = 7/(48 E_f²),   c_{3,3} = 11/(256 E_f³).

At leading order a 1×A loop puts only single fundamental quanta on its links, so every energy denominator is a
whole multiple of E_f and the group enters through E_f alone (every fundamental here is 2-dimensional); the
prefactors 1/2, 7/48 and 11/256 are the geometry of the strip. The groups first differ at the next order, where
other irreps and their energies enter: the first correction of the 1×1 loop on the 2×4 ladder is +185/456192
for Q8 and −5/12288 for 2T.

Predictions from the law (exact in Q(√2) and Q(√5)):

| Group | S | E_f | c_{1,1} | c_{2,2} | c_{3,3} |
| --- | --- | --- | --- | --- | --- |
| 2O | the 6 elements (1 ± i)/√2, (1 ± j)/√2, (1 ± k)/√2 | 6 − 3√2 = 1.7573593129 | (2 + √2)/12 = 0.2845177969 | (21 + 14√2)/864 = 0.0472210531 | (110 + 77√2)/27648 = 0.0079171891 |
| 2I | the 12 elements of real part φ/2 | 12 − 6φ = 9 − 3√5 = 2.2917960675 | (3 + √5)/24 = 0.2181694991 | (49 + 21√5)/3456 = 0.0277654594 | (99 + 44√5)/55296 = 0.0035696432 |

The 2O c_{1,1} and c_{2,2} are the registered W1 predictions of z2o (`PREREG.md`, `a5c0b994…`); c_{3,3} needs the
2×4 ladder. For 2I, the 12 elements of real part φ/2 form one conjugacy class, closed under inverses, and are
the elements nearest the identity, so χ_f = φ on S. The fundamental link energies are not monotone along
Q8 → 2T → 2O → 2I (6, 4, 1.757, 2.292): they depend on the size and rotation angle of each electric set,
E_f = |S| (1 − cos(θ/2)).

## 2. The 1×1 loop and the vacuum energy on the 2×3 strip

M = Σ_p (2 − χ_f(U_p)), and the two plaquettes of the open 2×3 strip are equivalent, so Feynman–Hellmann gives
dE₀/dr = ⟨M⟩ = 4 − 2⟨W(1×1)⟩, that is

    ⟨W(1×1)⟩ = 2 − ½ dE₀/dr,   c_{1,k} = −(k + 1) E_{k+1} / 2.

The sealed 2T series satisfies it exactly: E₀ = 4r − r²/8 + 5r⁴/24576 + O(r⁶) and c_{1,1} = 1/8 = −E₂,
c_{1,3} = −5/12288 = −2E₄. On the 2×3 strip the first correction of the 1×1 loop is therefore fixed by the r⁴
coefficient of the vacuum energy, with the opposite sign. On longer ladders the identity holds for the
plaquette average, Σ_p ⟨W_p⟩ = 2N_p − dE₀/dr.

## 3. Exact counts

By Burnside's lemma an open strip with k free group elements after gauge fixing (k = 2 on 2×3, k = 3 on 2×4) has
N_k = Σ_{[g]} |C(g)|^{k−1} gauge orbits. Flat connections on a periodic torus are commuting holonomy pairs up to
conjugation, N_flat = Σ_{[g]} k(C(g)), with k(C) the number of conjugacy classes of the centraliser.

Class data (class size × number of classes, centraliser, its order and class number):

| Group | Central classes | Non-central classes |
| --- | --- | --- |
| S3 | 1: S3, 6, 3 | 3: Z2, 2, 2; 2: Z3, 3, 3 |
| Q8 | 1 × 2: Q8, 8, 5 | 2 × 3: Z4, 4, 4 |
| 2T | 1 × 2: 2T, 24, 7 | 6: Z4, 4, 4; 4 × 4: Z6, 6, 6 |
| 2O | 1 × 2: 2O, 48, 8 | 6 (±i, ±j, ±k): Z8, 8, 8; 8 × 2: Z6, 6, 6; 6 × 2 (order 8): Z8, 8, 8; 12 (order 4 outside 2T): Z4, 4, 4 |
| 2I | 1 × 2: 2I, 120, 9 | 30: Z4, 4, 4; 20 × 2: Z6, 6, 6; 12 × 4: Z10, 10, 10 |

| Group | 2×3 open orbits | 2×4 open orbits | periodic flat vacua |
| --- | --- | --- | --- |
| S3 | 11 | 49 | 8 |
| Q8 | 28 | 176 | 22 |
| 2T | 76 | 1312 | 42 |
| 2O | 136 | 4888 | 56 |
| 2I | 296 | 29288 | 74 |

Every non-central centraliser is abelian, so N_flat = N_2 − Σ_{central} (|G| − k(G)). The open counts are the
orbit bases of the packages (S3 49 and Q8 176 on 2×4, 2T 76 and 1312, 2O 136 on 2×3), and S3's 8 is the
flat-vacuum count of every periodic S3 torus in The-matrix-law-essentials. All entries were also enumerated
directly with exact quaternion arithmetic (coordinates in Z[√2]/2 for 2O and Z[√5]/4 for 2I, closure checked),
the 2I 2×4 count from the Burnside sum. The 2I 2×3 strip has 296 orbit states and 14400 configurations; the
2O 2×4 ladder 4888 and 110592; the 2I 2×4 ladder 29288 and 1728000.

## 4. Exact counts by reflection sectors

The exact inertia counts at λ = 0 (S1 against the sealed brackets, X1, and any exact E₀(0) certificate) are the
expensive ones. At λ = 0 the left–right and up–down reflections of an open strip commute with H: they map gauge
orbits to gauge orbits, S is closed under inverses and conjugation, and χ_f is real. A rational change of basis
to reflection-even and reflection-odd combinations makes B(T) block diagonal with the same inertia (Sylvester's
law), so one exact count becomes up to four counts on smaller blocks; with entry sizes growing with the dimension
being eliminated, the saving exceeds the dimension ratio. A Wilson source on one plaquette breaks the left–right
reflection, so the λ ≠ 0 counts keep the full basis or the up–down sectors only. This is a method for a future
pre-registration (2O 2×4, 2I), registered with its sector checks (commutation with H on every orbit, exact
W-self-adjointness per block, block dimensions summing to the full basis).

# zt2: binary tetrahedral group 2T, second SU(2)-subgroup milestone

Track A of the SU(2) programme (Q8 → 2T → 2O → 2I). 2T, the 24 Hurwitz unit quaternions
{±1, ±i, ±j, ±k, (±1 ±i ±j ±k)/2}, with exact Eisenstein-integer (Z[ω]) characters: group and character
checks, electric spectrum, Gauss law, exact strong-coupling series of the Wilson loops on the open 2×3 strip
(two independent algorithms and bases) and the open 2×4 ladder (two bases), and exact integer-arithmetic
certificates of the 1×1 and 1×2 loops at finite coupling on the 2×3 strip. Every computation is a job of the
HEB-C HTTP-0 server (`hebc_ym_http0 serve`). Predictions were frozen in `PREREG.md` (`abef2496…`) before
any 2T code, including the W4 resolution rule learned from Q8. Rust 1.81, clippy `-D warnings`, integer
arithmetic only, json=0.

Hamiltonian: H/J_E = E + r Σ_p (2 − χ_f(U_p)), E = Σ_links Σ_{g∈S} (1 − L_g), S = the eight elements
(1 ± i ± j ± k)/2 (rotation angle 2π/3). Link energies: 0 (trivial), 4 (fundamental), 8 (3-dimensional),
10 (2′, 2″), 12 (1′, 1″). W = χ_f(U_C) = 2 Re of the loop's quaternion; W/2 is the (1/N) Re Tr normalisation.

## Result: G1–G4, R1, W1–W6 all PASS

| Check | Result |
| --- | --- |
| G1 group, characters | 24 elements, 7 classes, 7 irreps; axioms, row and column orthogonality in Z[ω], Σd² = 24 |
| G2 electric spectrum | characteristic polynomial of E_link on L²(2T) = x(x − 4)⁴(x − 8)⁹(x − 10)⁸(x − 12)² exactly |
| G3 Gauss law (2×3) | 76 orbit states × 6 vertices × 24 elements = 10944 local gauge transformations |
| G4 self-adjointness | every matrix exactly W-self-adjoint |
| R1 regression | sealed zwilson (S3) and zq8 (Q8) series reproduced exactly in both bases; S3 matrix byte-identical to `wilsonmat`, Q8 matrix shard identical to the sealed zq8 shard |
| W1 series | 2×3 (76-orbit elimination = 576-configuration spectral) and 2×4 (1312-orbit spectral = 13824-configuration spectral) identical to order 5; c_{A,A} = 1/8, 7/768, 11/16384 = (3/4)^A × S3, as derived; c_{A,A+1} = 0; E₁ = 2N_p, E₂ = −N_p/16, E₃ = E₅ = 0 |
| W2 certificates (2×3) | 8 values, widths 0.89–1.48 × 10⁻⁹, every lower end > 0 |
| W3 leading order | ⟨W⟩/(c_{A,A} r^A) at r = 1/27: [0.9999954, 0.9999957] (1×1), [0.99995, 1.00003] (1×2) |
| W4 convergence | both loops resolved (series deviation at r = 1/9: 4.0 × 10⁻⁵ and 9.3 × 10⁻⁵, above both interval widths) and both converge |
| W5 area law | ⟨W(1×2)⟩ < ⟨W(1×1)⟩ at every r; ratio at r = 1/27 = 0.0027005–0.0027007, predicted 0.0027006 |
| W6 even correction | D(1/9)/D(1/27) = [8.68, 9.35] for the 1×1 loop (r², as for Q8; S3 3.02, linear) |

### Exact series

2×3 strip (methods A and B identical):

| Loop | r | r² | r³ | r⁴ | r⁵ |
| --- | --- | --- | --- | --- | --- |
| 1×1 | **1/8** | 0 | −5/12288 | 0 | 2443/964689920 |
| 1×2 | 0 | **7/768** | 0 | −9/131072 | 0 |

2×4 ladder (methods A′ and B identical):

| Loop | r | r² | r³ | r⁴ | r⁵ |
| --- | --- | --- | --- | --- | --- |
| 1×1 | **1/8** | 0 | −5/12288 | 0 | 46159/17364418560 |
| 1×2 | 0 | **7/768** | 0 | −293/4718592 | 0 |
| 1×3 | 0 | 0 | **11/16384** | 0 | −90109/12079595520 |

Vacuum energy: E₀ = 4r − r²/8 + 5r⁴/24576 + O(r⁶) on 2×3 and 6r − 3r²/16 + 7r⁴/24576 + O(r⁶) on 2×4. The
leading coefficients and the first correction of the 1×1 loop (−5/12288) are the same on both lattices; the
geometry first shows at the next order. On the same 2×4 ladder the first correction is negative for 2T and
positive for Q8 (185/456192).

### Certified ⟨W(C)⟩ on the 2×3 strip (lower ends; exact rational intervals in `RESULT.hbp`)

| Loop | r = 1/27 | r = 1/9 | r = 1/3 | r = 1 |
| --- | --- | --- | --- | --- |
| 1×1 | 0.0046296082 | 0.0138883301 | 0.0416516060 | 0.1245956181 |
| 1×2 | 0.0000125022 | 0.0001125148 | 0.0010118841 | 0.0090465296 |

Each added plaquette multiplies the loop by 0.0027 at r = 1/27, 0.0081 at 1/9, 0.0243 at 1/3 and 0.0726 at 1.
At leading order the three groups line up exactly on the fundamental link energy E_f: c_{1,1} = 1/(2E_f)
(S3 1/6, Q8 1/12, 2T 1/8), and for every A, c_{A,A} scales as E_f^−A.

## How it is computed

- `fgroup`, 2T milestone: Z[ω] characters (exact pairs x + yω), 2T from exact Hamilton products of doubled
  Hurwitz coordinates, the coset map onto 2T/Q8 for the 1-dimensional characters, χ₃ = 4 Re(q)² − 1, and
  series method 2 (gauge-orbit basis with the spectral polynomial over all 58 possible nonzero electric
  eigenvalues, each solve checked exactly). Characteristic polynomials in big integers.
- Certified values by concavity of E₀ in a loop source λ = 10⁻⁸: (L₊ − U₀)/λ ≤ ⟨W⟩ ≤ (U₀ − L₋)/λ, each energy
  bracketed by `sturmrs` to one threshold unit (7.4 × 10⁻¹⁹ J_E) with J_E = 1.35 × 10¹⁵, the final brackets
  re-executed by `verify`; 836 counts.

## Replay and integrity chain

`run.sh` checks the prediction seal, builds the native/ym workspace and the driver with Rust 1.81, runs clippy,
starts `hebc_ym_http0 serve` on a fresh ROOT, runs everything through it and checks byte identity with
`RESULT.hbp`; `results-replay.yml` replays it in its own job on the exact head. `SHA256SUMS` covers every other
file, `ZT2.hbi` holds the chained `EMIT` rows (`pid[0] = sha16("EMITTER-SEED|SUN-ZT2-2026-10-30-CLOUD")`,
`pid[n] = sha16(pid[n-1]|name)`), and `ZT2.seal = sha256(SHA256SUMS + ZT2.hbi)`.

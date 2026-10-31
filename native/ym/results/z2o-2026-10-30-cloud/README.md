# z2o: binary octahedral group 2O, third SU(2)-subgroup milestone

Track A of the SU(2) programme (Q8 → 2T → 2O → 2I). 2O, the 48 unit quaternions: the 24 of 2T and the 24
elements (±u ±v)/√2 for distinct u, v in {1, i, j, k}, with exact characters in Z[√2]: group and character
checks, electric spectrum, Gauss law, exact strong-coupling series of the Wilson loops on the open 2×3 strip
(two independent algorithms and bases) over Q(√2), and exact integer-arithmetic certificates of the 1×1 and 1×2
loops at finite coupling on the 2×3 strip. Every computation is a job of the HEB-C HTTP-0 server
(`hebc_ym_http0 serve`). Predictions were frozen in `PREREG.md` (`a5c0b994…`) before any 2O code; addendum 1
(`1364b1a0…`, the counter for the certified values after the solver-scaling gate) and addendum 2 (`88cf2d17…`,
X1 and the square-root enclosure check) were sealed before any 2O count; addendum 3 (`c9da5a4a…`, the split CI
replay) was sealed before the CI replay it governs. Rust 1.81, clippy `-D warnings`,
exact integer and quadratic-field arithmetic only, json=0.

Hamiltonian: H/J_E = E + r Σ_p (2 − χ_f(U_p)), E = Σ_links Σ_{g∈S} (1 − L_g), S = the six elements
(1 ± i)/√2, (1 ± j)/√2, (1 ± k)/√2 (rotation angle π/2, one conjugacy class). Link energies: 0 (trivial),
6 − 3√2 (fundamental), 4 (3-dimensional), 6 (2 and 4), 8 (3′), 6 + 3√2 (fundamental′), 12 (sign).
W = χ_f(U_C) = 2 Re of the loop's quaternion; W/2 is the (1/N) Re Tr normalisation.

## Result: G1–G4, R1, S1, X1, W1–W6 all PASS

| Check | Result |
| --- | --- |
| G1 group, characters | 48 elements, 8 classes, 8 irreps; axioms, row and column orthogonality in Z[√2], Σd² = 48, one-dimensional characters multiplicative |
| G2 electric spectrum | characteristic polynomial of E_link on L²(2O) = x(x − 4)⁹(x − 6)²⁰(x − 8)⁹(x − 12)(x² − 12x + 18)⁴ exactly |
| G3 Gauss law (2×3) | 136 orbit states × 6 vertices × 48 elements = 39168 local gauge transformations |
| G4 self-adjointness | all 32 matrices exactly W-self-adjoint |
| R1 regression | sealed zwilson (S3), zq8 (Q8) and zt2 (2T) series reproduced exactly in both bases with quadratic-field arithmetic; S3 matrix byte-identical to `wilsonmat`, Q8 and 2T matrix shards identical to the sealed shards |
| S1 sturmq2 regression | exact and interval `sturmq2` reproduce the `sturmrs` counts at the final E₀(0) bracket thresholds of every certificate in zwilson, zq8 and zt2: 24 of 24, 0 at the lower end and 1 at the upper end |
| X1 exact against interval | on the 2O matrix at r = 1/27, exact `sturmq2` gives 0 and 1 at the two final E₀(0) thresholds, as the interval counts do |
| W1 series | 136-orbit elimination = 2304-configuration spectral (spectral degree 282), identical to order 5; c_{1,1} = (2 + √2)/12 and c_{2,2} = (21 + 14√2)/864, as derived; c_{A,A+1} = 0; E₁ = 4, E₂ = −(2 + √2)/12, E₃ = E₅ = 0 |
| W2 certificates | 8 values, widths 1.23–1.85 × 10⁻⁹, every lower end > 0 |
| W3 leading order | ⟨W⟩/(c_{A,A} r^A) at r = 1/27: [0.9999720, 0.9999722] (1×1), [0.999926, 0.999951] (1×2) |
| W4 convergence | both loops resolved (series deviation at r = 1/9: 2.5 × 10⁻⁴ and 5.6 × 10⁻⁴, above both interval widths) and both converge |
| W5 area law | ⟨W(1×2)⟩ < ⟨W(1×1)⟩ at every r; ratio at r = 1/27 = 0.0061467–0.0061469, predicted 0.0061470 |
| W6 even correction | D(1/9)/D(1/27) = [8.964, 9.023] for the 1×1 loop (r², as for Q8 and 2T) |

### Exact series

2×3 strip (methods A and B identical):

| Loop | r | r² | r³ | r⁴ | r⁵ |
| --- | --- | --- | --- | --- | --- |
| 1×1 | **(2 + √2)/12** | 0 | −3503/1233792 − 1285√2/616896 | 0 | (73057879 + 51422737√2)/789311029248 |
| 1×2 | 0 | **(21 + 14√2)/864** | 0 | −110059/103638528 − 312965√2/414554112 | 0 |

In decimals: c_{1,1} = 0.2845177969, c_{1,3} = −0.0057850342, c_{1,5} = 0.0001846935; c_{2,2} = 0.0472210531,
c_{2,4} = −0.0021296022.

Vacuum energy: E₀ = 4r − (2 + √2)r²/12 + (3503/2467584 + 1285√2/1233792) r⁴ + O(r⁶). The two plaquettes of
the strip are equivalent, so ⟨W(1×1)⟩ = 2 − ½ dE₀/dr, and the two series meet it exactly: c_{1,1} = −E₂ and
c_{1,3} = −2E₄.

Across the four groups the leading coefficients follow the fundamental link energy E_f: c_{1,1} = 1/(2E_f)
(S3 1/6, Q8 1/12, 2T 1/8, 2O (2 + √2)/12 with E_f = 6 − 3√2) and c_{2,2} = 7/(48E_f²) (S3 7/432, Q8 7/1728,
2T 7/768, 2O (21 + 14√2)/864). With the registered electric set 2O has the smallest fundamental link energy of
the four (E_f = 1.7573593129), and so the largest leading coefficients.

### Certified ⟨W(C)⟩ on the 2×3 strip (lower ends; exact rational intervals in `RESULT.hbp`)

| Loop | r = 1/27 | r = 1/9 | r = 1/3 | r = 1 |
| --- | --- | --- | --- | --- |
| 1×1 | 0.0105374012 | 0.0316051552 | 0.0946257617 | 0.2789119114 |
| 1×2 | 0.0000647703 | 0.0005826509 | 0.0052206287 | 0.0451875694 |

Ground-state energies E₀(0)/J_E (lower ends): 0.1477578685, 0.4409323198, 1.3017558706, 3.7183145481. Each
added plaquette multiplies the loop by 0.0061 at r = 1/27, 0.0184 at 1/9, 0.0552 at 1/3 and 0.1620 at 1.
Through order 5 the 1×1 series lies inside the certified interval at r = 1/27 and at r = 1/9; at r = 1/3 the
order-3 and order-5 partial sums fall just below and just above it.

## How it is computed

- `fgroup`, 2O lane: quadratic-field characters (exact pairs x + y√2), 2O from exact Hamilton products of the
  coordinates (x + y√2)/2 with closure checked, the QF exchange format for matrices whose diagonal carries
  Q(√2) values, and both series methods over Q(√2): method A gauge-orbit elimination (the reduced resolvent of
  the integer E is rational), method B configuration basis with the conjugation-closed candidate spectral
  polynomial, each solve checked exactly.
- `sturmq2`, inertia of B(T) = W (1000 (J_E E + J_B M + λ W_C) − T W) over Q(√2). Exact mode: LDLᵀ with
  symmetric pivoting and the exact sign test of a + b√2 (S1 and X1). Interval mode (addendum 1): the `sturmrs`
  algorithm with entries as outward-rounded integer intervals at scale 10⁸⁰ and √2 entering through the
  enclosure [L, L + 1]/10⁸⁰, L = isqrt(2 · 10¹⁶⁰), checked exactly before every count; every one of the 817
  interval counts of the certified values COMPLETE.
- Certified values by concavity of E₀ in a loop source λ = 4 × 10⁻⁹: (L₊ − U₀)/λ ≤ ⟨W⟩ ≤ (U₀ − L₋)/λ, each
  energy bracketed to one threshold unit (1.23 × 10⁻¹⁸ J_E) with J_E = 8.1 × 10¹⁴, the final brackets
  re-executed by `verify`; 867 counts in all (48 for S1, 2 for X1, 817 for the certified values).

## Runtime

- The run: 2026-10-31 02:10:12Z to 08:14:15Z, 21843 s on one 2.1 GHz cloud core; the driver writes its
  output when it exits.
- Exact counts: S1's eight Q8 counts 22.4–26.2 min each (176 states, entries up to 10062 bits), 3 h 16 min in
  all; S1's eight 2T counts 14–20 s (76 states, up to 4805 bits) and eight S3 counts about 2 s (49 states); the
  two X1 counts 80.5 and 81.8 min (136 states over Q(√2), entries up to 16889 bits). The addendum-1 timing count
  at E = 3 J_E reached 2599 bits in 136 s; at the final-bracket thresholds the entries reach 6.5 times as many
  bits.
- The 817 interval counts on 2O: 79 s in total, at most 175 ms each.
- Sampled peak resident memory (VmHWM): 54840 kB for an exact 2O count, 43668 kB for an exact Q8 count; at most
  753 MiB in use on the machine.
- A first start (2026-10-31 00:07Z, fresh root) was stopped by the execution environment's 2-hour limit on
  background tasks at the same point of S1 (after the Q8 counts at r = 1/9), before the driver wrote any
  output. The run above is the same driver and server, unchanged, on a fresh root.
- The replay (`run.sh` on a fresh root, the same machine, 2026-10-31 08:15:01Z to 14:40:15Z): 23114 s including
  the builds, every count again, output byte-identical to `RESULT.hbp`; S1's Q8 counts took 24.2–27.5 min each
  and the X1 counts 85.8 and 84.4 min.
- `RUNTIME.measure` records every exact count's threshold, result and elapsed time for both starts and the
  replay, the interval-count totals and the memory samples; it is runtime metadata, outside `RESULT.hbp`.
- On GitHub the full replay outgrows one hosted job: CI run 41211349149 cancelled it at the 360-minute limit
  (360 min 17 s) with every other package passing. The CI replay is therefore split (`ADDENDUM3.md`).

## Replay and integrity chain

`run.sh` checks the prediction and addendum seals, builds the native/ym workspace and the driver with Rust 1.81,
runs clippy, starts `hebc_ym_http0 serve` on a fresh ROOT, runs everything through it and checks byte identity
with `RESULT.hbp`: the full replay. On the exact head, `results-replay.yml` replays the package split
(`ADDENDUM3.md`). Six S1 part jobs run `run_s1.sh` (zwilson; zq8 bracket by bracket; zt2) and each recomputes its
MATRIX, COUNT and S1 lines through the server, byte-identical to the same lines of `RESULT.hbp`; the package job
runs `run.sh` with `Z2O_S1_FROM=RESULT.hbp`, which keeps those 84 lines in place and recomputes every other line,
the whole output byte-identical to `RESULT.hbp`. The gate requires all of them, so every line is recomputed by
exactly one job. The driver's split modes are environment-only; without them it runs the code of the official
run. `SHA256SUMS` covers every other file, `Z2O.hbi` holds the chained `EMIT` rows
(`pid[0] = sha16("EMITTER-SEED|SUN-Z2O-2026-10-30-CLOUD")`, `pid[n] = sha16(pid[n-1]|name)`), and
`Z2O.seal = sha256(SHA256SUMS + Z2O.hbi)`. The first seal (`fedcbf45…`, commit `85ca4f6`) covered the package
before the split; the package is resealed with addendum 3, `run_s1.sh` and the extended driver.

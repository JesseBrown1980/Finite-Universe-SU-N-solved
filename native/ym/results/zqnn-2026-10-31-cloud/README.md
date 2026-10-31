# zqnn: the qutrit variational network on the finite-group lane

The certified lane gives lower bounds on the vacuum energy E₀ from sturmrs inertia counts. This package
adds the upper side. A trained qutrit network writes down an explicit state whose energy U is an exact
rational, so U ≥ E₀. Each (group, coupling) cell on the open 2×3 strip therefore gets a two-sided enclosure
L ≤ E₀ ≤ U. On the 2×4 ladder the network gives U, compared with the sealed certificates of zwilson (S3) and
zq8 (Q8) and with the sealed 2T series.

The network is the one sealed as v7 in
[GGUF-files-operating-as-trained-quantum-neural-network-models](https://github.com/JesseBrown1980/GGUF-files-operating-as-trained-quantum-neural-network-models).
It has two layers of rotations with tan(angle/2) = b/6561, each b set by eight balanced trits. The second
layer is a correlated layer conjugated by a SUM ring. Every training run is an `fgqnn` job of
`hebc_ym_http0 serve`, and `verify` re-executes it byte for byte.

Predictions were frozen in `PREREG.md` (`d3e3b2b6…`) before any zqnn code. `ADDENDUM.md` (`9c2d5965…`) was
written before the driver: it covers the weak-coupling JE, the sealed ladder comparisons and five
development runs of the tool. Rust 1.81, clippy `-D warnings`, integer arithmetic only, json=0.

## Result: N1, N2, N3, N5, N6 PASS; N4 FAIL

| Check | Result |
| --- | --- |
| N1 exactness | every trained network has Σψ² = (∏(a² + b²))², and its inverse returns that square times the input |
| N2 verify | all 81 verify re-executions pass, the 45 network jobs among them |
| N3 enclosure (2×3) | L ≤ U in all 18 cells, and sturmrs finds an eigenvalue at or below U in every cell |
| N4 strong coupling | **FAIL** in 6 of the 21 cells with r ≤ 1: Q8 at r = 1/27 on both geometries, 2T at r = 1/27 and 1/9 on 2×3 and at r = 1/27 and 1/3 on 2×4; U ends just above the uniform-state energy 2rN_p |
| N5 weak coupling | at r = 3 and r = 9 the IDENTITY seed always ends below its start |
| N6 half the gap | at r = 1/27 and r = 9, for every group, U − L ≤ (S − L)/2 |

### Enclosures on the 2×3 strip (L certified, U from the network)

| Group | r | L | U | U − L | gap closed |
| --- | --- | --- | --- | --- | --- |
| S3 | 1/27 | 0.1479188174 | 0.1479254137 | 6.6 × 10⁻⁶ | 99.4% |
| S3 | 1/9 | 0.4423676644 | 0.4424328761 | 6.5 × 10⁻⁵ | 96.8% |
| S3 | 1/3 | 1.3142912178 | 1.3147280940 | 4.4 × 10⁻⁴ | 97.3% |
| S3 | 1 | 3.8187374278 | 3.8341571262 | 0.0154 | 91.0% |
| S3 | 3 | 10.0877381961 | 10.4586023167 | 0.371 | 80.2% |
| S3 | 9 | 17.4737214571 | 17.7966436395 | 0.323 | 90.8% |
| Q8 | 1/27 | 0.1480338359 | 0.1483886563 | 3.5 × 10⁻⁴ | 86.0% |
| Q8 | 1/9 | 0.4434156070 | 0.4441918972 | 7.8 × 10⁻⁴ | 74.1% |
| Q8 | 1/3 | 1.3240715708 | 1.3287247007 | 0.0047 | 53.0% |
| Q8 | 1 | 3.9164639347 | 3.9558440469 | 0.039 | 50.9% |
| Q8 | 3 | 11.2336098710 | 11.5921294016 | 0.359 | 52.3% |
| Q8 | 9 | 28.0276103283 | 31.1140629300 | 3.09 | 61.0% |
| 2T | 1/27 | 0.1479766808 | 0.1496187594 | 1.6 × 10⁻³ | 51.8% |
| 2T | 1/9 | 0.4429012656 | 0.4457153242 | 2.8 × 10⁻³ | 42.9% |
| 2T | 1/3 | 1.3194469550 | 1.3294556320 | 0.010 | 43.5% |
| 2T | 1 | 3.8752026095 | 3.9518876215 | 0.077 | 41.0% |
| 2T | 3 | 10.8908845405 | 11.4189112331 | 0.528 | 52.8% |
| 2T | 9 | 26.8734138410 | 29.9090259144 | 3.04 | 66.8% |

"Gap closed" is (S − U)/(S − L), with S the start energy of the seed that reached U. The certified L values
are exact rational lower ends, one threshold unit wide, in `RESULT.hbp`.

### The 2×4 ladder

| Group | r | U | sealed comparison | U − comparison |
| --- | --- | --- | --- | --- |
| S3 | 1/27, 1/3, 1 | 0.22211, 1.99025, 5.88732 | zwilson certified L | 2.3 × 10⁻⁴, 0.019, 0.159 |
| Q8 | 1/27, 1/3, 1 | 0.22287, 1.99576, 5.95610 | zq8 certified L | 8.2 × 10⁻⁴, 0.0097, 0.081 |
| 2T | 1/27, 1/3, 1 | 0.22566, 2.00031, 5.96114 | sealed series to r⁵ | 0.0037, 0.021, 0.148 |

Every ladder U lies above the sealed certified lower end, as it must.

### Reading

- **S3 is close to tight.** Its elements fit three-level registers, and the enclosure is 6.6 × 10⁻⁶ wide at
  r = 1/27.
- **Q8 and 2T sit in mismatched registers.** They have 8 and 24 elements, held in 9 and 27 qutrit codes, and
  the codes outside the group are discarded. A rational rotation cannot make the valid codes exactly
  uniform, and at the smallest couplings that error is larger than the strong-coupling gain. This is the
  N4 failure.
- **The next design.** A network for the SU(2) subgroups should act on group registers directly: rotations
  between group elements, and multiplications by electric-set elements in place of SUMs.

## Replay and integrity chain

`run.sh` reproduced `RESULT.hbp` byte for byte on a fresh build and a fresh server ROOT (about 7 minutes).
It checks the prediction and addendum seals, builds the native/ym workspace and the driver with Rust
1.81, runs clippy, starts `hebc_ym_http0 serve` on a fresh ROOT, runs everything through it and checks byte
identity with `RESULT.hbp`. `SHA256SUMS` covers every other file. `ZQNN.hbi` holds the chained `EMIT` rows
(`pid[0] = sha16("EMITTER-SEED|SUN-ZQNN-2026-10-31-CLOUD")`, `pid[n] = sha16(pid[n-1]|name)`), and
`ZQNN.seal = sha256(SHA256SUMS + ZQNN.hbi)`.

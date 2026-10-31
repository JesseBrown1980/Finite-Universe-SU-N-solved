# Finite-Universe SU(N)

Exact integer-arithmetic certificates of finite-lattice gauge theories for finite groups on the road to
SU(N), built under the HEB-C rules: Rust 1.81, clippy `-D warnings`, integer and exact rational arithmetic
only, json=0, every computation a job of the HEB-C HTTP-0 server (`hebc_ym_http0 serve`), and the
HBI → HBP → SHA → SH → HASH integrity chain on every result package.

Sequence:
1. S3 confinement: static-charge potential, string tension, Wilson loops, volume scaling.
2. SU(2), two tracks. Track A, finite subgroups: Q8, binary tetrahedral (2T), binary octahedral (2O),
   binary icosahedral (2I), on the same finite-group machinery, with exact arithmetic in Z, Z, Z[√2] and
   Z[φ] for their fundamental characters. Track B, continuous SU(2) in the representation basis with an
   explicit cutoff j_max and rigorous truncation-error enclosures. The subgroups are compared with each
   other and with track B.
3. 3+1-dimensional lattice gauge Hamiltonians.
4. SU(2) continuum scaling with controlled group discretisations.

## Results

| Package | What it measures | Result |
| --- | --- | --- |
| [`results/zconf-2026-10-30-cloud`](results/zconf-2026-10-30-cloud) | S3 static-charge potential and string tension, strong-coupling side, to O(r²) | V(d,r) = 3d − (5/72) d r² + O(r³); F1–F4 PASS; two independent methods give identical rationals |
| [`native/ym/results/zconfcert-2026-10-30-cloud`](native/ym/results/zconfcert-2026-10-30-cloud) | S3 static-charge energies certified at finite coupling (r = 1/27 … 1) through HTTP-0 | Q1–Q7 PASS; κ_eff within 1% of κ_g at r = 1/27; ladder/strip ratio 1.9999 (one against two plaquettes per string link); string tension certified positive in all 8 cells |
| [`native/ym/results/zwilson-2026-10-30-cloud`](native/ym/results/zwilson-2026-10-30-cloud) | S3 spatial Wilson loops 1×1, 1×2, 1×3: exact series and certified finite-coupling values | W1–W5 PASS; ⟨W(1×A)⟩ starts at r^A with 1/6, 7/432, 11/6912 (two bases identical); area law certified at every coupling |
| [`native/ym/results/zq8-2026-10-30-cloud`](native/ym/results/zq8-2026-10-30-cloud) | Q8, first SU(2) subgroup: group and Gauss-law checks, exact Wilson-loop series, certified finite-coupling loops | G1–G4, R1, W1–W3, W5, W6 PASS; W4 FAIL (1×3 convergence below resolution, kept as registered). Leading coefficients 1/12, 7/1728, 11/55296 = 2^−A × S3; centre symmetry removes every odd correction (c_{A,A+1} = 0, r² finite-coupling correction against S3's linear one) |
| [`native/ym/results/zt2-2026-10-30-cloud`](native/ym/results/zt2-2026-10-30-cloud) | 2T (binary tetrahedral, 24 elements, Z[ω] characters): group, spectrum and Gauss-law checks, exact series on 2×3 and 2×4, certified loops on 2×3 | G1–G4, R1, W1–W6 PASS; c_{A,A} = 1/8, 7/768, 11/16384 = (3/4)^A × S3; odd corrections vanish; W4 resolved and converging for both loops; across S3, Q8, 2T the leading coefficients follow c_{1,1} = 1/(2E_f) with E_f the fundamental link energy |
| [`native/ym/results/z2o-2026-10-30-cloud`](native/ym/results/z2o-2026-10-30-cloud) | 2O (binary octahedral, 48 elements, Z[√2] characters): group, spectrum and Gauss-law checks, exact series over Q(√2) on 2×3, certified loops on 2×3, the exact Q(√2) counter checked against every sealed bracket | G1–G4, R1, S1, X1, W1–W6 PASS; c_{1,1} = (2 + √2)/12 and c_{2,2} = (21 + 14√2)/864, as predicted from c_{1,1} = 1/(2E_f) and c_{2,2} = 7/(48E_f²) with E_f = 6 − 3√2, the law the four groups now share; odd corrections vanish; W4 resolved and converging for both loops |

The laws the packages share, their predictions for the 2O 2×4 ladder and 2I, and the exact counts that size
those packages: [`docs/FINITE_GROUP_LAWS.md`](docs/FINITE_GROUP_LAWS.md).

The `native/` directory is the HEB-C kernel and HTTP-0 toolchain, copied from the Yang–Mills repository
and extended here with the finite-group lane (`native/ym/crates/fgroup`); see
[`native/TOOLCHAIN.md`](native/TOOLCHAIN.md). `native/TOOLCHAIN.SHA256` lists every file and is checked on
every pull request. The two `native/ym/results` packages are the sealed originals from the
Yang–Mills repository, at their original paths; here they rebuild this copy of the toolchain and run every
job through its own HTTP-0 server.

Every package carries its frozen `PREREG.md`, `RESULT.hbp`, `run.sh`, `SHA256SUMS`, the chained `.hbi`
and the `.seal`. `.github/workflows/results-replay.yml` verifies the chain and replays every package
byte for byte on the exact head of each pull request.

Companion repository: [Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C).

# Finite-Universe SU(N)

Exact integer-arithmetic gauge-theory certificates for finite groups on the road to SU(N), built
under the HEB-C rules: Rust 1.81 (std only), clippy `-D warnings`, integer and exact rational
arithmetic only, json=0, and the HBI → HBP → SHA → SH → HASH integrity chain on every result package.

Sequence:
1. S3 confinement: static-charge potential, string tension, Wilson loops, volume scaling.
2. SU(2) finite subgroups: Q8, binary tetrahedral (2T), binary octahedral (2O), binary icosahedral (2I).
3. 3+1-dimensional lattice gauge Hamiltonians.
4. SU(2) continuum scaling with controlled group discretisations.

## Results

| Package | What it measures | Result |
| --- | --- | --- |
| [`results/zconf-2026-10-30-cloud`](results/zconf-2026-10-30-cloud) | S3 static-charge potential and string tension, strong-coupling side, to O(r²) | V(d,r) = 3d − (5/72) d r² + O(r³); F1–F4 PASS; two independent methods give identical rationals |

Every package carries its frozen `PREREG.md`, `RESULT.hbp`, `run.sh`, `SHA256SUMS`, the chained `.hbi`
and the `.seal`. `.github/workflows/results-replay.yml` verifies the chain and replays every package
byte for byte on the exact head of each pull request.

Companion repository: [Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C).

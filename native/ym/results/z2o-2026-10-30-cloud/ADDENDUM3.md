# z2o addendum 3: split CI replay (sealed before the CI run it governs)

The full replay does not fit one GitHub-hosted job. CI run 41211349149 on 85ca4f6 cancelled the z2o job at its
360-minute limit (360 min 17 s); the toolchain manifest and the other five package jobs passed. The run itself
took 21843 s and the local full replay 23114 s on one 2.1 GHz core, byte-identical. Under the reviewing seat's
decision (the z2o job at 360 minutes; if CI exceeds 6 h, S1 is split into mandatory replay shards, never
omitted), the CI replay is split three ways:

- S1 part jobs, in parallel: zwilson (part 0, all four brackets), zq8 bracket by bracket (parts 1:0, 1:1, 1:2,
  1:3) and zt2 (part 2). Each runs `run_s1.sh <part>`: it checks the seals, builds and lints as run.sh does,
  recomputes the part's MATRIX, COUNT and S1 lines through `hebc_ym_http0 serve` on a fresh ROOT
  (Z2O_PART=s1:<part>), reads the same lines from RESULT.hbp (Z2O_PART=s1rec:<part>, each line found exactly
  once), and requires the two byte-identical.
- The main job runs `run.sh` with Z2O_S1_FROM=RESULT.hbp. The driver puts the 84 recorded S1 lines in their
  place and recomputes every other line through the server; the whole output must equal RESULT.hbp byte for
  byte, which also fixes the place and order of the S1 lines.
- The gate passes only when the manifest, every package job and every S1 part job pass. Every line of RESULT.hbp
  is then recomputed by exactly one job.

The driver change is additive: with neither variable set, z2o/src/main.rs runs exactly the code of the official
run and of the local full replay (sha256 ceacdd006eb9ec12df5e7dd0c733986e7fadeb90b7f0ea5736e53e03bb0a51fd), and
run.sh without Z2O_S1_FROM remains the full replay. RESULT.hbp
(5aee1bccfc945084367bdbe2c28d655b4d418e2b51f1c5d03ad57939e81bde29) and every registered check, threshold and
prediction are unchanged. The package seal of 85ca4f6 (fedcbf45ff5afa7941b1c1e7078cc9ed3f0aca8ec5affa8057a1fa6d65698b08)
stays in the history; the package is resealed with this addendum, run_s1.sh and the extended driver.

PREREG.md (a5c0b994…), addendum 1 (1364b1a0…) and addendum 2 (88cf2d17…) are otherwise unchanged.

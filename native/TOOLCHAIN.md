# HEB-C toolchain copy

`native/` holds a copy of the HEB-C kernel and the `native/ym` workspace (HTTP-0 backend
`hebc_ym_http0`, `torusl` with `ymconf`, `sturmrs` with `ymwser`, `evolvers`, and the vendored
`num-bigint`, `num-integer`, `num-traits`, `autocfg` crates) from
[Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C)
at commit `90dcf73b3e459c0be24e3db07208c4182cbd5d10` (main after PR #58), copied byte for byte.
`TOOLCHAIN.SHA256` lists every file of the copy; the replay workflow checks it on every pull request, so any
later change to the toolchain in this repository updates the manifest in the same pull request and shows
up in its diff.

## Changes in this repository since the copy

- `ym/backend/src/main.rs`: `tool()` appends the platform executable suffix (the same one-line fix as
  Yang–Mills PR #59, so the backend finds `fgconf.exe` and the other tools on Windows); new job kinds
  `fgcheck`, `fgwilsonmat`, `fgwilsoncount` and `fgwilsonseries` for the finite-group lane.
- `ym/crates/fgroup` (new): finite-group data (S3 for regression, Q8 from Gaussian-integer matrices), the
  open-strip builder `fgconf`, the exact series solver `fgwser` and the acceptance checks `fgcheck`.
- `ym/Cargo.toml`, `ym/Cargo.lock`: `crates/fgroup` added to the workspace.

The S3 packages `native/ym/results/zconfcert-2026-10-30-cloud` (certified static-charge energies) and
`native/ym/results/zwilson-2026-10-30-cloud` (certified spatial Wilson loops) are copied unchanged, at their
original paths, so their `run.sh` files and seals are the sealed originals; here they rebuild this copy of the
toolchain, run every job through its own `hebc_ym_http0 serve`, and reproduce their `RESULT.hbp` byte for byte.

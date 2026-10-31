#!/bin/bash
# Split replay of z2o (ADDENDUM3.md), one S1 part: <si> or <si>:<bracket> (sealed package 0 zwilson, 1 zq8,
# 2 zt2; bracket 0..3). Checks the prediction and addendum seals, builds the native/ym workspace and the driver
# with Rust 1.81, runs clippy -D warnings on both, recomputes the part's MATRIX, COUNT and S1 lines through
# `hebc_ym_http0 serve` on a fresh ROOT, and checks them byte for byte against the same lines of RESULT.hbp.
# Never overwrites a record. json=0.
set -u
cd "$(dirname "$0")"
P=${1:-}
case "$P" in 0|2|1:0|1:1|1:2|1:3) ;; *) echo "REPLAY_S1|test=z2o|part=$P|result=INCONCLUSIVE|reason=part|json=0"; exit 2 ;; esac
sha256sum -c --quiet PREREG.md.sha256 ADDENDUM.md.sha256 ADDENDUM2.md.sha256 ADDENDUM3.md.sha256 || { echo "REPLAY_S1|test=z2o|part=$P|result=INCONCLUSIVE|reason=seal|json=0"; exit 2; }
W=$(cd ../.. && pwd)
(cd "$W" && cargo +1.81.0 build -q --release --offline --locked) || { echo "REPLAY_S1|test=z2o|part=$P|result=INCONCLUSIVE|reason=workspace_build|json=0"; exit 3; }
(cd "$W" && cargo +1.81.0 clippy -q --release --offline --locked --workspace -- -D warnings) || { echo "REPLAY_S1|test=z2o|part=$P|result=FAIL|reason=workspace_clippy|json=0"; exit 4; }
export CARGO_TARGET_DIR=$(mktemp -d)
cargo +1.81.0 build -q --release --locked --manifest-path z2o/Cargo.toml || { echo "REPLAY_S1|test=z2o|part=$P|result=INCONCLUSIVE|reason=build_failed|json=0"; exit 3; }
cargo +1.81.0 clippy -q --release --locked --manifest-path z2o/Cargo.toml -- -D warnings || { echo "REPLAY_S1|test=z2o|part=$P|result=FAIL|reason=clippy|json=0"; exit 4; }
ROOT=$(mktemp -d); PORT=18777; LOG=$(mktemp)
"$W/target/release/hebc_ym_http0" serve "$ROOT" "$PORT" 10000 > "$LOG" 2>&1 &
SRV=$!
for _ in $(seq 1 100); do grep -q 'server=READY' "$LOG" && break; sleep 0.1; done
OUT=$(mktemp); REC=$(mktemp)
Z2O_PART="s1:$P" "$CARGO_TARGET_DIR/release/z2o" "$PORT" > "$OUT"; live=$?
kill "$SRV" 2>/dev/null; wait "$SRV" 2>/dev/null
Z2O_PART="s1rec:$P" Z2O_S1_FROM=RESULT.hbp "$CARGO_TARGET_DIR/release/z2o" 0 > "$REC"; rec=$?
same=0; [ $live = 0 ] && [ $rec = 0 ] && [ -s "$OUT" ] && cmp -s "$OUT" "$REC" && same=1
echo "REPLAY_S1|test=z2o|part=$P|lines=$(wc -l < "$OUT")|live_exit=$live|recorded_exit=$rec|rustc=$(rustc +1.81.0 --version | tr ' ' '_')|src_sha256=$(sha256sum z2o/src/main.rs | cut -c1-64)|backend_sha256=$(sha256sum "$W/backend/src/main.rs" | cut -c1-64)|sturmq2_sha256=$(sha256sum "$W/crates/fgroup/src/bin/sturmq2.rs" | cut -c1-64)|byte_identical=$same|json=0"
[ $same = 1 ] || diff "$OUT" "$REC" | head -20
rm -rf "$CARGO_TARGET_DIR" "$ROOT" "$OUT" "$REC" "$LOG"
[ $same = 1 ]

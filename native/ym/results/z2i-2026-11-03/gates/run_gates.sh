#!/bin/bash
# z2i gates, timing only, at a threshold outside every registered window (PREREG.md, "Gates"). Below normal priority, detached, read-only on the store of nothing else.
# Part A: interval and exact counts on the 2I 2 x 3 matrix at r = 1/27, lambda = 0, at E = 3 J_E. Part B: method B at orders 1 and 2.
export MSYS_NO_PATHCONV=1
T=C:/tmp/su-n-z2i-target/release
O=/c/tmp/su-n-z2i-work/gates
cd "$O" || exit 9
now() { date -u +%Y-%m-%dT%H:%M:%SZ; }
part=$1
if [ "$part" = A ]; then
  echo "GATE_A_START|at=$(now)|json=0" >> GATES.hbp
  "$T/fgconf.exe" 5 2 3 0 1 2 m-r27.txt 30000000000000 0 > fgconf-r27.out 2>&1
  echo "GATE_MATRIX|at=$(now)|exit=$?|sha256=$(sha256sum m-r27.txt | cut -c1-64)|line=$(cut -c1-200 fgconf-r27.out | head -1)|json=0" >> GATES.hbp
  "$T/sturmq2.exe" m-r27.txt 810000000000000 30000000000000 2430000000000000000 interval > count-interval-3JE.out 2>&1
  echo "GATE_COUNT|at=$(now)|exit=$?|mode=interval|threshold=3JE|line=$(cut -c1-420 count-interval-3JE.out | head -1)|json=0" >> GATES.hbp
  "$T/sturmq2.exe" m-r27.txt 810000000000000 30000000000000 2430000000000000000 exact > count-exact-3JE.out 2>&1
  echo "GATE_COUNT|at=$(now)|exit=$?|mode=exact|threshold=3JE|line=$(cut -c1-420 count-exact-3JE.out | head -1)|json=0" >> GATES.hbp
  echo "GATE_A_DONE|at=$(now)|json=0" >> GATES.hbp
else
  echo "GATE_B_START|at=$(now)|json=0" >> GATES.hbp
  for order in 1 2; do
    t0=$(date +%s)
    "$T/fgwser.exe" 5 2 3 1 $order > method-b-order$order.out 2> method-b-order$order.err
    c=$?
    echo "GATE_METHOD_B|at=$(now)|order=$order|exit=$c|secs=$(( $(date +%s) - t0 ))|basis=$(grep FGSERIES_BASIS method-b-order$order.out | cut -c1-260)|json=0" >> GATES.hbp
  done
  echo "GATE_B_DONE|at=$(now)|json=0" >> GATES.hbp
fi

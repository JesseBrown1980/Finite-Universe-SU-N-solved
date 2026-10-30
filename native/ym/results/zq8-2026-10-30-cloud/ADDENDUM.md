# zq8 addendum 1 (sealed before any Q8 code is written or any count is run)

PREREG.md fixes JE = 1350000000000000. At r = 1 the registered bracket for E0(0) is the second-order
prediction 1000 JE (6 - 1/8) = 7.93 * 10^18 plus or minus m = 10 x |second-order term| = 1.69 * 10^18 T units,
so its upper end, 9.62 * 10^18, exceeds the backend's signed 64-bit spec field (9.22 * 10^18) and every r = 1
cell would be INCONCLUSIVE by construction.

Change: JE = 1080000000000000 (27 x 4 * 10^13), so JB = r JE stays an integer, and LAM = 10800000 (lambda =
10^-8 unchanged). At r = 1 the bracket is then 6.35 * 10^18 +- 1.35 * 10^18, upper end 7.70 * 10^18, and one
doubling still fits (9.05 * 10^18). One T unit becomes 9.3 * 10^-19 J_E; the expected <W> widths (about
1.2 * 10^-9 for 1x1) remain inside the registered W2 bound of 2 * 10^-9. Every prediction, threshold,
geometry and method is otherwise unchanged.

# z2o addendum 2 (sealed before any 2O count, series job or official run)

Two safeguards from the reviewing seat, added before the run they govern:

- X1, exact against interval on 2O itself: at r = 1/27, on the 2O 2 x 3 matrix with lambda = 0, exact sturmq2 must
  give the same counts as the interval mode at the two final E0(0) bracket thresholds (0 at the lower end,
  >= 1 at the upper end). A disagreement is a FAIL of X1; an exact or interval count that does not complete is
  INCONCLUSIVE.
- The interval mode checks its square-root enclosure exactly before counting (L^2 <= 2 * 10^160 < (L + 1)^2,
  L = isqrt(2 * 10^160)) and records it in its result line; a failed check refuses the count. Interval results
  record mode=interval; a count is accepted only with status=COMPLETE, anything else is INCONCLUSIVE.

PREREG.md (a5c0b994...) and addendum 1 (1364b1a0...) are otherwise unchanged.

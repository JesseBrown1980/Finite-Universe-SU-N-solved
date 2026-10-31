# z2o addendum 1 (sealed before any 2O count, series job or official run)

## Change: the counter for the certified values
PREREG.md registers sturmq2 as exact LDL^T over Q(sqrt2). The solver-scaling gate (proposed by the reviewing
seats) was run before any 2O count, at thresholds outside every registered window:
- sealed 2T matrix (76 orbit states, r = 1/27, JE = 1.35 * 10^15): exact sturmq2 1.3 s, 1.4 s and 19.5 s per
  count (largest entry 745, 773 and 4895 bits) against sturmrs 22-24 ms, with identical counts (0, 1, 1) at
  T = 1 and at E = 3 and 6 J_E;
- 2O 2 x 3 matrix (136 orbit states, r = 1/27, lambda = 0, JE = 8.1 * 10^14): exact sturmq2 136 s for one count
  (largest entry 2599 bits; negative = 1 at E = 3 J_E).
The certified values need about 800 counts, about 30 hours exact, beyond any exact-head replay.

The certified-value counts (W2-W6) therefore use sturmq2 in interval mode: the sturmrs algorithm unchanged
(entries as outward-rounded intervals of integers scaled by 10^80, Schur updates enclosed by floor and ceiling,
a pivot taken only when its interval excludes zero, the count COMPLETE only when every pivot is decided,
otherwise UNCERTAIN, which makes the cell INCONCLUSIVE), with each theta part entering through the exact
enclosure sqrt2 in [L, L + 1] / 10^80, L = isqrt(2 * 10^160). For an integer matrix it is the same computation as
sturmrs. Exact mode is kept for S1.

## S1, sharpened
S1 is run in both modes: exact sturmq2 and interval sturmq2 must each reproduce the sturmrs counts at every
final E0(0) bracket threshold of zwilson, zq8 and zt2 (0 at the lower end, >= 1 at the upper end).

Every other registered item (group, characters, G1-G4, R1, series methods, couplings, lambda, JE, bracketing,
widths, W1-W6) is unchanged.

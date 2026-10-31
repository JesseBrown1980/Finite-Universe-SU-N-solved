# zqnn addendum, written after the fgqnn tool was built and before the driver was written or run

1. Weak-coupling threshold field. The sturmrs threshold is a 64-bit spec field (at most 9.2 x 10^18). At
   JE = 1350000000000000 a bisection up to an energy near 35 J_E (r = 9) needs thresholds near 4.7 x 10^19,
   outside the field. The cells at r = 3 and r = 9 therefore use JE = 135000000000000 (one threshold unit
   = 7.4 x 10^-18 J_E); the cells at r <= 1 keep JE = 1350000000000000 as registered. No check changes.
2. Sealed comparisons on the 2 x 4 ladder. Besides the sealed series, the sealed certified vacuum energies
   of zwilson (S3, JE = 270000000000000) and zq8 (Q8, JE = 1080000000000000), the ENERGY0 rows of their
   RESULT.hbp (sha256 pinned), give certified lower ends at r = 1/27, 1/3, 1. They are reported as
   measured U - L on the ladder, not as a new check.
3. Development runs. While the tool was built, fgqnn ran directly (not through HTTP-0) five times to test
   it: GROUP 1, 2, 3 on 2 x 3 at r = 1/27, seed UNIFORM, SWEEPS 1; GROUP 3 on 2 x 4 at r = 1/3, seed
   UNIFORM, SWEEPS 1; GROUP 2 on 2 x 3 at r = 9, seed IDENTITY, SWEEPS 2. None of these is a registered
   cell configuration (every registered S2 cell runs SWEEPS 4, every S3 cell SWEEPS 2 at r = 1/27, 1/3, 1);
   no design or prediction changed after them.

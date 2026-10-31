# zqnn, cloud seat: the qutrit variational network on the finite SU(2)-subgroup lane (Finite-Universe-SU-N-solved).
Pre-registration, frozen before any zqnn code is written or run.

The certified lane gives lower bounds: sturmrs inertia counts bracket E0 from below on matrices it can
factor in replay time. This package adds the other side. An explicit variational state, built by a trained
qutrit network, gives an exact rational upper bound U >= E0 on the same Hamiltonian, so each (group,
coupling) cell gets a two-sided enclosure L <= E0 <= U. On the 2 x 4 ladder, where sturmrs is too slow in
replay time, the network gives U alone. The network is the one sealed as v7 in
JesseBrown1980/GGUF-files-operating-as-trained-quantum-neural-network-models (rotations with
tan(angle/2) = b/6561, b from eight balanced trits; a correlated layer conjugated by a SUM ring), carried
over to group registers. Rust 1.81, clippy -D warnings, integer arithmetic only, json=0, every computation a
job of the HTTP-0 backend (hebc_ym_http0 serve), HBI -> HBP -> SHA -> SH -> HASH with a separate seal.

## Groups, Hamiltonian, geometries
- Groups from the sealed fgroup lane: S3 (GROUP 1), Q8 (2), 2T (3); electric sets, fundamental characters
  and H/J_E = E + r M with M = sum_p (2 - chi_f(U_p)) exactly as in zq8 and zt2.
- Geometries: S2 = open 2 x 3 strip (2 plaquettes), S3 = open 2 x 4 ladder (3 plaquettes); comb tree gauge
  and the residual conjugation orbits of fgroup::basis.
- 2O is not in this package: its fundamental character takes the values +-sqrt 2, so M leaves the integers
  and the sturmrs exchange format. It is a separate package.

## The network (tool fgqnn in crates/fgroup, run by the new job kinds below)
- Register: each free link is k qutrits, k the least with 3^k >= |G| (S3 and Q8: k = 2; 2T: k = 3); code
  value v = sum_j d_j 3^j; codes v >= |G| are not group elements and are discarded.
- Gate R_ij(b) on levels i, j of one qutrit: [[a^2 - b^2, -2ab], [2ab, a^2 - b^2]], a^2 + b^2 on the third
  level, a = 6561, b = 2187 p0 + 729 p1 + ... + p7, p in {-1, 0, 1}; b = 0 is skipped.
- Layer 1: R_01(b), R_02(b') on every qutrit. Layer 2: SUM ring over all qutrits (q -> q + 1), R_01, R_02 on
  every qutrit, inverse ring; with its trits zero layer 2 is the identity.
- Amplitudes are big integers. The state psi on configurations is projected onto conjugation-invariant
  states, x_i = (|G| / w_i) sum over the orbit, and the energy is the exact rational
  U = (sum_i w_i x_i (E x)_i + r sum_i w_i M_i x_i^2) / sum_i w_i x_i^2 with E, M, w from fgroup::basis.
- Training: coordinate descent over the trits, each trit tries the other two values and keeps strict
  decreases, at most SWEEPS sweeps. Seeds: UNIFORM (layer 1 b = 2065, b' = 2718 on every qutrit, the
  sealed v6/v7 near-uniform qutrit; layer 2 zero) and IDENTITY (all trits zero: every link the identity
  element).

## Jobs
- fgqnn/GROUP/ROWS/COLS/RNUM/RDEN/SEED/SWEEPS: train and print the trained trits, the energy numerator and
  denominator, the exact norm and inverse checks; verify re-executes it byte for byte.
- Lower bounds: the sealed fgwilsonmat/fgwilsoncount jobs in vacuum (LAM = 0, loop 0..1), JE =
  1350000000000000, JB = r JE, bisected to one threshold unit 1/(1000 JE); final brackets re-executed by
  verify.

## Cells
- S2 (2 x 3) for S3, Q8, 2T at r = 1/27, 1/9, 1/3, 1, 3, 9: both seeds, SWEEPS = 4; certified L and the
  network U.
- S3 (2 x 4) for S3, Q8, 2T at r = 1/27, 1/3, 1: seed UNIFORM, SWEEPS = 2; U only, compared with the sealed
  vacuum series of zq8 and zt2 where they exist.

## Predictions (each an exact PASS)
- N1 exactness: every trained network has sum psi^2 equal to the square of the product of a^2 + b^2 over its
  applied gates, and the inverse network returns that square times the input.
- N2 every fgqnn job re-executes byte for byte under verify.
- N3 enclosure on S2: in every cell the certified lower end L is <= U, and the sturmrs count at the
  threshold ceil(1000 JE U) is at least 1 (an eigenvalue at or below U).
- N4 strong coupling: at every r <= 1 on both geometries, U is below 2 r N_p, the energy of the uniform
  state (E = 0, every plaquette at M = 2 on average).
- N5 weak coupling on S2: at r = 3 and r = 9 the IDENTITY seed ends below its start, the energy of the
  all-identity configuration.
- N6 on S2 the best network closes at least half of the gap from its seed: U - L <= (S - L) / 2, with S the
  start energy of the seed that reached U, at r = 1/27 and at r = 9, for every group.
Measured, not predicted: U - L in every S2 cell; the fraction of the second-order drop the network
captures at small r; U against the sealed series on S3.
Fail-closed: any failed job, failed verify or uncertified bracket makes the affected check INCONCLUSIVE.

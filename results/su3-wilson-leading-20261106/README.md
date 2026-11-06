# Continuous SU(3): four-dimensional leading Wilson-loop coefficient

**Status:** exploratory exact strong-coupling benchmark, **awaiting independent Opus review**. Prediction frozen on GitHub before the engine in commit [`dc1d29a837c079c141325b5d02593bfefc8635d5`](https://github.com/JesseBrown1980/Finite-Universe-SU-N-solved/commit/dc1d29a837c079c141325b5d02593bfefc8635d5). This is *not* a finite-coupling lattice simulation, glueball mass measurement, or proof of four-dimensional continuum QCD.

## Model / observable

Pure continuous SU(3) Wilson gauge action on an **open** 4D hypercubic lattice (4 vertices per axis), without quarks:

\[
\exp(-S)=\prod_p \exp\left[\frac{\beta}{6}(\mathrm{Tr}U_p+\mathrm{Tr}U_p^\dagger)\right],\quad W(R,T)=\frac13\mathrm{Re}\,\mathrm{Tr}U_C.
\]

This fixes the standard convention \(\beta=6/g^2\). `PREDICTION.hbp` records the nine leading-order tests for rectangular \(R\times T\) contours. The known Wilson strong-coupling result is

\[
\langle W(R,T)\rangle=(\beta/18)^{RT}+\text{higher orders in }\beta.
\]

The code constructs the oriented outer loop and every plaquette on its minimal planar surface as actual directed edges with **four-coordinate lattice endpoints**. The only selected term in the expansion has one opposite-oriented fundamental trace on each of the \(RT\) plaquettes, with coefficient \((\beta/6)^{RT}\). Integrating each once-paired link uses the genuine continuous Haar identity
\(\int_{SU(3)} U_{ab}\overline{U}_{cd}dU=\delta_{ac}\delta_{bd}/3\). Its exact union-find index contraction gives \(E=2RT+R+T\) distinct links and \(V=(R+1)(T+1)\) free color-index loops, whence

\[
\frac{3^V}{3\,6^{RT}\,3^E}=\frac1{18^{RT}}.
\]

The physical gauge theory has **no finite-group replacement**: the integrals use the continuous SU(3) Haar identity. Only the *optional jump diagnostic* reduces rational transfer coefficients modulo 1,000,081.

## Computed results

All nine predicted exact fractions (R,T in 1..3), all nine axis/translation checks, nine signed-area calibration checks, and five failure controls pass. The independent brute-force sum over **all 6,561 raw color-index assignments** on the 1×1 loop gives exactly 81 assignments satisfying the eight Haar Kronecker constraints, independently matching the contraction. The full loop's leading coefficient is 1/18; the 2×2 is 1/104976 and 3×3 is 1/198359290368.

A **three-sector diagonal leading-order transfer** has entries \(1,u^R,u^R\) for trialities 0,+1,-1, where \(u=\beta/18\). Cubing this (modular) transfer k times agrees with power evaluation and direct stepping on the short cases at \(T=3^k\), k=1,2,3,6,10,20. This is an exact jump for the **leading coefficient only**, *not* a demonstration that the full interacting four-dimensional SU(3) transfer operator is diagonal.

Local timing is captured **outside** the pinned result receipt with `--timing`: the 1×1 index brute sum and exact gluing each run locally; neither is an hours-to-seconds production benchmark. Full dimension-4 lattice Monte Carlo and finite-β transfer costs remain OPEN.

## Key caveats for independent reviewer

1. **Minimal-order coverage:** the signed-area/Stokes calibration proves minimal area within ordinary *integer-oriented plaquette 2-chains*. The SU(3) epsilon tensor also allows triple fundamental flux at a link. Independently check the usual strong-coupling selection rule excluding other order-\(RT\) baryonic networks; calibration alone does not prove that fact.
2. **Denominator:** the partition function has constant term 1; nonconstant corrections cannot change the first nonzero Wilson-loop coefficient. Confirm this along with the action normalization and use of normalized `ReTr`.
3. **Open vs periodic:** no wrapping surfaces or global-flux sectors are included; a torus of insufficient size would invalidate the minimal-surface assumption.
4. **No finite β:** the \(O(\beta^{RT+1})\) or stronger remainder here is **not certified**; no physical string tension, static potential at finite coupling, or glueball mass is claimed.
5. **Jump scope:** the three-sector transfer is a reconstructed leading-order area-law transfer. The modular power is a verification of that compressed *coefficient*, not the actual continuous SU(3) physical time evolution.
6. **Integrity vs provenance:** SHA values and the internal seal check bytes only. The outside prediction anchor is the immutable GitHub commit above; the completed result needs an external accepted root pin and Opus' independent check.

## Reproduce

From this directory using only Python 3 standard library:

```
python3 su3_wilson_leading.py > /tmp/SU3-RESULT.hbp
cmp /tmp/SU3-RESULT.hbp RESULT.hbp
python3 verify_su3.py
python3 su3_wilson_leading.py --timing | tail -1
```

The verifier compares every frozen prediction against the measured receipt, rejects missing or changed rows, verifies exact replay, verifies SHA256SUMS and HBI index entries, and checks the internal seal. The five planted faults must be rejected. All receipt lines use `json=0`. The timings are deliberately excluded from byte comparisons.

**Independent reference:** strong-coupling lattice-QCD literature, including Christof Gattringer and Christian B. Lang, *Quantum Chromodynamics on the Lattice*, section on strong-coupling Wilson loops, and CERN [Confinement XVI notes](https://indico.cern.ch/event/353906/contributions/2219959/attachments/1329809/1998182/confinement_16.pdf), with the Wilson coefficient \(u=\beta/18+O(\beta^2)\). This benchmark tests only that published leading coefficient; it is not a new prediction of full lattice QCD.

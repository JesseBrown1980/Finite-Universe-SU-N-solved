//! fgwser: exact strong-coupling series of rectangular Wilson loops for a finite gauge group on the open
//! ROWS x COLS strip. <W(C_A)> = sum_k c_{A,k} r^k, W = chi_f(U_C), C_A the 1 x A rectangle on rows 0..1,
//! columns 0..A, by Rayleigh-Schrodinger perturbation theory in H0 = E, V = M with exact rationals:
//!   E psi_n = -M psi_{n-1} + sum_{k=1..n} E_k psi_{n-k},  E_k = <psi_0|M|psi_{k-1}> / <psi_0|psi_0>,
//!   <psi_0|psi_n> = 0 (n >= 1), psi_0 = constant (the unique E = 0 state);
//! <W> = <psi|W|psi> / <psi|psi> to order ORDER (every <psi_j|W|psi_{k-j}> term and the norm).
//! METHOD 0: gauge-orbit basis under global conjugation with Gauss weights, reduced resolvent by exact
//!   Gauss-Jordan elimination of (W E + w w^T).
//! METHOD 1: full gauge-fixed configuration basis (|G|^free states, unit weights), reduced resolvent by the
//!   spectral polynomial m(x) = prod over every possible nonzero electric eigenvalue (all sums over the links
//!   of the irrep link energies): E^-1 b = -q(E) b / m(0), m(x) = m(0) + x q(x), with an exact residual check.
//! METHOD 2: gauge-orbit basis with Gauss weights and the same spectral polynomial (E is W-self-adjoint, so a
//!   W-orthogonal b lies in the range of E and q(E) b stays there), with the exact residual check.
//! Quadratic fields (2O: sqrt2): M and W take values x + y theta, so the perturbation vectors are pairs
//! p + q theta of rationals; E is an integer matrix, so every solve acts on the two rational components
//! separately. The spectral candidate set is closed under theta -> conjugate, so its polynomial is integer.
//! Values with q = 0 print exactly as rationals "n/d"; otherwise as "n/d+n'/d'<theta>".
//! usage: fgwser GROUP ROWS COLS METHOD ORDER        (GROUP 1 = S3, 2 = Q8, 3 = 2T, 4 = 2O; METHOD 0, 1, 2)        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use fgroup::{basis, group, strip, Alg, Basis, Field, Sys};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeSet;

#[derive(Clone, PartialEq, Debug)]
struct R {
    n: BigInt,
    d: BigInt,
}

impl R {
    fn new(n: BigInt, d: BigInt) -> R {
        let g = n.gcd(&d);
        let (mut n, mut d) = if g.is_zero() { (n, d) } else { (n / &g, d / &g) };
        if d.is_negative() {
            n = -n;
            d = -d;
        }
        R { n, d }
    }
    fn int(v: i64) -> R {
        R { n: BigInt::from(v), d: BigInt::one() }
    }
    fn zero() -> R {
        R::int(0)
    }
    fn add(&self, o: &R) -> R {
        R::new(&self.n * &o.d + &o.n * &self.d, &self.d * &o.d)
    }
    fn sub(&self, o: &R) -> R {
        R::new(&self.n * &o.d - &o.n * &self.d, &self.d * &o.d)
    }
    fn mul(&self, o: &R) -> R {
        R::new(&self.n * &o.n, &self.d * &o.d)
    }
    fn div(&self, o: &R) -> R {
        R::new(&self.n * &o.d, &self.d * &o.n)
    }
    fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    fn show(&self) -> String {
        format!("{}/{}", self.n, self.d)
    }
}

/// p + q theta with rational p, q (theta^2 = fp + fq theta)
#[derive(Clone, PartialEq, Debug)]
struct R2 {
    a: R,
    b: R,
}

impl R2 {
    fn rat(a: R) -> R2 {
        R2 { a, b: R::zero() }
    }
    fn alg(v: Alg) -> R2 {
        R2 { a: R::int(v.x), b: R::int(v.y) }
    }
    fn zero() -> R2 {
        R2::rat(R::zero())
    }
    fn add(&self, o: &R2) -> R2 {
        R2 { a: self.a.add(&o.a), b: self.b.add(&o.b) }
    }
    fn sub(&self, o: &R2) -> R2 {
        R2 { a: self.a.sub(&o.a), b: self.b.sub(&o.b) }
    }
    fn mul(&self, o: &R2, f: Field) -> R2 {
        let bb = self.b.mul(&o.b);
        R2 { a: self.a.mul(&o.a).add(&bb.mul(&R::int(f.p))), b: self.a.mul(&o.b).add(&self.b.mul(&o.a)).add(&bb.mul(&R::int(f.q))) }
    }
    fn mul_r(&self, r: &R) -> R2 {
        R2 { a: self.a.mul(r), b: self.b.mul(r) }
    }
    /// division through the conjugate: (p + q t)^-1 = ((p + q fq) - q t) / (p^2 + p q fq - q^2 fp)
    fn div(&self, o: &R2, f: Field) -> R2 {
        let conj = R2 { a: o.a.add(&o.b.mul(&R::int(f.q))), b: R::zero().sub(&o.b) };
        let norm = o.a.mul(&o.a).add(&o.a.mul(&o.b).mul(&R::int(f.q))).sub(&o.b.mul(&o.b).mul(&R::int(f.p)));
        let num = self.mul(&conj, f);
        R2 { a: num.a.div(&norm), b: num.b.div(&norm) }
    }
    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }
    fn show(&self, f: Field) -> String {
        if self.b.is_zero() {
            self.a.show()
        } else {
            format!("{}{}{}", self.a.show(), if self.b.n.is_negative() { "" } else { "+" }, format_args!("{}{}", self.b.show(), f.theta))
        }
    }
}

trait Solver {
    /// x with E x = b and <psi_0|x> = 0, for b with <psi_0|b> = 0
    fn solve(&self, b: &[R]) -> Option<Vec<R>>;
}

struct Elim {
    inv: Vec<Vec<R>>,
    w: Vec<i64>,
}

impl Solver for Elim {
    fn solve(&self, b: &[R]) -> Option<Vec<R>> {
        let n = b.len();
        let wb: Vec<R> = (0..n).map(|i| b[i].mul(&R::int(self.w[i]))).collect();
        Some((0..n).map(|i| (0..n).fold(R::zero(), |acc, j| acc.add(&self.inv[i][j].mul(&wb[j])))).collect())
    }
}

fn elim(bs: &Basis) -> Option<Elim> {
    let n = bs.reps.len();
    let mut m: Vec<Vec<R>> = vec![vec![R::zero(); n]; n];
    for (i, row) in m.iter_mut().enumerate() {
        for (&j, &v) in &bs.e[i] {
            row[j] = row[j].add(&R::int(bs.w[i] * v));
        }
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = cell.add(&R::int(bs.w[i] * bs.w[j]));
        }
    }
    let mut inv: Vec<Vec<R>> = (0..n).map(|i| (0..n).map(|j| R::int(i64::from(i == j))).collect()).collect();
    for c in 0..n {
        let p = (c..n).find(|&r| !m[r][c].is_zero())?;
        m.swap(c, p);
        inv.swap(c, p);
        let pv = m[c][c].clone();
        for j in 0..n {
            m[c][j] = m[c][j].div(&pv);
            inv[c][j] = inv[c][j].div(&pv);
        }
        for r in 0..n {
            if r != c && !m[r][c].is_zero() {
                let f = m[r][c].clone();
                for j in 0..n {
                    let (a1, b1) = (m[c][j].mul(&f), inv[c][j].mul(&f));
                    m[r][j] = m[r][j].sub(&a1);
                    inv[r][j] = inv[r][j].sub(&b1);
                }
            }
        }
    }
    Some(Elim { inv, w: bs.w.clone() })
}

struct Spectral<'a> {
    e: &'a [std::collections::BTreeMap<usize, i64>],
    /// coefficients a_0..a_D of m(x) = prod (x - lambda)
    a: Vec<BigInt>,
}

impl Solver for Spectral<'_> {
    fn solve(&self, b: &[R]) -> Option<Vec<R>> {
        let n = b.len();
        // common denominator: b = B / den
        let den = b.iter().fold(BigInt::one(), |acc, x| acc.lcm(&x.d));
        let bb: Vec<BigInt> = b.iter().map(|x| &x.n * (&den / &x.d)).collect();
        let d = self.a.len() - 1;
        let mut y: Vec<BigInt> = bb.iter().map(|v| v * &self.a[d]).collect();
        for k in (1..d).rev() {
            let ey: Vec<BigInt> = (0..n).map(|i| self.e[i].iter().fold(BigInt::zero(), |acc, (&j, &v)| acc + &y[j] * v)).collect();
            y = ey.into_iter().zip(bb.iter()).map(|(v, bi)| v + bi * &self.a[k]).collect();
        }
        let scale = &den * &self.a[0];
        let x: Vec<R> = y.into_iter().map(|v| R::new(-v, scale.clone())).collect();
        // exact residual check E x = b
        for (row, bi) in self.e.iter().zip(b.iter()) {
            let ex = row.iter().fold(R::zero(), |acc, (&j, &v)| acc.add(&x[j].mul(&R::int(v))));
            if ex != *bi {
                return None;
            }
        }
        Some(x)
    }
}

fn main() {
    let a: Vec<i64> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    let ok_args = a.len() == 5 && (2..=3).contains(&a[1]) && (2..=4).contains(&a[2]) && (0..=2).contains(&a[3]) && (1..=6).contains(&a[4]);
    let Some(g) = a.first().and_then(|&id| group(id)).filter(|_| ok_args) else {
        eprintln!("usage: fgwser GROUP ROWS COLS METHOD ORDER");
        std::process::exit(1);
    };
    let (rows, cols, orbit, elimination, order) = (a[1] as usize, a[2] as usize, a[3] != 1, a[3] == 0, a[4] as usize);
    let s = strip(rows, cols);
    let sys = Sys { g: &g, s: &s };
    let links: Vec<usize> = (0..s.edges.len()).collect();
    let bs = basis(&sys, orbit, &links);
    let n = bs.reps.len();
    let mut adjoint = true;
    for i in 0..n {
        for (&j, &v) in &bs.e[i] {
            adjoint &= bs.w[i] * v == bs.w[j] * bs.e[j].get(&i).copied().unwrap_or(0);
        }
    }
    let method = match a[3] {
        0 => "A_orbit_elimination",
        1 => "B_config_spectral",
        _ => "A2_orbit_spectral",
    };
    // spectral polynomial over every possible nonzero electric eigenvalue: all sums over the links of the irrep
    // link energies (algebraic integers x + y theta), closed under conjugation; its coefficients are integers
    let f0 = g.field;
    let link_e: BTreeSet<(i64, i64)> = (0..g.chars.len()).filter_map(|k| g.e_irrep(k)).map(|e| (e.x, e.y)).collect();
    let mut sums: BTreeSet<(i64, i64)> = BTreeSet::from([(0, 0)]);
    for _ in 0..links.len() {
        sums = sums.iter().flat_map(|s| link_e.iter().map(move |e| (s.0 + e.0, s.1 + e.1))).collect();
    }
    let conj = |v: (i64, i64)| (v.0 + v.1 * f0.q, -v.1);
    let closed: BTreeSet<(i64, i64)> = sums.iter().flat_map(|&v| [v, conj(v)]).collect();
    let lambdas: Vec<(i64, i64)> = closed.into_iter().filter(|&v| v != (0, 0)).collect();
    // multiply out prod (x - lambda) with coefficients x + y theta in big integers
    let mut poly2: Vec<(BigInt, BigInt)> = vec![(BigInt::one(), BigInt::zero())];
    for &(lx, ly) in &lambdas {
        let mut np = vec![(BigInt::zero(), BigInt::zero()); poly2.len() + 1];
        for (i, c) in poly2.iter().enumerate() {
            np[i + 1].0 += &c.0;
            np[i + 1].1 += &c.1;
            // c * (-lambda) = -(c0 lx + c1 ly p) - (c0 ly + c1 lx + c1 ly q) theta
            np[i].0 -= &c.0 * lx + &c.1 * ly * f0.p;
            np[i].1 -= &c.0 * ly + &c.1 * lx + &c.1 * ly * f0.q;
        }
        poly2 = np;
    }
    let poly_rational = poly2.iter().all(|c| c.1.is_zero());
    let poly: Vec<BigInt> = poly2.into_iter().map(|c| c.0).collect();
    let solver: Box<dyn Solver> = if elimination {
        match elim(&bs) {
            Some(e) => Box::new(e),
            None => {
                println!("FGSERIES|result=FAIL|reason=singular|json=0");
                std::process::exit(2);
            }
        }
    } else {
        if !poly_rational {
            println!("FGSERIES|result=FAIL|reason=spectral_polynomial_not_rational|json=0");
            std::process::exit(2);
        }
        Box::new(Spectral { e: &bs.e, a: poly })
    };
    let f = g.field;
    let mag: Vec<R2> = bs.reps.iter().map(|q| R2::alg(sys.magnetic(&sys.raw(q)))).collect();
    let dot = |x: &[R2], y: &[R2]| -> R2 { (0..n).fold(R2::zero(), |acc, i| acc.add(&x[i].mul(&y[i], f).mul_r(&R::int(bs.w[i])))) };
    let n0 = R::int(bs.w.iter().sum());
    let mut psi: Vec<Vec<R2>> = vec![vec![R2::rat(R::int(1)); n]];
    let mut ek: Vec<R2> = vec![R2::zero()];
    let mut solved = true;
    for k in 1..=order {
        let mp: Vec<R2> = (0..n).map(|i| mag[i].mul(&psi[k - 1][i], f)).collect();
        let e = dot(&psi[0], &mp);
        ek.push(R2 { a: e.a.div(&n0), b: e.b.div(&n0) });
        let mut b: Vec<R2> = mp.iter().map(|v| R2::zero().sub(v)).collect();
        for j in 1..=k {
            for (bi, pv) in b.iter_mut().zip(psi[k - j].iter()) {
                *bi = bi.add(&ek[j].mul(pv, f));
            }
        }
        solved &= dot(&psi[0], &b).is_zero();
        let ba: Vec<R> = b.iter().map(|v| v.a.clone()).collect();
        let bb: Vec<R> = b.iter().map(|v| v.b.clone()).collect();
        match (solver.solve(&ba), solver.solve(&bb)) {
            (Some(xa), Some(xb)) => psi.push(xa.into_iter().zip(xb).map(|(a, b)| R2 { a, b }).collect()),
            _ => {
                solved = false;
                break;
            }
        }
    }
    let ortho = solved && (1..=order).all(|k| dot(&psi[0], &psi[k]).is_zero());
    println!(
        "FGSERIES_BASIS|group={}|rows={rows}|cols={cols}|method={method}|states={n}|e_w_self_adjoint={}|rhs_orthogonal_and_solved={}|orthogonality={}|spectral_degree={}|order={order}|json=0",
        g.name,
        if adjoint { "PASS" } else { "FAIL" },
        if solved { "PASS" } else { "FAIL" },
        if ortho { "PASS" } else { "FAIL" },
        if elimination { 0 } else { lambdas.len() }
    );
    if !(solved && ortho) {
        std::process::exit(2);
    }
    let energies: Vec<String> = ek.iter().skip(1).map(|v| v.show(f)).collect();
    println!("FGSERIES_ENERGY|group={}|E_k={}|json=0", g.name, energies.join(","));
    let dser: Vec<R2> = (0..=order).map(|k| (0..=k).fold(R2::zero(), |acc, j| acc.add(&dot(&psi[j], &psi[k - j])))).collect();
    for x1 in 1..cols {
        let wv: Vec<R2> = bs.reps.iter().map(|q| R2::alg(sys.wilson(&sys.raw(q), 0, x1))).collect();
        let nser: Vec<R2> = (0..=order)
            .map(|k| {
                (0..=k).fold(R2::zero(), |acc, j| {
                    let wy: Vec<R2> = (0..n).map(|i| wv[i].mul(&psi[k - j][i], f)).collect();
                    acc.add(&dot(&psi[j], &wy))
                })
            })
            .collect();
        let mut c: Vec<R2> = Vec::new();
        for k in 0..=order {
            let mut v = nser[k].clone();
            for j in 1..=k {
                v = v.sub(&dser[j].mul(&c[k - j], f));
            }
            c.push(v.div(&dser[0], f));
        }
        let cs: Vec<String> = c.iter().map(|v| v.show(f)).collect();
        println!("FGSERIES|group={}|loop=1x{x1}|method={method}|c_k={}|json=0", g.name, cs.join(","));
    }
}

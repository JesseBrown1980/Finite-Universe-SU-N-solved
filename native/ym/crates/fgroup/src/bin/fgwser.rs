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
//! usage: fgwser GROUP ROWS COLS METHOD ORDER        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use fgroup::{basis, group, strip, Basis, Sys};
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
    let ok_args = a.len() == 5 && (2..=3).contains(&a[1]) && (2..=4).contains(&a[2]) && (0..=1).contains(&a[3]) && (1..=6).contains(&a[4]);
    let Some(g) = a.first().and_then(|&id| group(id)).filter(|_| ok_args) else {
        eprintln!("usage: fgwser GROUP ROWS COLS METHOD ORDER");
        std::process::exit(1);
    };
    let (rows, cols, orbit, order) = (a[1] as usize, a[2] as usize, a[3] == 0, a[4] as usize);
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
    let method = if orbit { "A_orbit_elimination" } else { "B_config_spectral" };
    // spectral polynomial over every possible nonzero electric eigenvalue
    let mut sums: BTreeSet<i64> = BTreeSet::from([0]);
    let link_e: BTreeSet<i64> = (0..g.chars.len()).filter_map(|k| g.e_irrep(k)).collect();
    for _ in 0..links.len() {
        sums = sums.iter().flat_map(|s| link_e.iter().map(move |e| s + e)).collect();
    }
    let lambdas: Vec<i64> = sums.into_iter().filter(|&v| v != 0).collect();
    let solver: Box<dyn Solver> = if orbit {
        match elim(&bs) {
            Some(e) => Box::new(e),
            None => {
                println!("FGSERIES|result=FAIL|reason=singular|json=0");
                std::process::exit(2);
            }
        }
    } else {
        let mut poly = vec![BigInt::one()];
        for &l in &lambdas {
            // multiply by (x - l)
            let mut np = vec![BigInt::zero(); poly.len() + 1];
            for (i, c) in poly.iter().enumerate() {
                np[i + 1] += c;
                np[i] -= c * l;
            }
            poly = np;
        }
        Box::new(Spectral { e: &bs.e, a: poly })
    };
    let mag: Vec<R> = bs.reps.iter().map(|q| R::int(sys.magnetic(&sys.raw(q)))).collect();
    let dot = |x: &[R], y: &[R]| -> R { (0..n).fold(R::zero(), |acc, i| acc.add(&x[i].mul(&y[i]).mul(&R::int(bs.w[i])))) };
    let n0 = R::int(bs.w.iter().sum());
    let mut psi: Vec<Vec<R>> = vec![vec![R::int(1); n]];
    let mut ek: Vec<R> = vec![R::zero()];
    let mut solved = true;
    for k in 1..=order {
        let mp: Vec<R> = (0..n).map(|i| mag[i].mul(&psi[k - 1][i])).collect();
        ek.push(dot(&psi[0], &mp).div(&n0));
        let mut b: Vec<R> = mp.iter().map(|v| R::zero().sub(v)).collect();
        for j in 1..=k {
            for (bi, pv) in b.iter_mut().zip(psi[k - j].iter()) {
                *bi = bi.add(&ek[j].mul(pv));
            }
        }
        solved &= dot(&psi[0], &b).is_zero();
        match solver.solve(&b) {
            Some(x) => psi.push(x),
            None => {
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
        if orbit { 0 } else { lambdas.len() }
    );
    if !(solved && ortho) {
        std::process::exit(2);
    }
    let energies: Vec<String> = ek.iter().skip(1).map(R::show).collect();
    println!("FGSERIES_ENERGY|group={}|E_k={}|json=0", g.name, energies.join(","));
    let dser: Vec<R> = (0..=order).map(|k| (0..=k).fold(R::zero(), |acc, j| acc.add(&dot(&psi[j], &psi[k - j])))).collect();
    for x1 in 1..cols {
        let wv: Vec<R> = bs.reps.iter().map(|q| R::int(sys.wilson(&sys.raw(q), 0, x1))).collect();
        let nser: Vec<R> = (0..=order)
            .map(|k| {
                (0..=k).fold(R::zero(), |acc, j| {
                    let wy: Vec<R> = (0..n).map(|i| wv[i].mul(&psi[k - j][i])).collect();
                    acc.add(&dot(&psi[j], &wy))
                })
            })
            .collect();
        let mut c: Vec<R> = Vec::new();
        for k in 0..=order {
            let mut v = nser[k].clone();
            for j in 1..=k {
                v = v.sub(&dser[j].mul(&c[k - j]));
            }
            c.push(v.div(&dser[0]));
        }
        let cs: Vec<String> = c.iter().map(R::show).collect();
        println!("FGSERIES|group={}|loop=1x{x1}|method={method}|c_k={}|json=0", g.name, cs.join(","));
    }
}

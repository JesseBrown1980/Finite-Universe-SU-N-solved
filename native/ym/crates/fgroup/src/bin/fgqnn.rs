//! fgqnn: the qutrit variational network (v7 of GGUF-files-operating-as-trained-quantum-neural-network-models)
//! on a finite-group gauge theory on an open ROWS x COLS strip. Each free link is k qutrits (3^k >= |G|);
//! codes >= |G| are discarded. Gate R_ij(b) on one qutrit: [[a^2-b^2, -2ab], [2ab, a^2-b^2]] on levels i, j and
//! a^2+b^2 on the third, a = 6561, b from eight balanced trits. Layer 1: R_01, R_02 on every qutrit. Layer 2:
//! SUM ring, R_01, R_02 on every qutrit, inverse ring (the identity when its trits are zero). The state is
//! projected onto conjugation orbits, x_i = (|G| / w_i) * (orbit sum), and the energy is the exact rational
//! (RDEN sum_i w_i x_i (E x)_i + RNUM sum_i w_i M_i x_i^2) / (RDEN sum_i w_i x_i^2) with E, M, w from
//! fgroup::basis, H/J_E = E + r M. Training: coordinate descent over the trits, strict decreases kept.
//! usage: fgqnn GROUP ROWS COLS RNUM RDEN SEED SWEEPS     (SEED 0 = UNIFORM, 1 = IDENTITY)        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use fgroup::{basis, group, strip, Basis, Sys};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

const A: i64 = 6561;
const TPB: usize = 8;

fn b8(p: &[i64]) -> i64 {
    p.iter().fold(0i64, |acc, &x| acc * 3 + x)
}

fn trits8(b: i64) -> [i64; TPB] {
    let mut out = [0i64; TPB];
    let mut x = b;
    for t in out.iter_mut().rev() {
        let r = x.rem_euclid(3);
        let d = if r == 2 { -1 } else { r };
        *t = d;
        x = (x - d) / 3;
    }
    out
}

struct St {
    v: Vec<BigInt>,
    scale: BigInt,
}

fn pw(q: usize) -> usize {
    3usize.pow(q as u32)
}

impl St {
    fn zero(nq: usize) -> St {
        let mut v = vec![BigInt::zero(); pw(nq)];
        v[0] = BigInt::one();
        St { v, scale: BigInt::one() }
    }
    fn rotate(&mut self, q: usize, i: usize, j: usize, b: i64) {
        if b == 0 {
            return;
        }
        let k = 3 - i - j;
        let (c, s, d) = (BigInt::from(A * A - b * b), BigInt::from(2 * A * b), BigInt::from(A * A + b * b));
        let st = pw(q);
        for base in 0..self.v.len() {
            if (base / st) % 3 != 0 {
                continue;
            }
            let (ii, jj, kk) = (base + i * st, base + j * st, base + k * st);
            let ni = &self.v[ii] * &c - &self.v[jj] * &s;
            let nj = &self.v[ii] * &s + &self.v[jj] * &c;
            self.v[kk] = &self.v[kk] * &d;
            self.v[ii] = ni;
            self.v[jj] = nj;
        }
        self.scale = &self.scale * &d;
    }
    fn sum(&mut self, c: usize, t: usize, s: i64) {
        let (sc, st) = (pw(c), pw(t));
        let mut out = vec![BigInt::zero(); self.v.len()];
        for (x, z) in self.v.iter().enumerate() {
            let dc = ((x / sc) % 3) as i64;
            let dt = (x / st) % 3;
            let nt = (dt as i64 + s * dc).rem_euclid(3) as usize;
            out[x - dt * st + nt * st] = z.clone();
        }
        self.v = out;
    }
}

struct Net<'a> {
    p: &'a [i64],
    nq: usize,
}

impl Net<'_> {
    fn bs(&self, l: usize, q: usize) -> (i64, i64) {
        let i = (l * self.nq + q) * 2 * TPB;
        (b8(&self.p[i..i + TPB]), b8(&self.p[i + TPB..i + 2 * TPB]))
    }
    fn layer2_active(&self) -> bool {
        self.p[self.nq * 2 * TPB..].iter().any(|&x| x != 0)
    }
    fn rot(&self, l: usize, s: &mut St, sign: i64) {
        if sign > 0 {
            for q in 0..self.nq {
                let (b1, b2) = self.bs(l, q);
                s.rotate(q, 0, 1, b1);
                s.rotate(q, 0, 2, b2);
            }
        } else {
            for q in (0..self.nq).rev() {
                let (b1, b2) = self.bs(l, q);
                s.rotate(q, 0, 2, -b2);
                s.rotate(q, 0, 1, -b1);
            }
        }
    }
    fn ring(&self, s: &mut St, sign: i64) {
        if sign > 0 {
            for q in 0..self.nq {
                s.sum(q, (q + 1) % self.nq, 1);
            }
        } else {
            for q in (0..self.nq).rev() {
                s.sum(q, (q + 1) % self.nq, -1);
            }
        }
    }
    fn layer1(&self) -> St {
        let mut s = St::zero(self.nq);
        self.rot(0, &mut s, 1);
        s
    }
    fn layer2(&self, s: &mut St) {
        if self.layer2_active() {
            self.ring(s, 1);
            self.rot(1, s, 1);
            self.ring(s, -1);
        }
    }
    fn inverse(&self, s: &mut St) {
        if self.layer2_active() {
            self.ring(s, 1);
            self.rot(1, s, -1);
            self.ring(s, -1);
        }
        self.rot(0, s, -1);
    }
}

struct Model<'a> {
    b: &'a Basis,
    n: i64,
    /// orbit of every register index, None for an invalid code
    orbit: Vec<Option<usize>>,
    mag: Vec<i64>,
    rnum: i64,
    rden: i64,
}

impl Model<'_> {
    fn energy(&self, s: &St) -> Option<(BigInt, BigInt)> {
        let m = self.b.reps.len();
        let mut sums = vec![BigInt::zero(); m];
        for (x, z) in s.v.iter().enumerate() {
            if let Some(o) = self.orbit[x] {
                if !z.is_zero() {
                    sums[o] += z;
                }
            }
        }
        let xs: Vec<BigInt> = sums.iter().enumerate().map(|(o, v)| v * BigInt::from(self.n / self.b.w[o])).collect();
        let mut den = BigInt::zero();
        let mut em = BigInt::zero();
        let mut ee = BigInt::zero();
        for o in 0..m {
            if xs[o].is_zero() {
                continue;
            }
            let w = BigInt::from(self.b.w[o]);
            let sq = &xs[o] * &xs[o] * &w;
            em += &sq * BigInt::from(self.mag[o]);
            den += &sq;
            let mut ex = BigInt::zero();
            for (&j, &v) in &self.b.e[o] {
                if !xs[j].is_zero() {
                    ex += &xs[j] * BigInt::from(v);
                }
            }
            ee += &xs[o] * ex * w;
        }
        if den.is_zero() {
            return None;
        }
        Some((ee * BigInt::from(self.rden) + em * BigInt::from(self.rnum), den * BigInt::from(self.rden)))
    }
}

fn less(x: &(BigInt, BigInt), y: &(BigInt, BigInt)) -> bool {
    &x.0 * &y.1 < &y.0 * &x.1
}

fn milli(x: &(BigInt, BigInt)) -> BigInt {
    (&x.0 * BigInt::from(1000)).div_floor(&x.1)
}

fn main() {
    let args: Vec<i64> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    if args.len() != 7 {
        eprintln!("usage: fgqnn GROUP ROWS COLS RNUM RDEN SEED SWEEPS");
        std::process::exit(1);
    }
    let (gid, rows, cols, rnum, rden, seed, sweeps) = (args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
    let ok = (2..=3).contains(&rows) && (2..=4).contains(&cols) && rnum >= 0 && rden > 0 && (0..=1).contains(&seed) && (1..=8).contains(&sweeps);
    let Some(g) = group(gid).filter(|_| ok) else {
        eprintln!("FG_QNN|result=FAIL|reason=out_of_range|json=0");
        std::process::exit(1);
    };
    let s = strip(rows as usize, cols as usize);
    let sys = Sys { g: &g, s: &s };
    let links: Vec<usize> = (0..s.edges.len()).collect();
    let b = basis(&sys, true, &links);
    let n = g.n;
    let k = (1..).find(|&k| pw(k) >= n).unwrap_or(1);
    let nf = s.free.len();
    let nq = nf * k;
    let dim = pw(nq);
    let orbit: Vec<Option<usize>> = (0..dim)
        .map(|x| {
            let q: Vec<u8> = (0..nf).map(|i| (0..k).fold(0usize, |acc, j| acc + ((x / pw(i * k + j)) % 3) * pw(j))).map(|v| v as u8).collect();
            if q.iter().all(|&v| (v as usize) < n) {
                b.index.get(&sys.key(&q)).copied()
            } else {
                None
            }
        })
        .collect();
    let mag: Vec<i64> = b.reps.iter().map(|q| sys.magnetic(&sys.raw(q))).collect();
    let model = Model { b: &b, n: n as i64, orbit, mag, rnum, rden };
    let count = 2 * nq * 2 * TPB;
    let mut p = vec![0i64; count];
    if seed == 0 {
        let (u1, u2) = (trits8(2065), trits8(2718));
        for q in 0..nq {
            let i = q * 2 * TPB;
            p[i..i + TPB].copy_from_slice(&u1);
            p[i + TPB..i + 2 * TPB].copy_from_slice(&u2);
        }
    }
    let seed_p = p.clone();
    let half = nq * 2 * TPB;
    let eval_from = |p: &[i64], l1: &St| {
        let net = Net { p, nq };
        let mut s = St { v: l1.v.clone(), scale: l1.scale.clone() };
        net.layer2(&mut s);
        model.energy(&s)
    };
    let mut cache = Net { p: &p, nq }.layer1();
    let start = eval_from(&p, &cache).expect("physical seed");
    let mut best = start.clone();
    let mut evals = 1usize;
    let mut used = 0usize;
    for _ in 0..sweeps {
        used += 1;
        let mut changed = false;
        for i in 0..count {
            let keep = p[i];
            let mut best_v = keep;
            for v in -1..=1 {
                if v == keep {
                    continue;
                }
                p[i] = v;
                evals += 1;
                let e = if i < half {
                    let l1 = Net { p: &p, nq }.layer1();
                    eval_from(&p, &l1)
                } else {
                    eval_from(&p, &cache)
                };
                if let Some(e) = e {
                    if less(&e, &best) {
                        best = e;
                        best_v = v;
                    }
                }
            }
            p[i] = best_v;
            if best_v != keep {
                changed = true;
                if i < half {
                    cache = Net { p: &p, nq }.layer1();
                }
            }
        }
        if !changed {
            break;
        }
    }
    // exact norm and inverse of the trained network
    let net = Net { p: &p, nq };
    let mut st = net.layer1();
    net.layer2(&mut st);
    let d = st.scale.clone();
    let norm_ok = st.v.iter().fold(BigInt::zero(), |a, z| a + z * z) == &d * &d;
    let fwd_scale = st.scale.clone();
    net.inverse(&mut st);
    let dd = &fwd_scale * &fwd_scale;
    let inv_ok = st.v.iter().enumerate().all(|(x, z)| if x == 0 { *z == dd } else { z.is_zero() }) && st.scale == dd;
    let again = { let mut s2 = net.layer1(); net.layer2(&mut s2); model.energy(&s2) };
    let again_ok = again.as_ref().is_some_and(|e| e.0 == best.0 && e.1 == best.1);
    let changed = p.iter().zip(seed_p.iter()).filter(|(a, b)| a != b).count();
    let bs: Vec<String> = p.chunks(TPB).map(|c| b8(c).to_string()).collect();
    let trits: String = p.iter().map(|&t| match t { -1 => '-', 0 => '0', _ => '+' }).collect();
    let sign_ok = !best.0.is_negative() && best.1.is_positive();
    println!(
        "FG_QNN|group={}|rows={rows}|cols={cols}|r={rnum}/{rden}|seed={}|qutrits_per_link={k}|qutrits={nq}|params={count}|states={}|sweeps={used}|evaluations={evals}|params_changed={changed}|start_milli_floor={}|energy_milli_floor={}|norm_exact={}|inverse_exact={}|recomputed_equal={}|nonnegative={}|json=0",
        g.name,
        if seed == 0 { "UNIFORM" } else { "IDENTITY" },
        b.reps.len(),
        milli(&start),
        milli(&best),
        u8::from(norm_ok),
        u8::from(inv_ok),
        u8::from(again_ok),
        u8::from(sign_ok)
    );
    println!("FG_QNN_START|numerator={}|denominator={}|json=0", start.0, start.1);
    println!("FG_QNN_ENERGY|numerator={}|denominator={}|json=0", best.0, best.1);
    println!("FG_QNN_B|layer1={}|layer2={}|json=0", bs[..2 * nq].join(","), bs[2 * nq..].join(","));
    println!("FG_QNN_TRITS|trits={trits}|json=0");
    std::process::exit(if norm_ok && inv_ok && again_ok && sign_ok { 0 } else { 2 });
}

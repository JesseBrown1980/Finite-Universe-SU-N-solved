//! sturmq2: exact inertia count over a real quadratic field Q(theta), theta^2 = P + Q theta, for the
//! finite-group Hamiltonians whose diagonal carries algebraic values (2O: theta = sqrt2).
//! B(T) = W (1000 (JE Z + JB M) - T W), with Z the integer electric matrix and M the diagonal (x + y theta),
//! built exactly symmetric from the exchange file (upper triples mirrored, as in sturmrs; the builder has checked
//! w_i Z_ij = w_j Z_ji before writing). Reads the integer format (y = 0) and the "QF P Q" format of fgconf.
//! Inertia by exact LDL^T with symmetric pivoting: a nonzero diagonal pivot when one exists (fewest nonzero
//! entries in its row, then lowest index), otherwise an explicit 2 x 2 block on the first nonzero off-diagonal
//! pair (determinant -a^2 < 0: one negative and one positive eigenvalue); a remaining all-zero block counts as
//! zero eigenvalues. Sylvester's law of inertia makes the counts those of B. Signs of a + b theta are exact:
//! a + b theta = u/2 + (b/2) sqrt D with u = 2a + bQ, D = Q^2 + 4P > 0, compared through u^2 against D b^2.
//! No intervals and no floats; rationals in big integers.
//! MODE interval (z2o addendum 1): the sturmrs algorithm unchanged (entries as outward-rounded intervals of
//! integers scaled by S = 10^80, Schur updates enclosed by floor and ceiling, a pivot only when its interval
//! excludes zero, COMPLETE only when every pivot is decided, otherwise UNCERTAIN), with each theta part entering
//! through sqrt2 in [L, L + 1] / S, L = isqrt(2 S^2), or phi in [L, L + 1] / S, L = floor((S + isqrt(5 S^2)) / 2)
//! (z2i; those two fields). For an integer matrix it is the same computation as sturmrs.
//! MODE signcheck (z2i S2): the exact sign of a + b phi against a rational enclosure of phi, on a grid of small pairs
//! and on consecutive Fibonacci pairs, where a + b phi is smallest in size.
//! usage: sturmq2 MATRIX JE JB T [exact|interval]        sturmq2 signcheck        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use num_bigint::BigInt;
use num_integer::Integer;
use std::collections::HashMap;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;
use std::time::Instant;

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
    fn int(v: BigInt) -> R {
        R { n: v, d: BigInt::one() }
    }
    fn zero() -> R {
        R::int(BigInt::zero())
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
    fn bits(&self) -> u64 {
        self.n.bits().max(self.d.bits())
    }
}

#[derive(Clone, Copy)]
struct Fld {
    p: i64,
    q: i64,
}

/// a + b theta
#[derive(Clone, PartialEq, Debug)]
struct E {
    a: R,
    b: R,
}

impl E {
    fn zero() -> E {
        E { a: R::zero(), b: R::zero() }
    }
    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }
    fn add(&self, o: &E) -> E {
        E { a: self.a.add(&o.a), b: self.b.add(&o.b) }
    }
    fn sub(&self, o: &E) -> E {
        E { a: self.a.sub(&o.a), b: self.b.sub(&o.b) }
    }
    fn mul(&self, o: &E, f: Fld) -> E {
        let bb = self.b.mul(&o.b);
        E { a: self.a.mul(&o.a).add(&bb.mul(&R::int(BigInt::from(f.p)))), b: self.a.mul(&o.b).add(&self.b.mul(&o.a)).add(&bb.mul(&R::int(BigInt::from(f.q)))) }
    }
    fn div(&self, o: &E, f: Fld) -> E {
        let conj = E { a: o.a.add(&o.b.mul(&R::int(BigInt::from(f.q)))), b: R::zero().sub(&o.b) };
        let norm = o.a.mul(&o.a).add(&o.a.mul(&o.b).mul(&R::int(BigInt::from(f.q)))).sub(&o.b.mul(&o.b).mul(&R::int(BigInt::from(f.p))));
        let num = self.mul(&conj, f);
        E { a: num.a.div(&norm), b: num.b.div(&norm) }
    }
    /// exact sign of a + b theta, theta the positive root of t^2 = P + Q t
    fn sign(&self, f: Fld) -> Ordering {
        // u = 2a + bQ, v = b: value = (u + v sqrt D) / 2
        let u = self.a.add(&self.a).add(&self.b.mul(&R::int(BigInt::from(f.q))));
        let v = self.b.clone();
        let su = u.n.sign();
        let sv = v.n.sign();
        use num_bigint::Sign::{Minus, NoSign, Plus};
        match (su, sv) {
            (NoSign, NoSign) => Ordering::Equal,
            (Plus, Plus) | (Plus, NoSign) | (NoSign, Plus) => Ordering::Greater,
            (Minus, Minus) | (Minus, NoSign) | (NoSign, Minus) => Ordering::Less,
            _ => {
                // compare u^2 with D v^2
                let d = R::int(BigInt::from(f.q * f.q + 4 * f.p));
                let lhs = u.mul(&u);
                let rhs = d.mul(&v).mul(&v);
                let c = (&lhs.n * &rhs.d).cmp(&(&rhs.n * &lhs.d));
                match (su, c) {
                    (Plus, Ordering::Greater) | (Minus, Ordering::Less) => Ordering::Greater,
                    (Plus, Ordering::Less) | (Minus, Ordering::Greater) => Ordering::Less,
                    _ => Ordering::Equal,
                }
            }
        }
    }
    fn bits(&self) -> u64 {
        self.a.bits().max(self.b.bits())
    }
}

struct Input {
    f: Fld,
    b: Vec<Vec<E>>,
    entries: usize,
    algebraic: bool,
}

fn load(path: &str, je: i64, jb: i64, tn: i64) -> Result<Input, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut toks = text.split_ascii_whitespace().peekable();
    let mut f = Fld { p: 2, q: 0 };
    let mut algebraic = false;
    if toks.peek() == Some(&"QF") {
        toks.next();
        let p = toks.next().and_then(|t| t.parse().ok()).ok_or("bad QF P")?;
        let q = toks.next().and_then(|t| t.parse().ok()).ok_or("bad QF Q")?;
        f = Fld { p, q };
        algebraic = true;
        if f.q * f.q + 4 * f.p <= 0 {
            return Err("field is not real".into());
        }
    }
    let int = |t: Option<&str>| -> Result<i64, String> { t.and_then(|v| v.parse::<i64>().ok()).ok_or_else(|| "bad integer token".to_string()) };
    let n = usize::try_from(int(toks.next())?).map_err(|_| "bad n")?;
    let w: Vec<i64> = (0..n).map(|_| int(toks.next())).collect::<Result<_, _>>()?;
    if w.iter().any(|&v| v <= 0) {
        return Err("non-positive weight".into());
    }
    let mut mag: Vec<(i64, i64)> = Vec::with_capacity(n);
    for _ in 0..n {
        let t = toks.next().ok_or("truncated diagonal")?;
        let v = match t.split_once(':') {
            Some((x, y)) if algebraic => (x.parse::<i64>().map_err(|_| "bad X")?, y.parse::<i64>().map_err(|_| "bad Y")?),
            Some(_) => return Err("algebraic diagonal without QF header".into()),
            None => (t.parse::<i64>().map_err(|_| "bad diagonal")?, 0),
        };
        mag.push(v);
    }
    let rest: Vec<i64> = toks.map(|t| t.parse::<i64>().map_err(|_| format!("bad token {t}"))).collect::<Result<_, _>>()?;
    if rest.len() % 3 != 0 {
        return Err("malformed triples".into());
    }
    let mut b = vec![vec![E::zero(); n]; n];
    let mut seen = std::collections::HashSet::new();
    for t in rest.chunks_exact(3) {
        let (i, j, z) = (usize::try_from(t[0]).map_err(|_| "bad i")?, usize::try_from(t[1]).map_err(|_| "bad j")?, t[2]);
        if i > j || j >= n || !seen.insert((i, j)) {
            return Err("bad entry".into());
        }
        let v = R::int(BigInt::from(1000) * je * w[i] * z);
        b[i][j] = b[i][j].add(&E { a: v.clone(), b: R::zero() });
        if i != j {
            b[j][i] = b[j][i].add(&E { a: v, b: R::zero() });
        }
    }
    for k in 0..n {
        let da = BigInt::from(1000) * jb * mag[k].0 * w[k] - BigInt::from(tn) * w[k];
        let db = BigInt::from(1000) * jb * mag[k].1 * w[k];
        b[k][k] = b[k][k].add(&E { a: R::int(da), b: R::int(db) });
    }
    Ok(Input { f, b, entries: seen.len(), algebraic })
}

/// (negative, positive, zero, two_by_two_blocks, max_entry_bits)
fn inertia(mut a: Vec<Vec<E>>, f: Fld) -> (usize, usize, usize, usize, u64) {
    let n = a.len();
    let mut active: Vec<usize> = (0..n).collect();
    let (mut neg, mut pos, mut zero, mut blocks, mut bits) = (0usize, 0usize, 0usize, 0usize, 0u64);
    while !active.is_empty() {
        // a nonzero diagonal pivot: fewest nonzero entries in its active row, then lowest index
        let mut best: Option<(usize, usize)> = None;
        for &k in &active {
            if !a[k][k].is_zero() {
                let nz = active.iter().filter(|&&j| !a[k][j].is_zero()).count();
                if best.map_or(true, |(_, c)| nz < c) {
                    best = Some((k, nz));
                }
            }
        }
        if let Some((k, _)) = best {
            let d = a[k][k].clone();
            match d.sign(f) {
                Ordering::Less => neg += 1,
                _ => pos += 1,
            }
            active.retain(|&x| x != k);
            let col: Vec<(usize, E)> = active.iter().filter(|&&i| !a[i][k].is_zero()).map(|&i| (i, a[i][k].div(&d, f))).collect();
            for (i, fac) in &col {
                for &j in &active {
                    if !a[k][j].is_zero() {
                        let upd = fac.mul(&a[k][j], f);
                        a[*i][j] = a[*i][j].sub(&upd);
                        bits = bits.max(a[*i][j].bits());
                    }
                }
            }
            continue;
        }
        // no nonzero diagonal: a 2 x 2 block on the first nonzero off-diagonal pair
        let mut pair: Option<(usize, usize)> = None;
        'outer: for (x, &i) in active.iter().enumerate() {
            for &j in &active[x + 1..] {
                if !a[i][j].is_zero() {
                    pair = Some((i, j));
                    break 'outer;
                }
            }
        }
        let Some((i, j)) = pair else {
            zero += active.len();
            break;
        };
        let av = a[i][j].clone();
        neg += 1;
        pos += 1;
        blocks += 1;
        active.retain(|&x| x != i && x != j);
        let rows: Vec<(usize, E, E)> = active.iter().map(|&r| (r, a[r][i].clone(), a[r][j].clone())).collect();
        for (r, ari, arj) in &rows {
            for &s in &active {
                let ajs = a[j][s].clone();
                let ais = a[i][s].clone();
                if (ari.is_zero() || ajs.is_zero()) && (arj.is_zero() || ais.is_zero()) {
                    continue;
                }
                let upd = ari.mul(&ajs, f).add(&arj.mul(&ais, f)).div(&av, f);
                a[*r][s] = a[*r][s].sub(&upd);
                bits = bits.max(a[*r][s].bits());
            }
        }
    }
    (neg, pos, zero, blocks, bits)
}

type Iv = (BigInt, BigInt);

fn cdiv(a: &BigInt, b: &BigInt) -> BigInt {
    -((-a).div_floor(b))
}

/// enclosure of x y / d (the sturmrs rule)
fn quotient(x: &Iv, y: &Iv, d: &Iv) -> Result<Iv, String> {
    if !d.0.is_positive() && !d.1.is_negative() {
        return Err("pivot interval contains zero".into());
    }
    let p = [&x.0 * &y.0, &x.0 * &y.1, &x.1 * &y.0, &x.1 * &y.1];
    let mut lo = p.iter().min().cloned().unwrap_or_default();
    let mut hi = p.iter().max().cloned().unwrap_or_default();
    let (mut dlo, mut dhi) = (d.0.clone(), d.1.clone());
    if dhi.is_negative() {
        let (l, h) = (-hi, -lo);
        lo = l;
        hi = h;
        let (a, b) = (-dhi, -dlo);
        dlo = a;
        dhi = b;
    }
    Ok(if !lo.is_negative() {
        (lo.div_floor(&dhi), cdiv(&hi, &dlo))
    } else if !hi.is_positive() {
        (lo.div_floor(&dlo), cdiv(&hi, &dhi))
    } else {
        (lo.div_floor(&dlo), cdiv(&hi, &dlo))
    })
}

/// interval rows as in sturmrs (scale S = 10^80); theta parts through sqrt2 in [L, L + 1] / S
type IntervalInput = (Vec<HashMap<usize, Iv>>, usize, Option<(i64, i64)>);

fn load_interval(path: &str, je: i64, jb: i64, tn: i64) -> Result<IntervalInput, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut toks = text.split_ascii_whitespace().peekable();
    let mut field = None;
    if toks.peek() == Some(&"QF") {
        toks.next();
        let p: i64 = toks.next().and_then(|t| t.parse().ok()).ok_or("bad QF P")?;
        let q: i64 = toks.next().and_then(|t| t.parse().ok()).ok_or("bad QF Q")?;
        if (p, q) != (2, 0) && (p, q) != (1, 1) {
            return Err("interval mode supports the sqrt2 and phi fields only".into());
        }
        field = Some((p, q));
    }
    let algebraic = field.is_some();
    let int = |t: Option<&str>| -> Result<i64, String> { t.and_then(|v| v.parse::<i64>().ok()).ok_or_else(|| "bad integer token".to_string()) };
    let n = usize::try_from(int(toks.next())?).map_err(|_| "bad n")?;
    let w: Vec<i64> = (0..n).map(|_| int(toks.next())).collect::<Result<_, _>>()?;
    if w.iter().any(|&v| v <= 0) {
        return Err("non-positive weight".into());
    }
    let mut mag: Vec<(i64, i64)> = Vec::with_capacity(n);
    for _ in 0..n {
        let t = toks.next().ok_or("truncated diagonal")?;
        let v = match t.split_once(':') {
            Some((x, y)) if algebraic => (x.parse::<i64>().map_err(|_| "bad X")?, y.parse::<i64>().map_err(|_| "bad Y")?),
            Some(_) => return Err("algebraic diagonal without QF header".into()),
            None => (t.parse::<i64>().map_err(|_| "bad diagonal")?, 0),
        };
        mag.push(v);
    }
    let rest: Vec<i64> = toks.map(|t| t.parse::<i64>().map_err(|_| format!("bad token {t}"))).collect::<Result<_, _>>()?;
    if rest.len() % 3 != 0 {
        return Err("malformed triples".into());
    }
    let scale = num_traits::pow(BigInt::from(10), 80);
    let l = if field == Some((1, 1)) {
        phi_lower(&scale)?
    } else {
        let two_s2 = BigInt::from(2) * &scale * &scale;
        let l = two_s2.sqrt();
        // the enclosure sqrt2 in [L, L + 1] / S holds exactly: L^2 <= 2 S^2 < (L + 1)^2
        if !(&l * &l <= two_s2 && (&l + 1) * (&l + 1) > two_s2) {
            return Err("sqrt2 enclosure check failed".into());
        }
        l
    };
    let mut a: Vec<HashMap<usize, Iv>> = vec![HashMap::new(); n];
    let mut entries = 0usize;
    for t in rest.chunks_exact(3) {
        let (i, j, z) = (usize::try_from(t[0]).map_err(|_| "bad i")?, usize::try_from(t[1]).map_err(|_| "bad j")?, t[2]);
        if i > j || j >= n || a[i].contains_key(&j) {
            return Err("bad entry".into());
        }
        entries += 1;
        if z != 0 {
            let v = BigInt::from(1000) * je * w[i] * z * &scale;
            a[i].insert(j, (v.clone(), v.clone()));
            a[j].insert(i, (v.clone(), v));
        }
    }
    for k in 0..n {
        let dv = (BigInt::from(1000) * jb * mag[k].0 * w[k] - BigInt::from(tn) * w[k]) * &scale;
        let b = BigInt::from(1000) * jb * mag[k].1 * w[k];
        let (blo, bhi) = if b.is_negative() { (&b * (&l + 1), &b * &l) } else { (&b * &l, &b * (&l + 1)) };
        let old = a[k].get(&k).cloned().unwrap_or((BigInt::zero(), BigInt::zero()));
        a[k].insert(k, (old.0 + &dv + blo, old.1 + &dv + bhi));
    }
    Ok((a, entries, field))
}

/// L = floor(phi S), checked exactly: phi S is the positive root of x^2 - S x - S^2, so L^2 - S L - S^2 <= 0 < (L + 1)^2 - S (L + 1) - S^2
fn phi_lower(scale: &BigInt) -> Result<BigInt, String> {
    let five_s2 = BigInt::from(5) * scale * scale;
    let l: BigInt = (scale + five_s2.sqrt()) / 2;
    let value = |x: &BigInt| x * x - scale * x - scale * scale;
    if value(&l).is_positive() || !value(&(&l + 1)).is_positive() {
        return Err("phi enclosure check failed".into());
    }
    Ok(l)
}

/// z2i S2: E::sign in the phi field against an enclosure of phi: the sign the enclosure decides must equal it, and every
/// nonzero pair must be decided. Grid |a|, |b| <= 64; then both signs of consecutive Fibonacci pairs (F_{n+1}, -F_n) up to 2^4096.
fn signcheck() -> Result<String, String> {
    let f = Fld { p: 1, q: 1 };
    let bits = 8300u32;
    let scale = num_traits::pow(BigInt::from(2), bits as usize);
    let l = phi_lower(&scale)?;
    let (mut disagree, mut undecided) = (0u64, 0u64);
    let mut judge = |a: &BigInt, b: &BigInt| {
        let (lo, hi) = if b.is_negative() { (a * &scale + b * (&l + 1), a * &scale + b * &l) } else { (a * &scale + b * &l, a * &scale + b * (&l + 1)) };
        let want = if a.is_zero() && b.is_zero() {
            Some(Ordering::Equal)
        } else if lo.is_positive() {
            Some(Ordering::Greater)
        } else if hi.is_negative() {
            Some(Ordering::Less)
        } else {
            None
        };
        let got = E { a: R::int(a.clone()), b: R::int(b.clone()) }.sign(f);
        match want {
            Some(w) if w == got => {}
            Some(_) => disagree += 1,
            None => undecided += 1,
        }
    };
    let bound = 64i64;
    let mut grid = 0u64;
    for a in -bound..=bound {
        for b in -bound..=bound {
            judge(&BigInt::from(a), &BigInt::from(b));
            grid += 1;
        }
    }
    let limit = num_traits::pow(BigInt::from(2), 4096);
    let (mut lo, mut hi) = (BigInt::one(), BigInt::one());
    let mut fib = 0u64;
    let mut max_bits = 0u64;
    while hi <= limit {
        judge(&hi, &(-&lo));
        judge(&(-&hi), &lo);
        fib += 2;
        max_bits = hi.bits();
        let next = &hi + &lo;
        lo = hi;
        hi = next;
    }
    let ok = disagree == 0 && undecided == 0;
    Ok(format!(
        "STURM-Q2-S2|result={}|field=QF_1_1|grid_bound={bound}|grid_pairs={grid}|fibonacci_pairs={fib}|fibonacci_max_bits={max_bits}|enclosure_bits={bits}|disagreements={disagree}|undecided={undecided}|json=0",
        if ok { "PASS" } else { "FAIL" }
    ))
}

/// the sturmrs interval LDL^T: (negative, positive, pivots, complete)
fn certify(a: &mut [HashMap<usize, Iv>]) -> Result<(usize, usize, usize, bool), String> {
    let n = a.len();
    let mut active: Vec<bool> = vec![true; n];
    let (mut neg, mut pos) = (0usize, 0usize);
    for step in 0..n {
        let mut best: Option<(usize, (usize, BigInt))> = None;
        for (i, row) in a.iter().enumerate() {
            if !active[i] {
                continue;
            }
            let Some(d) = row.get(&i) else { continue };
            if !d.0.is_positive() && !d.1.is_negative() {
                continue;
            }
            let key = (row.len(), &d.1 - &d.0);
            if best.as_ref().map_or(true, |(_, k)| key < *k) {
                best = Some((i, key));
            }
        }
        let Some((k, _)) = best else { return Ok((neg, pos, step, false)) };
        let d = a[k][&k].clone();
        if d.1.is_negative() {
            neg += 1;
        } else {
            pos += 1;
        }
        let nb: Vec<(usize, Iv)> = a[k].iter().filter(|(&j, _)| j != k).map(|(&j, v)| (j, v.clone())).collect();
        for x in 0..nb.len() {
            for y in x..nb.len() {
                let (i, vi) = (&nb[x].0, &nb[x].1);
                let (j, vj) = (&nb[y].0, &nb[y].1);
                let c = quotient(vi, vj, &d)?;
                let o = a[*i].get(j).cloned().unwrap_or((BigInt::zero(), BigInt::zero()));
                let v: Iv = (o.0 - &c.1, o.1 - &c.0);
                if v.0.is_zero() && v.1.is_zero() {
                    a[*i].remove(j);
                    a[*j].remove(i);
                } else {
                    a[*i].insert(*j, v.clone());
                    a[*j].insert(*i, v);
                }
            }
        }
        let keys: Vec<usize> = a[k].keys().copied().filter(|&j| j != k).collect();
        for j in keys {
            a[j].remove(&k);
        }
        a[k].clear();
        active[k] = false;
    }
    Ok((neg, pos, n, true))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let run = || -> Result<String, String> {
        if args.len() == 2 && args[1] == "signcheck" {
            return signcheck();
        }
        if args.len() != 5 && args.len() != 6 {
            return Err("usage: sturmq2 MATRIX JE JB T [exact|interval]".into());
        }
        let p = |i: usize| args[i].parse::<i64>().map_err(|_| format!("bad arg {}", args[i]));
        let (je, jb, tn) = (p(2)?, p(3)?, p(4)?);
        let mode = args.get(5).map_or("exact", |s| s.as_str());
        let t0 = Instant::now();
        match mode {
            "exact" => {
                let input = load(&args[1], je, jb, tn)?;
                let n = input.b.len();
                let (neg, pos, zero, blocks, bits) = inertia(input.b, input.f);
                let ms = t0.elapsed().as_millis();
                Ok(format!(
                    "STURM-Q2|result=PASS|status=COMPLETE|mode=exact|field={}|JE={je}|JB={jb}|threshold={tn}/1000|states={n}|upper_entries={}|negative={neg}|positive={pos}|zero={zero}|pivot_2x2_blocks={blocks}|max_entry_bits={bits}|elapsed_ms={ms}|json=0",
                    if input.algebraic { format!("QF_{}_{}", input.f.p, input.f.q) } else { "Q".to_string() },
                    input.entries
                ))
            }
            "interval" => {
                let (mut a, entries, field) = load_interval(&args[1], je, jb, tn)?;
                let (fname, enclosure) = match field {
                    Some((1, 1)) => ("QF_1_1", "phi_enclosure=L^2-SL-S^2<=0<(L+1)^2-S(L+1)-S^2_VERIFIED"),
                    Some(_) => ("QF_2_0", "sqrt2_enclosure=L^2<=2S^2<(L+1)^2_VERIFIED"),
                    None => ("Q", "sqrt2_enclosure=L^2<=2S^2<(L+1)^2_VERIFIED"),
                };
                let n = a.len();
                let (neg, pos, done, complete) = certify(&mut a)?;
                let ms = t0.elapsed().as_millis();
                Ok(format!(
                    "STURM-Q2|result={}|status={}|mode=interval|field={}|JE={je}|JB={jb}|threshold={tn}/1000|scale_digits=80|{enclosure}|states={n}|upper_entries={entries}|pivots={done}/{n}|negative={neg}|positive={pos}|elapsed_ms={ms}|json=0",
                    if complete { "PASS" } else { "FAIL" },
                    if complete { "COMPLETE" } else { "UNCERTAIN" },
                    fname
                ))
            }
            _ => Err("mode must be exact or interval".into()),
        }
    };
    match run() {
        Ok(line) => println!("{line}"),
        Err(e) => {
            println!("STURM-Q2|result=FAIL|status=UNCERTAIN|reason={}|json=0", e.replace('|', "/"));
            std::process::exit(1);
        }
    }
}

//! fgcheck: acceptance checks G1-G3 of a finite group for the lattice gauge lane.
//! G1 group axioms (closure, associativity, identity, inverses), class count, character-table row and column
//!    orthogonality, sum of d^2 = |G|, 1-dimensional characters multiplicative, fundamental = listed irrep.
//! G2 the |G| x |G| integer matrix of E_link = sum_{g in S} (1 - L_g) on L^2(G) has characteristic polynomial
//!    prod_R (x - E_R)^(d_R^2), with E_R = |S| - sum_{g in S} chi_R(g)/d_R from the characters.
//! G3 Gauss law on the open ROWS x COLS strip: for every orbit representative, every vertex v and every h in G,
//!    the gauge-transformed raw configuration (tree links included) fixes into the same orbit, M and every
//!    Wilson loop on rows 0..1 are unchanged, and the electric move targets are the same orbits with the same
//!    multiplicities.
//! usage: fgcheck GROUP ROWS COLS        (GROUP 1 = S3, 2 = Q8)        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use fgroup::{basis, group, strip, FGroup, Sys};
use std::collections::BTreeMap;

fn pf(b: bool) -> &'static str {
    if b {
        "PASS"
    } else {
        "FAIL"
    }
}

/// characteristic polynomial det(x I - A), coefficients from x^0 up, by Faddeev-LeVerrier (exact division)
fn charpoly(a: &[Vec<i128>]) -> Option<Vec<i128>> {
    let n = a.len();
    let mut c = vec![0i128; n + 1];
    c[n] = 1;
    let mut m = vec![vec![0i128; n]; n]; // M_0 = 0
    for k in 1..=n {
        // M_k = A M_{k-1} + c_{n-k+1} I
        let mut nm = vec![vec![0i128; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0i128;
                for (t, row) in m.iter().enumerate() {
                    s = s.checked_add(a[i][t].checked_mul(row[j])?)?;
                }
                nm[i][j] = s;
            }
            nm[i][i] = nm[i][i].checked_add(c[n - k + 1])?;
        }
        // c_{n-k} = -tr(A M_k) / k
        let mut tr = 0i128;
        for i in 0..n {
            for (t, row) in nm.iter().enumerate() {
                tr = tr.checked_add(a[i][t].checked_mul(row[i])?)?;
            }
        }
        if tr % k as i128 != 0 {
            return None;
        }
        c[n - k] = -tr / k as i128;
        m = nm;
    }
    Some(c)
}

fn polymul(p: &[i128], q: &[i128]) -> Vec<i128> {
    let mut r = vec![0i128; p.len() + q.len() - 1];
    for (i, a) in p.iter().enumerate() {
        for (j, b) in q.iter().enumerate() {
            r[i + j] += a * b;
        }
    }
    r
}

fn g1(g: &FGroup) -> (bool, String) {
    let n = g.n;
    let closure = g.mul.iter().all(|r| r.iter().all(|&x| (x as usize) < n));
    let mut assoc = closure;
    if closure {
        for a in 0..n as u8 {
            for b in 0..n as u8 {
                for c in 0..n as u8 {
                    assoc &= g.m(g.m(a, b), c) == g.m(a, g.m(b, c));
                }
            }
        }
    }
    let identity = closure && (0..n as u8).all(|a| g.m(0, a) == a && g.m(a, 0) == a);
    let inverses = closure && (0..n as u8).all(|a| (a as usize) < n && g.i(a) as usize != 255 && g.m(a, g.i(a)) == 0 && g.m(g.i(a), a) == 0);
    let ok = closure && assoc && identity && inverses;
    // classes
    let mut class = vec![usize::MAX; n];
    let mut nclass = 0;
    if ok {
        for a in 0..n {
            if class[a] == usize::MAX {
                for c in 0..n as u8 {
                    class[g.m(g.m(c, a as u8), g.i(c)) as usize] = nclass;
                }
                nclass += 1;
            }
        }
    }
    let k = g.chars.len();
    let rows = (0..k).all(|r| (0..k).all(|s| (0..n).map(|x| g.chars[r][x] * g.chars[s][g.i(x as u8) as usize]).sum::<i64>() == if r == s { n as i64 } else { 0 }));
    let mut cols = ok;
    if ok {
        for a in 0..n {
            for b in 0..n {
                let s: i64 = (0..k).map(|r| g.chars[r][a] * g.chars[r][g.i(b as u8) as usize]).sum();
                let csize = (0..n).filter(|&x| class[x] == class[a]).count() as i64;
                let want = if class[a] == class[b] { n as i64 / csize } else { 0 };
                cols &= s == want;
            }
        }
    }
    let dims_ok = (0..k).all(|r| g.chars[r][0] == g.dims[r]) && g.dims.iter().map(|d| d * d).sum::<i64>() == n as i64;
    let mult = (0..k).filter(|&r| g.dims[r] == 1).all(|r| (0..n).all(|a| (0..n).all(|b| g.chars[r][g.m(a as u8, b as u8) as usize] == g.chars[r][a] * g.chars[r][b])));
    let fund = g.chars[g.fund] == g.chi_f && (0..n).all(|a| g.chi_f[a] == g.chi_f[g.i(a as u8) as usize]);
    let electric = g.electric.iter().all(|&s| g.electric.contains(&g.i(s))) && g.electric.iter().all(|&s| (0..n as u8).all(|c| g.electric.contains(&g.m(g.m(c, s), g.i(c)))));
    let all = ok && nclass == k && rows && cols && dims_ok && mult && fund && electric;
    (
        all,
        format!(
            "FG_G1|group={}|order={n}|closure={}|associativity={}|identity={}|inverses={}|classes={nclass}|irreps={k}|row_orthogonality={}|column_orthogonality={}|sum_d2={}|one_dim_multiplicative={}|fundamental_real_listed={}|electric_set_inverse_and_conjugation_closed={}|result={}|json=0",
            g.name,
            pf(closure),
            pf(assoc),
            pf(identity),
            pf(inverses),
            pf(rows),
            pf(cols),
            pf(dims_ok),
            pf(mult),
            pf(fund),
            pf(electric),
            pf(all)
        ),
    )
}

fn g2(g: &FGroup) -> (bool, String) {
    let n = g.n;
    // (E f)(x) = sum_{s in S} (f(x) - f(s x)); matrix A[x][y]
    let mut a = vec![vec![0i128; n]; n];
    for (x, row) in a.iter_mut().enumerate() {
        row[x] += g.electric.len() as i128;
        for &s in &g.electric {
            row[g.m(s, x as u8) as usize] -= 1;
        }
    }
    let cp = charpoly(&a);
    let mut want = vec![1i128];
    let mut eig: Vec<String> = Vec::new();
    let mut ok = true;
    for r in 0..g.chars.len() {
        match g.e_irrep(r) {
            Some(e) => {
                eig.push(format!("{e}^{}", g.dims[r] * g.dims[r]));
                for _ in 0..g.dims[r] * g.dims[r] {
                    want = polymul(&want, &[-i128::from(e), 1]);
                }
            }
            None => ok = false,
        }
    }
    let same = cp.as_ref() == Some(&want);
    let cps = cp.map_or("overflow".to_string(), |c| c.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(","));
    let fe = g.e_irrep(g.fund).unwrap_or(-1);
    (
        ok && same,
        format!("FG_G2|group={}|electric_set_size={}|charpoly_low_to_high={cps}|irrep_eigenvalues_with_multiplicity={}|fundamental_link_energy={fe}|charpoly_equals_character_prediction={}|json=0", g.name, g.electric.len(), eig.join(","), pf(ok && same)),
    )
}

fn g3(g: &FGroup, rows: usize, cols: usize) -> (bool, String) {
    let s = strip(rows, cols);
    let sys = Sys { g, s: &s };
    let links: Vec<usize> = (0..s.edges.len()).collect();
    let b = basis(&sys, true, &links);
    let (mut fix_ok, mut mag_ok, mut wil_ok, mut move_ok) = (true, true, true, true);
    let mut checks = 0u64;
    let targets = |l: &[u8]| -> BTreeMap<u64, i64> {
        let mut m = BTreeMap::new();
        for e in 0..s.edges.len() {
            for &t in &g.electric {
                let mut u = l.to_vec();
                u[e] = g.m(t, u[e]);
                *m.entry(sys.key(&sys.fix(&u))).or_insert(0) += 1;
            }
        }
        m
    };
    for q in &b.reps {
        let raw = sys.raw(q);
        let k0 = sys.key(q);
        let m0 = sys.magnetic(&raw);
        let w0: Vec<i64> = (1..cols).map(|x1| sys.wilson(&raw, 0, x1)).collect();
        let t0 = targets(&raw);
        for v in 0..s.sites {
            for h in 0..g.n as u8 {
                let gl = sys.gauge(&raw, v, h);
                fix_ok &= sys.key(&sys.fix(&gl)) == k0;
                mag_ok &= sys.magnetic(&gl) == m0;
                wil_ok &= (1..cols).map(|x1| sys.wilson(&gl, 0, x1)).collect::<Vec<i64>>() == w0;
                move_ok &= targets(&gl) == t0;
                checks += 1;
            }
        }
    }
    let all = fix_ok && mag_ok && wil_ok && move_ok;
    (
        all,
        format!(
            "FG_G3|group={}|rows={rows}|cols={cols}|orbit_states={}|gauge_checks={checks}|same_orbit={}|magnetic_invariant={}|wilson_invariant={}|electric_moves_equivariant={}|result={}|json=0",
            g.name,
            b.reps.len(),
            pf(fix_ok),
            pf(mag_ok),
            pf(wil_ok),
            pf(move_ok),
            pf(all)
        ),
    )
}

fn main() {
    let a: Vec<i64> = std::env::args().skip(1).filter_map(|x| x.parse().ok()).collect();
    let (Some(g), true) = (a.first().and_then(|&id| group(id)), a.len() == 3 && (2..=3).contains(&a[1]) && (2..=4).contains(&a[2])) else {
        eprintln!("usage: fgcheck GROUP ROWS COLS");
        std::process::exit(1);
    };
    let (o1, l1) = g1(&g);
    println!("{l1}");
    let (o2, l2) = g2(&g);
    println!("{l2}");
    let (o3, l3) = if o1 { g3(&g, a[1] as usize, a[2] as usize) } else { (false, "FG_G3|result=FAIL|reason=group_invalid|json=0".to_string()) };
    println!("{l3}");
    println!("FG_CHECK|group={}|G1={}|G2={}|G3={}|json=0", g.name, pf(o1), pf(o2), pf(o3));
    std::process::exit(if o1 && o2 && o3 { 0 } else { 2 });
}

//! fgconf: finite-group gauge Hamiltonian on an open ROWS x COLS strip in the exchange format read by sturmrs
//! (weights, magnetic diagonal, electric triples), vacuum sector (CH = 0) or Wilson mode (CH = 2: diagonal
//! JB M + LAM chi_f(U_C) for the rectangle on rows 0..1, columns XA..XB). For S3 the output is byte-identical
//! to ymconf's. When every diagonal value is a rational integer the standard integer format is written (read
//! by sturmrs and sturmq2); when the group's field enters the diagonal (2O: sqrt2), the file starts with the
//! line "QF P Q" (theta^2 = P + Q theta) and each diagonal value is written as "X:Y" for X + Y theta (read by
//! sturmq2). The electric part is always an integer matrix.
//! The written operator is checked exactly W-self-adjoint (w_i Z_ij = w_j Z_ji).
//! usage: fgconf GROUP ROWS COLS XA XB CH OUT_MATRIX [JB LAM]        (GROUP 1 = S3, 2 = Q8, 3 = 2T, 4 = 2O)        json=0
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use fgroup::{basis, group, strip, Alg, Sys};
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let p: Vec<i64> = args.iter().skip(1).take(6).filter_map(|x| x.parse().ok()).collect();
    let (true, Some(out)) = (p.len() == 6, args.get(7)) else {
        eprintln!("usage: fgconf GROUP ROWS COLS XA XB CH OUT_MATRIX [JB LAM]");
        std::process::exit(1);
    };
    let wl: Option<(i64, i64)> = match (args.get(8).and_then(|v| v.parse().ok()), args.get(9).and_then(|v| v.parse().ok())) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    };
    let (rows, cols, xa, xb, ch) = (p[1], p[2], p[3], p[4], p[5]);
    let geo = (2..=3).contains(&rows) && (2..=5).contains(&cols);
    let wil = ch == 2 && 0 <= xa && xa < xb && xb < cols && wl.is_some_and(|(jb, _)| jb >= 0);
    let vac = ch == 0 && xa == 0 && xb == 0;
    let Some(g) = group(p[0]).filter(|_| geo && (wil || vac)) else {
        eprintln!("FG_MATRIX|result=FAIL|reason=out_of_range|json=0");
        std::process::exit(1);
    };
    let (rows, cols, xa, xb) = (rows as usize, cols as usize, xa as usize, xb as usize);
    let s = strip(rows, cols);
    let sys = Sys { g: &g, s: &s };
    let links: Vec<usize> = (0..s.edges.len()).collect();
    let b = basis(&sys, true, &links);
    let n = b.reps.len();
    let mags: Vec<Alg> = b
        .reps
        .iter()
        .map(|q| {
            let raw = sys.raw(q);
            match wl.filter(|_| ch == 2) {
                Some((jb, lam)) => sys.magnetic(&raw).scale(jb).plus(sys.wilson(&raw, xa, xb).scale(lam)),
                None => sys.magnetic(&raw),
            }
        })
        .collect();
    let rational = mags.iter().all(|m| m.y == 0);
    let mut adjoint = true;
    for i in 0..n {
        for (&j, &v) in &b.e[i] {
            adjoint &= b.w[i] * v == b.w[j] * b.e[j].get(&i).copied().unwrap_or(0);
        }
    }
    let join = |v: &[i64]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
    let mag_text = if rational {
        mags.iter().map(|m| m.x.to_string()).collect::<Vec<_>>().join(" ")
    } else {
        mags.iter().map(|m| format!("{}:{}", m.x, m.y)).collect::<Vec<_>>().join(" ")
    };
    let header = if rational { String::new() } else { format!("QF {} {}\n", g.field.p, g.field.q) };
    let mut text = format!("{header}{n}\n{}\n{mag_text}\n", join(&b.w));
    let mut entries = 0usize;
    for (i, row) in b.e.iter().enumerate() {
        for (&j, &v) in row.range(i..) {
            text.push_str(&format!("{i} {j} {v}\n"));
            entries += 1;
        }
    }
    let ok = adjoint && std::fs::File::create(out).and_then(|mut f| f.write_all(text.as_bytes())).is_ok();
    println!(
        "FG_MATRIX|group={}|rows={rows}|cols={cols}|xa={xa}|xb={xb}|mode={}|{}links={}|plaquettes={}|free_links={}|states={n}|upper_entries={entries}|{}w_self_adjoint={}|json=0",
        g.name,
        if ch == 2 { "wilson" } else { "vacuum" },
        wl.filter(|_| ch == 2).map_or(String::new(), |(jb, lam)| format!("wilson_loop=1x{}|diag=JB*M+LAM*W|JB={jb}|LAM={lam}|", xb - xa)),
        s.edges.len(),
        s.plaq.len(),
        s.free.len(),
        if rational { String::new() } else { format!("diagonal_field=QF_{}_{}|", g.field.p, g.field.q) },
        if adjoint { "PASS" } else { "FAIL" }
    );
    std::process::exit(if ok { 0 } else { 2 });
}

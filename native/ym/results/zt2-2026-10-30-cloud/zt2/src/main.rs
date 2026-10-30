//! zt2 driver: binary tetrahedral group 2T, second SU(2)-subgroup milestone: open 2 x 3 strip (series in two
//! independent algorithms and bases, certified loops) and open 2 x 4 ladder (series in two bases). Every
//! computation is a job of the HEB-C HTTP-0 backend (hebc_ym_http0 serve) over HTTP/1.1. See PREREG.md.
//! Rust 1.81, integer arithmetic only (i128 and big integers), json=0.
//! usage: zt2 PORT
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::net::TcpStream;

const JE: i128 = 1_350_000_000_000_000;
const UNIT: i128 = 1000 * JE;
const LAM: i128 = 13_500_000; // lambda = 10^-8
const M_FLOOR: i128 = 1_350_000_000_000; // 10^-6 J_E
const NP: i128 = 2; // the 2 x 3 strip carries the certified values
const COLS: usize = 3;
const RS: [(i128, i128); 4] = [(1, 27), (1, 9), (1, 3), (1, 1)];
const LOOPS: [usize; 2] = [1, 2];
const ORDER: usize = 5;
const G: i64 = 3; // 2T
const ZWILSON_RESULT: &str = "../zwilson-2026-10-30-cloud/RESULT.hbp";
const ZWILSON_RESULT_SHA: &str = "8fe1c2a4d8fa0a7058d290f0ddc020a22b0422ca9bb888db1f32a6e1a794c818";
const ZQ8_RESULT: &str = "../zq8-2026-10-30-cloud/RESULT.hbp";
const ZQ8_RESULT_SHA: &str = "b6c129eb5a9982e5014622fb0fead32e0012372b9a82b3af166c88cd07f1d50b";
const ZQ8_MATRIX_SPEC: &str = "fgwilsonmat/2/2/4/0/1/40000000000000/0";
const ZQ8_MATRIX_HBP: &str = "2b79db3ff77fc1f491d1d4ad661a0b93639b25efd9715c446ded5137301bf038";

fn http(port: u16, method: &str, path: &str) -> Result<(u16, String), String> {
    let mut s = TcpStream::connect(("127.0.0.1", port)).map_err(|e| e.to_string())?;
    s.write_all(format!("{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Length: 0\r\n\r\n").as_bytes()).map_err(|e| e.to_string())?;
    let mut buf = String::new();
    s.read_to_string(&mut buf).map_err(|e| e.to_string())?;
    let (head, body) = buf.split_once("\r\n\r\n").ok_or("malformed response")?;
    let code = head.split_ascii_whitespace().nth(1).and_then(|c| c.parse().ok()).ok_or("no status")?;
    Ok((code, body.to_string()))
}

fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split('|').find_map(|f| f.strip_prefix(key))
}

#[derive(Clone, Debug)]
struct Q {
    n: BigInt,
    d: BigInt,
}

fn q(n: i128, d: i128) -> Q {
    Q { n: BigInt::from(n), d: BigInt::from(d) }
}

impl Q {
    fn parse(s: &str) -> Option<Q> {
        let (a, b) = s.split_once('/')?;
        Some(Q { n: a.parse().ok()?, d: b.parse().ok()? })
    }
    fn mul(&self, o: &Q) -> Q {
        Q { n: &self.n * &o.n, d: &self.d * &o.d }
    }
    fn div(&self, o: &Q) -> Q {
        let (n, d) = (&self.n * &o.d, &self.d * &o.n);
        if d.is_negative() {
            Q { n: -n, d: -d }
        } else {
            Q { n, d }
        }
    }
    fn sub(&self, o: &Q) -> Q {
        Q { n: &self.n * &o.d - &o.n * &self.d, d: &self.d * &o.d }
    }
    fn add(&self, o: &Q) -> Q {
        Q { n: &self.n * &o.d + &o.n * &self.d, d: &self.d * &o.d }
    }
    fn lt(&self, o: &Q) -> bool {
        &self.n * &o.d < &o.n * &self.d
    }
    fn eq(&self, o: &Q) -> bool {
        &self.n * &o.d == &o.n * &self.d
    }
    fn abs(&self) -> Q {
        Q { n: self.n.abs(), d: self.d.clone() }
    }
    fn pow(&self, k: usize) -> Q {
        (0..k).fold(q(1, 1), |a, _| a.mul(self))
    }
    fn dec(&self) -> String {
        let s = BigInt::from(10).pow(15);
        let v = num_integer::Integer::div_floor(&(&self.n * &s), &self.d);
        let neg = v.is_negative();
        let a = v.abs();
        let (ip, fp) = (num_integer::Integer::div_floor(&a, &s), num_integer::Integer::mod_floor(&a, &s));
        format!("{}{}.{:0>15}", if neg { "-" } else { "" }, ip, fp.to_string())
    }
}

/// SHA-256 (FIPS 180-4), for checking the sealed zwilson result this package compares against
fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let mut m = data.to_vec();
    let bits = (data.len() as u64).wrapping_mul(8);
    m.push(0x80);
    while m.len() % 64 != 56 {
        m.push(0);
    }
    m.extend_from_slice(&bits.to_be_bytes());
    for chunk in m.chunks_exact(64) {
        let mut w = [0u32; 64];
        for t in 0..16 {
            w[t] = u32::from_be_bytes([chunk[4 * t], chunk[4 * t + 1], chunk[4 * t + 2], chunk[4 * t + 3]]);
        }
        for t in 16..64 {
            let s0 = w[t - 15].rotate_right(7) ^ w[t - 15].rotate_right(18) ^ (w[t - 15] >> 3);
            let s1 = w[t - 2].rotate_right(17) ^ w[t - 2].rotate_right(19) ^ (w[t - 2] >> 10);
            w[t] = w[t - 16].wrapping_add(s0).wrapping_add(w[t - 7]).wrapping_add(s1);
        }
        let mut v = h;
        for t in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[t]).wrapping_add(w[t]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for (a, b) in h.iter_mut().zip(v.iter()) {
            *a = a.wrapping_add(*b);
        }
    }
    h.iter().fold(String::new(), |mut acc, x| {
        use std::fmt::Write as _;
        let _ = write!(acc, "{x:08x}");
        acc
    })
}

struct Client {
    port: u16,
    out: Vec<String>,
    cache: BTreeMap<String, Option<i64>>,
}

impl Client {
    fn job(&mut self, spec: &str) -> Result<(String, String), String> {
        let (code, body) = http(self.port, "POST", &format!("/job/{spec}"))?;
        if code != 200 {
            return Err(body.trim().replace('|', "/"));
        }
        let hs = field(body.trim(), "hbp_sha256=").unwrap_or("").to_string();
        let (_, shard) = http(self.port, "GET", &format!("/shard/{spec}"))?;
        Ok((hs, shard))
    }
    fn verify(&mut self, spec: &str) -> bool {
        let ok = matches!(http(self.port, "GET", &format!("/verify/{spec}")), Ok((200, ref b)) if b.contains("result=PASS|action=verify"));
        self.out.push(format!("VERIFY|spec={spec}|result={}|json=0", if ok { "PASS" } else { "FAIL" }));
        ok
    }
    fn matrix(&mut self, x1: usize, jb: i128, lam: i128) -> bool {
        let spec = format!("fgwilsonmat/{G}/2/{COLS}/0/{x1}/{jb}/{lam}");
        if let Some(v) = self.cache.get(&spec) {
            return v.is_some();
        }
        let r = self.job(&spec);
        let ok = matches!(&r, Ok((_, s)) if s.contains("w_self_adjoint=PASS"));
        match &r {
            Ok((hs, shard)) => {
                let fm = shard.lines().find(|l| l.starts_with("FG_MATRIX|")).unwrap_or("FG_MATRIX|missing");
                self.out.push(format!("MATRIX|spec={spec}|ok={}|states={}|w_self_adjoint={}|hbp_sha256={hs}|json=0", u8::from(ok), field(fm, "states=").unwrap_or("-"), field(fm, "w_self_adjoint=").unwrap_or("-")));
            }
            Err(e) => self.out.push(format!("MATRIX|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0")),
        }
        self.cache.insert(spec, if ok { Some(1) } else { None });
        ok
    }
    fn count(&mut self, x1: usize, jb: i128, lam: i128, t: i128) -> Option<i64> {
        let spec = format!("fgwilsoncount/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{t}");
        if let Some(v) = self.cache.get(&spec) {
            return *v;
        }
        let v = match self.job(&spec) {
            Ok((hs, shard)) => {
                let line = shard.lines().find(|l| l.starts_with("STURM-RS|")).unwrap_or("");
                match (field(line, "status="), field(line, "negative=").and_then(|v| v.parse::<i64>().ok())) {
                    (Some("COMPLETE"), Some(n)) => {
                        self.out.push(format!("COUNT|spec={spec}|negative={n}|hbp_sha256={hs}|json=0"));
                        Some(n)
                    }
                    _ => {
                        self.out.push(format!("COUNT|spec={spec}|result=INCONCLUSIVE|reason=uncertain|json=0"));
                        None
                    }
                }
            }
            Err(e) => {
                self.out.push(format!("COUNT|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0"));
                None
            }
        };
        self.cache.insert(spec, v);
        v
    }
    fn bisect(&mut self, x1: usize, jb: i128, lam: i128, mut lo: i128, mut hi: i128) -> Option<(i128, i128)> {
        while hi - lo > 1 {
            let mid = lo + (hi - lo) / 2;
            if self.count(x1, jb, lam, mid)? == 0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let a = self.verify(&format!("fgwilsoncount/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{lo}"));
        let b = self.verify(&format!("fgwilsoncount/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{hi}"));
        (a && b).then_some((lo, hi))
    }
}

fn series_of(shard: &str) -> BTreeMap<usize, Vec<Q>> {
    let mut m = BTreeMap::new();
    for l in shard.lines() {
        if let (Some(lp), Some(cs)) = (field(l, "loop=1x"), field(l, "c_k=")) {
            let v: Option<Vec<Q>> = cs.split(',').map(Q::parse).collect();
            if let (Ok(a), Some(v)) = (lp.parse::<usize>(), v) {
                m.insert(a, v);
            }
        }
    }
    m
}

/// one exact series job: (hbp sha, verify, basis lines ok, loop -> c_k, energies)
type SeriesOut = (bool, BTreeMap<usize, Vec<Q>>, Vec<Q>);

fn series_job(c: &mut Client, spec: &str) -> SeriesOut {
    match c.job(spec) {
        Ok((hs, shard)) => {
            for l in shard.lines().filter(|l| l.starts_with("FGSERIES")) {
                c.out.push(l.to_string());
            }
            let basis_ok = shard.lines().any(|l| l.starts_with("FGSERIES_BASIS|") && l.contains("e_w_self_adjoint=PASS") && l.contains("rhs_orthogonal_and_solved=PASS") && l.contains("orthogonality=PASS"));
            let e: Vec<Q> = shard.lines().find_map(|l| field(l, "E_k=")).map(|e| e.split(',').filter_map(Q::parse).collect()).unwrap_or_default();
            let ver = c.verify(spec);
            c.out.push(format!("SERIES_JOB|spec={spec}|hbp_sha256={hs}|verify={}|json=0", if ver { "PASS" } else { "FAIL" }));
            (basis_ok && ver, series_of(&shard), e)
        }
        Err(e) => {
            c.out.push(format!("SERIES_JOB|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0"));
            (false, BTreeMap::new(), Vec::new())
        }
    }
}

fn sealed_series(path: &str, sha: &str, prefix: &str, method: &str) -> (bool, Vec<String>) {
    let bytes = std::fs::read(path).unwrap_or_default();
    let ok = sha256_hex(&bytes) == sha;
    let text = String::from_utf8(bytes).unwrap_or_default();
    let v = text.lines().filter(|l| l.starts_with(prefix) && l.contains(method)).filter_map(|l| field(l, "c_k=").map(str::to_string)).collect();
    (ok, v)
}

fn main() {
    let port: u16 = std::env::args().nth(1).and_then(|p| p.parse().ok()).unwrap_or(0);
    let mut c = Client { port, out: Vec::new(), cache: BTreeMap::new() };
    c.out.push(format!("ZT2|group=2T|electric_set=(1+-i+-j+-k)/2|certified_geometry=S2_open_2x3|series_geometries=S2_open_2x3,S3_open_2x4|JE={JE}|lambda=LAM/JE={LAM}/{JE}|r=1/27,1/9,1/3,1|json=0"));

    // G1-G3 on S2
    let mut g123 = (false, false, false);
    let check = format!("fgcheck/{G}/2/{COLS}");
    match c.job(&check) {
        Ok((hs, shard)) => {
            for l in shard.lines().filter(|l| l.starts_with("FG_")) {
                c.out.push(l.to_string());
            }
            let fc = shard.lines().find(|l| l.starts_with("FG_CHECK|")).unwrap_or("");
            let g2line = shard.lines().find(|l| l.starts_with("FG_G2|")).unwrap_or("");
            let reg = field(g2line, "charpoly_low_to_high=") == Some(&registered_charpoly());
            let ver = c.verify(&check);
            g123 = (field(fc, "G1=") == Some("PASS") && ver, field(fc, "G2=") == Some("PASS") && reg && ver, field(fc, "G3=") == Some("PASS") && ver);
            c.out.push(format!("CHECK_JOB|spec={check}|hbp_sha256={hs}|charpoly_is_x(x-4)^4(x-8)^9(x-10)^8(x-12)^2={}|verify={}|json=0", u8::from(reg), if ver { "PASS" } else { "FAIL" }));
        }
        Err(e) => c.out.push(format!("CHECK_JOB|result=INCONCLUSIVE|reason={e}|json=0")),
    }

    // R1: S3 and Q8 regressions against the sealed packages
    let mut r1 = true;
    for (group, order, path, sha, prefix, ma, mb) in [
        (1, 4, ZWILSON_RESULT, ZWILSON_RESULT_SHA, "WSERIES|loop=", "method=A_orbit", "method=B_config"),
        (2, 5, ZQ8_RESULT, ZQ8_RESULT_SHA, "FGSERIES|group=Q8|loop=", "method=A_orbit_elimination", "method=B_config_spectral"),
    ] {
        for (method, tag) in [(0, ma), (1, mb)] {
            let (sealed_ok, want) = sealed_series(path, sha, prefix, tag);
            let spec = format!("fgwilsonseries/{group}/2/4/{method}/{order}");
            match c.job(&spec) {
                Ok((hs, shard)) => {
                    let got: Vec<String> = shard.lines().filter(|l| l.starts_with("FGSERIES|")).filter_map(|l| field(l, "c_k=").map(str::to_string)).collect();
                    let same = sealed_ok && got == want && got.len() == 3;
                    let ver = c.verify(&spec);
                    r1 &= same && ver;
                    c.out.push(format!("REGRESSION_SERIES|spec={spec}|identical_to_sealed={}|sealed_result_sha_ok={}|hbp_sha256={hs}|verify={}|json=0", u8::from(same), u8::from(sealed_ok), if ver { "PASS" } else { "FAIL" }));
                }
                Err(e) => {
                    r1 = false;
                    c.out.push(format!("REGRESSION_SERIES|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0"));
                }
            }
        }
    }
    {
        let (jb, lam) = (40_000_000_000_000i128, 10_800_000i128);
        let old = c.job(&format!("wilsonmat/2/4/0/2/{jb}/{lam}"));
        let new = c.job(&format!("fgwilsonmat/1/2/4/0/2/{jb}/{lam}"));
        let sha = |r: &Result<(String, String), String>| -> String {
            r.as_ref().ok().and_then(|(_, s)| s.lines().find_map(|l| l.strip_prefix("MATRIX_SHA|sha256=").map(|v| v[..64.min(v.len())].to_string()))).unwrap_or_default()
        };
        let (a, b) = (sha(&old), sha(&new));
        let same = !a.is_empty() && a == b;
        r1 &= same;
        c.out.push(format!("REGRESSION_S3_MATRIX|wilsonmat_sha256={a}|fgwilsonmat_sha256={b}|byte_identical={}|json=0", u8::from(same)));
        let q8m = c.job(ZQ8_MATRIX_SPEC);
        let hs = q8m.as_ref().map(|(h, _)| h.clone()).unwrap_or_default();
        let same8 = hs == ZQ8_MATRIX_HBP;
        r1 &= same8;
        c.out.push(format!("REGRESSION_Q8_MATRIX|spec={ZQ8_MATRIX_SPEC}|hbp_sha256={hs}|sealed_zq8_hbp_sha256={ZQ8_MATRIX_HBP}|identical={}|json=0", u8::from(same8)));
    }

    // W1: exact 2T series, S2 (elimination vs configuration spectral) and S3 (orbit spectral vs configuration spectral)
    let predicted_lead = [q(1, 8), q(7, 768), q(11, 16384)];
    let mut w1 = true;
    let mut lead: BTreeMap<usize, Q> = BTreeMap::new();
    let mut s2full: BTreeMap<usize, Vec<Q>> = BTreeMap::new();
    for (cols, ma, np) in [(3usize, 0usize, 2i128), (4, 2, 3)] {
        let a = series_job(&mut c, &format!("fgwilsonseries/{G}/2/{cols}/{ma}/{ORDER}"));
        let b = series_job(&mut c, &format!("fgwilsonseries/{G}/2/{cols}/1/{ORDER}"));
        w1 &= a.0 && b.0;
        for x in 1..cols {
            let (Some(sa), Some(sb)) = (a.1.get(&x), b.1.get(&x)) else {
                w1 = false;
                continue;
            };
            let same = sa.len() == ORDER + 1 && sb.len() == ORDER + 1 && sa.iter().zip(sb).all(|(u, v)| u.eq(v));
            let below = (0..x).all(|k| sa[k].n.is_zero());
            let lead_ok = sa[x].eq(&predicted_lead[x - 1]);
            let next_zero = sa[x + 1].n.is_zero();
            w1 &= same && below && lead_ok && next_zero;
            if cols == COLS {
                lead.insert(x, sa[x].clone());
                s2full.insert(x, sa.clone());
            }
            c.out.push(format!(
                "SERIES_CHECK|geometry=2x{cols}|loop=1x{x}|methods_identical={}|zero_below_order_A={}|c_AA={}/{}|c_AA_equals_(3/4)^A_S3={}|c_A_A+1={}/{}|centre_symmetry_zero={}|normalised_W/2_c_AA={}/{}|json=0",
                u8::from(same),
                u8::from(below),
                sa[x].n,
                sa[x].d,
                u8::from(lead_ok),
                sa[x + 1].n,
                sa[x + 1].d,
                u8::from(next_zero),
                sa[x].n,
                &sa[x].d * 2
            ));
        }
        let e_ok = a.2.len() == ORDER && b.2.len() == ORDER && a.2.iter().zip(&b.2).all(|(u, v)| u.eq(v)) && a.2[0].eq(&q(2 * np, 1)) && a.2[1].eq(&q(-np, 16)) && a.2[2].n.is_zero() && a.2[4].n.is_zero();
        w1 &= e_ok;
        c.out.push(format!("SERIES_ENERGY_CHECK|geometry=2x{cols}|E1={}|E2=-{np}/16|E3=0|E5=0|methods_identical_and_registered={}|json=0", 2 * np, u8::from(e_ok)));
    }

    // certified values on S2
    let mut wv: BTreeMap<(usize, usize), (Q, Q)> = BTreeMap::new();
    let mut w2 = true;
    let mut inconclusive: Vec<String> = Vec::new();
    for (ri, &(rn, rd)) in RS.iter().enumerate() {
        let jb = JE * rn / rd;
        let second = -1000 * NP * jb * jb / (16 * JE);
        let pred = 2000 * NP * jb + second;
        let mut m = (10 * second.abs()).max(M_FLOOR);
        let e0 = (|| -> Option<(i128, i128)> {
            if !c.matrix(1, jb, 0) {
                return None;
            }
            for _ in 0..=20 {
                let (lo, hi) = (pred - m, pred + m);
                if c.count(1, jb, 0, lo)? == 0 && c.count(1, jb, 0, hi)? >= 1 {
                    return c.bisect(1, jb, 0, lo, hi);
                }
                m *= 2;
            }
            None
        })();
        let Some((l0, u0)) = e0 else {
            w2 = false;
            inconclusive.push(format!("E0:r={rn}/{rd}"));
            continue;
        };
        c.out.push(format!("ENERGY0|r={rn}/{rd}|T_lo={l0}|T_hi={u0}|E0_lo={}|json=0", q(l0, UNIT).dec()));
        for &a in &LOOPS {
            let mut ends = Vec::new();
            for lam in [LAM, -LAM] {
                let got = (|| -> Option<(i128, i128)> {
                    if !c.matrix(a, jb, lam) {
                        return None;
                    }
                    let (lo, hi) = (l0 - 2000 * LAM - 1, u0 + 2000 * LAM + 1);
                    if c.count(a, jb, lam, lo)? == 0 && c.count(a, jb, lam, hi)? >= 1 {
                        c.bisect(a, jb, lam, lo, hi)
                    } else {
                        None
                    }
                })();
                ends.push(got);
            }
            let (Some((lp, _)), Some((lm, _))) = (ends[0], ends[1]) else {
                w2 = false;
                inconclusive.push(format!("1x{a}:r={rn}/{rd}"));
                continue;
            };
            let (lo, hi) = (q(lp - u0, 1000 * LAM), q(u0 - lm, 1000 * LAM));
            let width_t = u0 - lm - (lp - u0);
            // width <= 2.5e-9: width_T / (1000 LAM) <= 25 / 10^10
            let width_ok = width_t * 10_000_000_000 <= 25 * 1000 * LAM;
            let pos = lp - u0 > 0;
            w2 &= width_ok && pos;
            c.out.push(format!(
                "WILSON|loop=1x{a}|r={rn}/{rd}|W_lo={}/{}|W_hi={}/{}|W_lo_dec={}|W_hi_dec={}|normalised_W/2_lo_dec={}|width_T={width_t}|width_le_2.5e-9={}|positive={}|json=0",
                lo.n,
                lo.d,
                hi.n,
                hi.d,
                lo.dec(),
                hi.dec(),
                lo.div(&q(2, 1)).dec(),
                u8::from(width_ok),
                u8::from(pos)
            ));
            wv.insert((a, ri), (lo, hi));
        }
    }

    // W3, W4 (with the resolution rule), series comparison
    let rq = |ri: usize| q(RS[ri].0, RS[ri].1);
    let one = q(1, 1);
    let mut w3 = true;
    let (mut w4_resolved, mut w4_fail) = (0usize, false);
    for &a in &LOOPS {
        let (Some(ca), Some(full)) = (lead.get(&a), s2full.get(&a)) else {
            w3 = false;
            w4_fail = true;
            continue;
        };
        let ratio = |ri: usize| -> Option<(Q, Q)> {
            let (lo, hi) = wv.get(&(a, ri))?;
            let p = ca.mul(&rq(ri).pow(a));
            Some((lo.div(&p), hi.div(&p)))
        };
        let (Some((r27l, r27h)), Some((r9l, r9h))) = (ratio(0), ratio(1)) else {
            w3 = false;
            w4_fail = true;
            continue;
        };
        let ok3 = !r27l.lt(&q(4, 5)) && !q(5, 4).lt(&r27h);
        w3 &= ok3;
        let delta = full[a + 2].div(ca).abs().mul(&q(1, 81));
        let resolved = r27h.sub(&r27l).lt(&delta) && r9h.sub(&r9l).lt(&delta);
        let maxd = |l: &Q, h: &Q| -> Q {
            let (x, y) = (l.sub(&one).abs(), h.sub(&one).abs());
            if x.lt(&y) {
                y
            } else {
                x
            }
        };
        let mind = |l: &Q, h: &Q| -> Q {
            if !one.lt(l) && !h.lt(&one) {
                q(0, 1)
            } else {
                let (x, y) = (l.sub(&one).abs(), h.sub(&one).abs());
                if x.lt(&y) {
                    x
                } else {
                    y
                }
            }
        };
        let conv = maxd(&r27l, &r27h).lt(&mind(&r9l, &r9h));
        if resolved {
            w4_resolved += 1;
            w4_fail |= !conv;
        }
        c.out.push(format!(
            "LEADING|loop=1x{a}|ratio_r=1/27=[{},{}]|ratio_r=1/9=[{},{}]|within_0.8_1.25={}|series_delta_at_1/9={}|resolved={}|converging={}|w4={}|json=0",
            r27l.dec(),
            r27h.dec(),
            r9l.dec(),
            r9h.dec(),
            u8::from(ok3),
            delta.dec(),
            u8::from(resolved),
            u8::from(conv),
            if !resolved {
                "INCONCLUSIVE_BY_RESOLUTION"
            } else if conv {
                "PASS"
            } else {
                "FAIL"
            }
        ));
        for (ri, &(rn, rd)) in RS.iter().enumerate() {
            let r = rq(ri);
            let t1 = ca.mul(&r.pow(a));
            let t2 = t1.add(&full[a + 1].mul(&r.pow(a + 1))).add(&full[a + 2].mul(&r.pow(a + 2)));
            if let Some((lo, hi)) = wv.get(&(a, ri)) {
                let inside = !t2.lt(lo) && !hi.lt(&t2);
                c.out.push(format!("SERIES_COMPARE|loop=1x{a}|r={rn}/{rd}|leading={}|through_order_A+2={}|inside={}|json=0", t1.dec(), t2.dec(), u8::from(inside)));
            }
        }
    }
    // W5 area law (1x2 / 1x1)
    let mut w5 = true;
    for (ri, &(rn, rd)) in RS.iter().enumerate() {
        let (Some((al, ah)), Some((bl, bh))) = (wv.get(&(1, ri)), wv.get(&(2, ri))) else {
            w5 = false;
            continue;
        };
        let dec = bh.lt(al);
        w5 &= dec;
        let mut line = format!("AREA|loops=1x2/1x1|r={rn}/{rd}|ratio=[{},{}]|decreasing={}", bl.div(ah).dec(), bh.div(al).dec(), u8::from(dec));
        if ri == 0 {
            if let (Some(ca), Some(cb)) = (lead.get(&1), lead.get(&2)) {
                let pred = cb.div(ca).mul(&rq(0));
                let (lo, hi) = (bl.div(ah).div(&pred), bh.div(al).div(&pred));
                let ok = !lo.lt(&q(4, 5)) && !q(5, 4).lt(&hi);
                w5 &= ok;
                line.push_str(&format!("|predicted_ratio={}|within_0.8_1.25={}", pred.dec(), u8::from(ok)));
            } else {
                w5 = false;
            }
        }
        c.out.push(format!("{line}|json=0"));
    }
    // W6 even correction of the 1x1 loop
    let w6 = match (wv.get(&(1, 0)), wv.get(&(1, 1))) {
        (Some((l27, h27)), Some((l9, h9))) => {
            let d = |l: &Q, h: &Q, ri: usize| -> (Q, Q) {
                let p = rq(ri).div(&q(8, 1));
                (l.div(&p).sub(&one), h.div(&p).sub(&one))
            };
            let (d27l, d27h) = d(l27, h27, 0);
            let (d9l, d9h) = d(l9, h9, 1);
            let pos = q(0, 1).lt(&d27l) && q(0, 1).lt(&d9l);
            let neg = d27h.lt(&q(0, 1)) && d9h.lt(&q(0, 1));
            let (rlo, rhi) = if pos {
                (d9l.div(&d27h), d9h.div(&d27l))
            } else if neg {
                (d9h.div(&d27l), d9l.div(&d27h))
            } else {
                (q(0, 1), q(0, 1))
            };
            let inside = |a: i128, b: i128| !rlo.lt(&q(a, 1)) && !q(b, 1).lt(&rhi);
            let ok = (pos || neg) && (inside(6, 12) || inside(54, 108));
            c.out.push(format!("EVEN_CORRECTION|loop=1x1|D_1/27=[{},{}]|D_1/9=[{},{}]|ratio=[{},{}]|in_6_12_or_54_108={}|json=0", d27l.dec(), d27h.dec(), d9l.dec(), d9h.dec(), rlo.dec(), rhi.dec(), u8::from(ok)));
            ok
        }
        _ => false,
    };
    let inc = !inconclusive.is_empty();
    let v = |b: bool| if inc { "INCONCLUSIVE" } else if b { "PASS" } else { "FAIL" };
    let p = |b: bool| if b { "PASS" } else { "FAIL" };
    let w4 = if inc {
        "INCONCLUSIVE"
    } else if w4_fail {
        "FAIL"
    } else if w4_resolved == 0 {
        "INCONCLUSIVE"
    } else {
        "PASS"
    };
    c.out.push(format!(
        "ZT2_RESULT|G1={}|G2={}|G3={}|G4={}|R1={}|W1={}|W2={}|W3={}|W4={w4}|W4_resolved_loops={w4_resolved}|W5={}|W6={}|inconclusive={}|json=0",
        p(g123.0),
        p(g123.1),
        p(g123.2),
        p(c.out.iter().filter(|l| l.starts_with("MATRIX|")).all(|l| l.contains("|w_self_adjoint=PASS|")) && c.out.iter().any(|l| l.starts_with("MATRIX|"))),
        p(r1),
        p(w1),
        v(w2),
        v(w3),
        v(w5),
        v(w6),
        if inc { inconclusive.join(",") } else { "none".into() }
    ));
    for l in &c.out {
        println!("{l}");
    }
}

/// x (x - 4)^4 (x - 8)^9 (x - 10)^8 (x - 12)^2, coefficients low to high, as written by fgcheck
fn registered_charpoly() -> String {
    let mut p: Vec<i128> = vec![1];
    let mut mulby = |r: i128| {
        let mut n = vec![0i128; p.len() + 1];
        for (i, c) in p.iter().enumerate() {
            n[i + 1] += c;
            n[i] -= c * r;
        }
        p = n;
    };
    mulby(0);
    for (root, mult) in [(4, 4), (8, 9), (10, 8), (12, 2)] {
        for _ in 0..mult {
            mulby(root);
        }
    }
    p.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",")
}

#[cfg(test)]
mod tests {
    #[test]
    fn sha256_known_vectors() {
        assert_eq!(super::sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(super::sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
        let long = vec![b'a'; 1000];
        assert_eq!(super::sha256_hex(&long), "41edece42d63e8d9bf515a9ba6932e1c20cbc9f5a5d134645adb5db1b9737ea3");
    }
}

//! z2o driver: binary octahedral group 2O on the open 2 x 3 strip, third SU(2)-subgroup milestone. Characters,
//! the plaquette term and the loops take values in Q(sqrt2). Every computation is a job of the HEB-C HTTP-0
//! backend (hebc_ym_http0 serve) over HTTP/1.1. See PREREG.md and ADDENDUM.md.
//! Rust 1.81, integer arithmetic only (i128, big integers, exact Q(sqrt2) comparisons), json=0.
//! usage: z2o PORT
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::net::TcpStream;

const JE: i128 = 810_000_000_000_000;
const UNIT: i128 = 1000 * JE;
const LAM: i128 = 3_240_000; // lambda = 4 * 10^-9
const M_FLOOR: i128 = 810_000_000_000; // 10^-6 J_E
const NP: i128 = 2;
const COLS: usize = 3;
const RS: [(i128, i128); 4] = [(1, 27), (1, 9), (1, 3), (1, 1)];
const LOOPS: [usize; 2] = [1, 2];
const ORDER: usize = 5;
const G: i64 = 4; // 2O
const SEALED: [(&str, &str, i64, usize, i128, &str); 3] = [
    ("../zwilson-2026-10-30-cloud/RESULT.hbp", "8fe1c2a4d8fa0a7058d290f0ddc020a22b0422ca9bb888db1f32a6e1a794c818", 1, 4, 270_000_000_000_000, "S3"),
    ("../zq8-2026-10-30-cloud/RESULT.hbp", "b6c129eb5a9982e5014622fb0fead32e0012372b9a82b3af166c88cd07f1d50b", 2, 4, 1_080_000_000_000_000, "Q8"),
    ("../zt2-2026-10-30-cloud/RESULT.hbp", "45040a1708c015442c99f0eebd78240f0de2686884427afa11eda96667e11290", 3, 3, 1_350_000_000_000_000, "2T"),
];
const ZQ8_MATRIX_SPEC: &str = "fgwilsonmat/2/2/4/0/1/40000000000000/0";
const ZQ8_MATRIX_HBP: &str = "2b79db3ff77fc1f491d1d4ad661a0b93639b25efd9715c446ded5137301bf038";
const ZT2_MATRIX_SPEC: &str = "fgwilsonmat/3/2/3/0/1/50000000000000/0";
const ZT2_MATRIX_HBP: &str = "1250860aaa1a6e0d9098628aa44c8a8924b7fd447444a133f9230dc875a5c456";

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

/// a + b sqrt2 with rational a, b (big integers, unreduced, positive denominators)
#[derive(Clone, Debug)]
struct Q2 {
    a: Q,
    b: Q,
}

fn q2(a: Q, b: Q) -> Q2 {
    Q2 { a, b }
}

fn qr(n: i128, d: i128) -> Q2 {
    q2(q(n, d), q(0, 1))
}

impl Q2 {
    /// "n/d" or "n/d+n'/d'r2" or "n/d-n'/d'r2" as printed by fgwser
    fn parse(s: &str) -> Option<Q2> {
        match s.strip_suffix("r2") {
            None => Some(q2(Q::parse(s)?, q(0, 1))),
            Some(t) => {
                let k = t.char_indices().skip(1).filter(|&(_, c)| c == '+' || c == '-').map(|(k, _)| k).last()?;
                let (x, y) = t.split_at(k);
                let y = y.strip_prefix('+').unwrap_or(y);
                Some(q2(Q::parse(x)?, Q::parse(y)?))
            }
        }
    }
    fn add(&self, o: &Q2) -> Q2 {
        q2(self.a.add(&o.a), self.b.add(&o.b))
    }
    fn sub(&self, o: &Q2) -> Q2 {
        q2(self.a.sub(&o.a), self.b.sub(&o.b))
    }
    fn mul(&self, o: &Q2) -> Q2 {
        q2(self.a.mul(&o.a).add(&self.b.mul(&o.b).mul(&q(2, 1))), self.a.mul(&o.b).add(&self.b.mul(&o.a)))
    }
    fn div(&self, o: &Q2) -> Q2 {
        // (a + b r)^-1 = (a - b r) / (a^2 - 2 b^2)
        let conj = q2(o.a.clone(), q(0, 1).sub(&o.b));
        let norm = o.a.mul(&o.a).sub(&o.b.mul(&o.b).mul(&q(2, 1)));
        let num = self.mul(&conj);
        q2(num.a.div(&norm), num.b.div(&norm))
    }
    fn pow(&self, k: usize) -> Q2 {
        (0..k).fold(qr(1, 1), |acc, _| acc.mul(self))
    }
    /// exact sign of a + b sqrt2
    fn sign(&self) -> Ordering {
        let sa = (&self.a.n * &self.a.d).sign();
        let sb = (&self.b.n * &self.b.d).sign();
        use num_bigint::Sign::{Minus, NoSign, Plus};
        match (sa, sb) {
            (NoSign, NoSign) => Ordering::Equal,
            (Plus, Plus) | (Plus, NoSign) | (NoSign, Plus) => Ordering::Greater,
            (Minus, Minus) | (Minus, NoSign) | (NoSign, Minus) => Ordering::Less,
            _ => {
                let a2 = self.a.mul(&self.a);
                let b2 = self.b.mul(&self.b).mul(&q(2, 1));
                let c = (&a2.n * &b2.d).cmp(&(&b2.n * &a2.d));
                match (sa, c) {
                    (Plus, Ordering::Greater) | (Minus, Ordering::Less) => Ordering::Greater,
                    (Plus, Ordering::Less) | (Minus, Ordering::Greater) => Ordering::Less,
                    _ => Ordering::Equal,
                }
            }
        }
    }
    fn lt(&self, o: &Q2) -> bool {
        self.sub(o).sign() == Ordering::Less
    }
    fn eq(&self, o: &Q2) -> bool {
        self.sub(o).sign() == Ordering::Equal
    }
    fn abs(&self) -> Q2 {
        if self.sign() == Ordering::Less {
            qr(0, 1).sub(self)
        } else {
            self.clone()
        }
    }
    fn is_zero(&self) -> bool {
        self.sign() == Ordering::Equal
    }
    /// 15 decimals, display only: floor of 10^15 a plus floor of 10^15 b sqrt2 (integer square root)
    fn dec(&self) -> String {
        let s = BigInt::from(10).pow(15);
        let x = num_integer::Integer::div_floor(&(&self.a.n * &s), &self.a.d);
        let y = &self.b.n * &s; // b sqrt2 10^15 = y sqrt2 / bd
        let ys = (BigInt::from(2) * &y * &y).sqrt();
        let ys = if y.is_negative() { -ys - 1 } else { ys };
        let v = x + num_integer::Integer::div_floor(&ys, &self.b.d);
        let neg = v.is_negative();
        let a = v.abs();
        let (ip, fp) = (num_integer::Integer::div_floor(&a, &s), num_integer::Integer::mod_floor(&a, &s));
        format!("{}{}.{:0>15}", if neg { "-" } else { "" }, ip, fp.to_string())
    }
    fn show(&self) -> String {
        if self.b.n.is_zero() {
            format!("{}/{}", self.a.n, self.a.d)
        } else {
            format!("{}/{}{}{}/{}r2", self.a.n, self.a.d, if self.b.n.is_negative() { "" } else { "+" }, self.b.n, self.b.d)
        }
    }
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
    fn matrix_spec(&mut self, spec: &str) -> bool {
        if let Some(v) = self.cache.get(spec) {
            return v.is_some();
        }
        let r = self.job(spec);
        let ok = matches!(&r, Ok((_, s)) if s.contains("w_self_adjoint=PASS"));
        match &r {
            Ok((hs, shard)) => {
                let fm = shard.lines().find(|l| l.starts_with("FG_MATRIX|")).unwrap_or("FG_MATRIX|missing");
                self.out.push(format!("MATRIX|spec={spec}|ok={}|states={}|w_self_adjoint={}|hbp_sha256={hs}|json=0", u8::from(ok), field(fm, "states=").unwrap_or("-"), field(fm, "w_self_adjoint=").unwrap_or("-")));
            }
            Err(e) => self.out.push(format!("MATRIX|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0")),
        }
        self.cache.insert(spec.to_string(), if ok { Some(1) } else { None });
        ok
    }
    /// sturmq2 negative count through fgwilsoncountq2 (mode 0 exact, 1 interval), None unless COMPLETE
    fn count_spec(&mut self, spec: &str) -> Option<i64> {
        if let Some(v) = self.cache.get(spec) {
            return *v;
        }
        let v = match self.job(spec) {
            Ok((hs, shard)) => {
                let line = shard.lines().find(|l| l.starts_with("STURM-Q2|")).unwrap_or("");
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
        self.cache.insert(spec.to_string(), v);
        v
    }
    fn matrix(&mut self, x1: usize, jb: i128, lam: i128) -> bool {
        self.matrix_spec(&format!("fgwilsonmat/{G}/2/{COLS}/0/{x1}/{jb}/{lam}"))
    }
    fn count(&mut self, x1: usize, jb: i128, lam: i128, t: i128) -> Option<i64> {
        self.count_spec(&format!("fgwilsoncountq2/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{t}/1"))
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
        let a = self.verify(&format!("fgwilsoncountq2/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{lo}/1"));
        let b = self.verify(&format!("fgwilsoncountq2/{G}/2/{COLS}/0/{x1}/{jb}/{lam}/{JE}/{hi}/1"));
        (a && b).then_some((lo, hi))
    }
}

fn series_of(shard: &str) -> BTreeMap<usize, Vec<Q2>> {
    let mut m = BTreeMap::new();
    for l in shard.lines() {
        if let (Some(lp), Some(cs)) = (field(l, "loop=1x"), field(l, "c_k=")) {
            let v: Option<Vec<Q2>> = cs.split(',').map(Q2::parse).collect();
            if let (Ok(a), Some(v)) = (lp.parse::<usize>(), v) {
                m.insert(a, v);
            }
        }
    }
    m
}

type SeriesOut = (bool, BTreeMap<usize, Vec<Q2>>, Vec<Q2>);

fn series_job(c: &mut Client, spec: &str) -> SeriesOut {
    match c.job(spec) {
        Ok((hs, shard)) => {
            for l in shard.lines().filter(|l| l.starts_with("FGSERIES")) {
                c.out.push(l.to_string());
            }
            let basis_ok = shard.lines().any(|l| l.starts_with("FGSERIES_BASIS|") && l.contains("e_w_self_adjoint=PASS") && l.contains("rhs_orthogonal_and_solved=PASS") && l.contains("orthogonality=PASS"));
            let e: Vec<Q2> = shard.lines().find_map(|l| field(l, "E_k=")).map(|e| e.split(',').filter_map(Q2::parse).collect()).unwrap_or_default();
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

fn read_sealed(path: &str, sha: &str) -> (bool, String) {
    let bytes = std::fs::read(path).unwrap_or_default();
    let ok = sha256_hex(&bytes) == sha;
    (ok, String::from_utf8(bytes).unwrap_or_default())
}

fn main() {
    let port: u16 = std::env::args().nth(1).and_then(|p| p.parse().ok()).unwrap_or(0);
    let mut c = Client { port, out: Vec::new(), cache: BTreeMap::new() };
    c.out.push(format!("Z2O|group=2O|electric_set=(1+-i)/r2,(1+-j)/r2,(1+-k)/r2|geometry=S2_open_2x3|field=Q(sqrt2)|JE={JE}|lambda=LAM/JE={LAM}/{JE}|r=1/27,1/9,1/3,1|counter=sturmq2_interval|json=0"));

    // G1-G3
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
            c.out.push(format!("CHECK_JOB|spec={check}|hbp_sha256={hs}|charpoly_is_x(x-4)^9(x-6)^20(x-8)^9(x-12)(x^2-12x+18)^4={}|verify={}|json=0", u8::from(reg), if ver { "PASS" } else { "FAIL" }));
        }
        Err(e) => c.out.push(format!("CHECK_JOB|result=INCONCLUSIVE|reason={e}|json=0")),
    }

    // R1: S3, Q8, 2T series and matrices against the sealed packages
    let mut r1 = true;
    let series_regress: [(usize, i64, usize, usize, &str, &str); 3] = [
        (0, 1, 4, 4, "WSERIES|loop=", "method=A_orbit"),
        (1, 2, 4, 5, "FGSERIES|group=Q8|loop=", "method=A_orbit_elimination"),
        (2, 3, 3, 5, "FGSERIES|group=2T|loop=", "method=A_orbit_elimination"),
    ];
    for (si, group, cols, order, prefix, tag) in series_regress {
        let (sealed_ok, text) = read_sealed(SEALED[si].0, SEALED[si].1);
        let want: Vec<String> = text.lines().filter(|l| l.starts_with(prefix) && l.contains(tag)).filter_map(|l| field(l, "c_k=").map(str::to_string)).collect();
        for method in 0..2 {
            let spec = format!("fgwilsonseries/{group}/2/{cols}/{method}/{order}");
            match c.job(&spec) {
                Ok((hs, shard)) => {
                    let got: Vec<String> = shard.lines().filter(|l| l.starts_with("FGSERIES|")).filter_map(|l| field(l, "c_k=").map(str::to_string)).collect();
                    let same = sealed_ok && !want.is_empty() && got == want;
                    let ver = c.verify(&spec);
                    r1 &= same && ver;
                    c.out.push(format!("REGRESSION_SERIES|group={}|spec={spec}|identical_to_sealed={}|sealed_result_sha_ok={}|hbp_sha256={hs}|verify={}|json=0", SEALED[si].5, u8::from(same), u8::from(sealed_ok), if ver { "PASS" } else { "FAIL" }));
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
        for (spec, want) in [(ZQ8_MATRIX_SPEC, ZQ8_MATRIX_HBP), (ZT2_MATRIX_SPEC, ZT2_MATRIX_HBP)] {
            let hs = c.job(spec).map(|(h, _)| h).unwrap_or_default();
            let same = hs == want;
            r1 &= same;
            c.out.push(format!("REGRESSION_MATRIX_SHARD|spec={spec}|hbp_sha256={hs}|sealed_hbp_sha256={want}|identical={}|json=0", u8::from(same)));
        }
    }

    // S1: sturmq2 exact and interval against every sealed E0(0) bracket
    let mut s1 = true;
    for (path, sha, group, cols, je, name) in SEALED {
        let (sealed_ok, text) = read_sealed(path, sha);
        s1 &= sealed_ok;
        let brackets: Vec<(i128, i128, i128, i128)> = text
            .lines()
            .filter(|l| l.starts_with("ENERGY0|"))
            .filter_map(|l| {
                let (rn, rd) = field(l, "r=")?.split_once('/')?;
                Some((rn.parse().ok()?, rd.parse().ok()?, field(l, "T_lo=")?.parse().ok()?, field(l, "T_hi=")?.parse().ok()?))
            })
            .collect();
        s1 &= brackets.len() == 4;
        for (rn, rd, tlo, thi) in brackets {
            let jb = je * rn / rd;
            let mspec = format!("fgwilsonmat/{group}/2/{cols}/0/1/{jb}/0");
            if !c.matrix_spec(&mspec) {
                s1 = false;
                continue;
            }
            for mode in [0, 1] {
                let lo = c.count_spec(&format!("fgwilsoncountq2/{group}/2/{cols}/0/1/{jb}/0/{je}/{tlo}/{mode}"));
                let hi = c.count_spec(&format!("fgwilsoncountq2/{group}/2/{cols}/0/1/{jb}/0/{je}/{thi}/{mode}"));
                let ok = lo == Some(0) && hi.is_some_and(|v| v >= 1);
                s1 &= ok;
                c.out.push(format!("S1|group={name}|r={rn}/{rd}|mode={}|negative_at_T_lo={}|negative_at_T_hi={}|reproduces_sturmrs={}|json=0", if mode == 0 { "exact" } else { "interval" }, lo.map_or(-1, |v| v), hi.map_or(-1, |v| v), u8::from(ok)));
            }
        }
    }

    // W1: exact 2O series, elimination (136 orbits) and configuration spectral (2304 configurations)
    let predicted_lead = [q2(q(1, 6), q(1, 12)), q2(q(7, 288), q(7, 432))];
    let mut w1 = true;
    let a = series_job(&mut c, &format!("fgwilsonseries/{G}/2/{COLS}/0/{ORDER}"));
    let b = series_job(&mut c, &format!("fgwilsonseries/{G}/2/{COLS}/1/{ORDER}"));
    w1 &= a.0 && b.0;
    let mut lead: BTreeMap<usize, Q2> = BTreeMap::new();
    let mut full: BTreeMap<usize, Vec<Q2>> = BTreeMap::new();
    for x in 1..COLS {
        let (Some(sa), Some(sb)) = (a.1.get(&x), b.1.get(&x)) else {
            w1 = false;
            continue;
        };
        let same = sa.len() == ORDER + 1 && sb.len() == ORDER + 1 && sa.iter().zip(sb).all(|(u, v)| u.eq(v));
        let below = (0..x).all(|k| sa[k].is_zero());
        let lead_ok = sa[x].eq(&predicted_lead[x - 1]);
        let next_zero = sa[x + 1].is_zero();
        w1 &= same && below && lead_ok && next_zero;
        lead.insert(x, sa[x].clone());
        full.insert(x, sa.clone());
        c.out.push(format!(
            "SERIES_CHECK|loop=1x{x}|methods_identical={}|zero_below_order_A={}|c_AA={}|c_AA_dec={}|c_AA_equals_predicted={}|c_A_A+1={}|centre_symmetry_zero={}|json=0",
            u8::from(same),
            u8::from(below),
            sa[x].show(),
            sa[x].dec(),
            u8::from(lead_ok),
            sa[x + 1].show(),
            u8::from(next_zero)
        ));
    }
    let e2 = q2(q(-1, 6), q(-1, 12));
    let e_ok = a.2.len() == ORDER && b.2.len() == ORDER && a.2.iter().zip(&b.2).all(|(u, v)| u.eq(v)) && a.2[0].eq(&qr(4, 1)) && a.2[1].eq(&e2) && a.2[2].is_zero() && a.2[4].is_zero();
    w1 &= e_ok;
    c.out.push(format!("SERIES_ENERGY_CHECK|E1=4|E2=-1/6-1/12r2|E3=0|E5=0|methods_identical_and_registered={}|json=0", u8::from(e_ok)));

    // certified values (interval sturmq2)
    let mut wv: BTreeMap<(usize, usize), (Q, Q)> = BTreeMap::new();
    let mut w2 = true;
    let mut inconclusive: Vec<String> = Vec::new();
    let margin = (2 * 1000 * LAM * 99 + 69) / 70 + 1; // 2 sqrt2 lambda, sqrt2 < 99/70, rounded outward
    let mut x1: Option<bool> = None;
    for (ri, &(rn, rd)) in RS.iter().enumerate() {
        let jb = JE * rn / rd;
        let second = -(1000 * NP * jb * jb / (24 * JE)) * 239 / 70;
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
        c.out.push(format!("ENERGY0|r={rn}/{rd}|T_lo={l0}|T_hi={u0}|E0_lo={}|json=0", Q2 { a: q(l0, UNIT), b: q(0, 1) }.dec()));
        if ri == 0 {
            // X1 (addendum 2): exact sturmq2 against the interval counts on the 2O matrix itself
            let ex_lo = c.count_spec(&format!("fgwilsoncountq2/{G}/2/{COLS}/0/1/{jb}/0/{JE}/{l0}/0"));
            let ex_hi = c.count_spec(&format!("fgwilsoncountq2/{G}/2/{COLS}/0/1/{jb}/0/{JE}/{u0}/0"));
            x1 = match (ex_lo, ex_hi) {
                (Some(a), Some(b)) => Some(a == 0 && b >= 1),
                _ => None,
            };
            c.out.push(format!("X1|r={rn}/{rd}|exact_negative_at_T_lo={}|exact_negative_at_T_hi={}|interval_negative_at_T_lo=0|interval_negative_at_T_hi>=1|agree={}|json=0", ex_lo.map_or(-1, |v| v), ex_hi.map_or(-1, |v| v), x1.map_or("INCONCLUSIVE".to_string(), |b| u8::from(b).to_string())));
        }
        for &aa in &LOOPS {
            let mut ends = Vec::new();
            for lam in [LAM, -LAM] {
                let got = (|| -> Option<(i128, i128)> {
                    if !c.matrix(aa, jb, lam) {
                        return None;
                    }
                    let (lo, hi) = (l0 - margin, u0 + margin);
                    if c.count(aa, jb, lam, lo)? == 0 && c.count(aa, jb, lam, hi)? >= 1 {
                        c.bisect(aa, jb, lam, lo, hi)
                    } else {
                        None
                    }
                })();
                ends.push(got);
            }
            let (Some((lp, _)), Some((lm, _))) = (ends[0], ends[1]) else {
                w2 = false;
                inconclusive.push(format!("1x{aa}:r={rn}/{rd}"));
                continue;
            };
            let (lo, hi) = (q(lp - u0, 1000 * LAM), q(u0 - lm, 1000 * LAM));
            let width_t = u0 - lm - (lp - u0);
            let width_ok = width_t * 1_000_000_000 <= 3 * 1000 * LAM;
            let pos = lp - u0 > 0;
            w2 &= width_ok && pos;
            c.out.push(format!(
                "WILSON|loop=1x{aa}|r={rn}/{rd}|W_lo={}/{}|W_hi={}/{}|W_lo_dec={}|W_hi_dec={}|width_T={width_t}|width_le_3e-9={}|positive={}|json=0",
                lo.n,
                lo.d,
                hi.n,
                hi.d,
                Q2 { a: lo.clone(), b: q(0, 1) }.dec(),
                Q2 { a: hi.clone(), b: q(0, 1) }.dec(),
                u8::from(width_ok),
                u8::from(pos)
            ));
            wv.insert((aa, ri), (lo, hi));
        }
    }

    // W3, W4, W5, W6 with exact Q(sqrt2) comparisons
    let rq = |ri: usize| qr(RS[ri].0, RS[ri].1);
    let lift = |x: &Q| Q2 { a: x.clone(), b: q(0, 1) };
    let one = qr(1, 1);
    let mut w3 = true;
    let (mut w4_resolved, mut w4_fail) = (0usize, false);
    for &aa in &LOOPS {
        let (Some(ca), Some(fs)) = (lead.get(&aa), full.get(&aa)) else {
            w3 = false;
            w4_fail = true;
            continue;
        };
        let ratio = |ri: usize| -> Option<(Q2, Q2)> {
            let (lo, hi) = wv.get(&(aa, ri))?;
            let p = ca.mul(&rq(ri).pow(aa));
            Some((lift(lo).div(&p), lift(hi).div(&p)))
        };
        let (Some((r27l, r27h)), Some((r9l, r9h))) = (ratio(0), ratio(1)) else {
            w3 = false;
            w4_fail = true;
            continue;
        };
        let ok3 = !r27l.lt(&qr(4, 5)) && !qr(5, 4).lt(&r27h);
        w3 &= ok3;
        let delta = fs[aa + 2].div(ca).abs().mul(&qr(1, 81));
        let resolved = r27h.sub(&r27l).lt(&delta) && r9h.sub(&r9l).lt(&delta);
        let maxd = |l: &Q2, h: &Q2| -> Q2 {
            let (x, y) = (l.sub(&one).abs(), h.sub(&one).abs());
            if x.lt(&y) {
                y
            } else {
                x
            }
        };
        let mind = |l: &Q2, h: &Q2| -> Q2 {
            if !one.lt(l) && !h.lt(&one) {
                qr(0, 1)
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
            "LEADING|loop=1x{aa}|ratio_r=1/27=[{},{}]|ratio_r=1/9=[{},{}]|within_0.8_1.25={}|series_delta_at_1/9={}|resolved={}|converging={}|w4={}|json=0",
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
            let t1 = ca.mul(&r.pow(aa));
            let t2 = t1.add(&fs[aa + 1].mul(&r.pow(aa + 1))).add(&fs[aa + 2].mul(&r.pow(aa + 2)));
            if let Some((lo, hi)) = wv.get(&(aa, ri)) {
                let inside = !t2.lt(&lift(lo)) && !lift(hi).lt(&t2);
                c.out.push(format!("SERIES_COMPARE|loop=1x{aa}|r={rn}/{rd}|leading={}|through_order_A+2={}|inside={}|json=0", t1.dec(), t2.dec(), u8::from(inside)));
            }
        }
    }
    let mut w5 = true;
    for (ri, &(rn, rd)) in RS.iter().enumerate() {
        let (Some((al, ah)), Some((bl, bh))) = (wv.get(&(1, ri)), wv.get(&(2, ri))) else {
            w5 = false;
            continue;
        };
        let dec = bh.lt(al);
        w5 &= dec;
        let (rl, rh) = (lift(bl).div(&lift(ah)), lift(bh).div(&lift(al)));
        let mut line = format!("AREA|loops=1x2/1x1|r={rn}/{rd}|ratio=[{},{}]|decreasing={}", rl.dec(), rh.dec(), u8::from(dec));
        if ri == 0 {
            if let (Some(ca), Some(cb)) = (lead.get(&1), lead.get(&2)) {
                let pred = cb.div(ca).mul(&rq(0));
                let (lo, hi) = (rl.div(&pred), rh.div(&pred));
                let ok = !lo.lt(&qr(4, 5)) && !qr(5, 4).lt(&hi);
                w5 &= ok;
                line.push_str(&format!("|predicted_ratio={}|within_0.8_1.25={}", pred.dec(), u8::from(ok)));
            } else {
                w5 = false;
            }
        }
        c.out.push(format!("{line}|json=0"));
    }
    let w6 = match (wv.get(&(1, 0)), wv.get(&(1, 1)), lead.get(&1)) {
        (Some((l27, h27)), Some((l9, h9)), Some(c11)) => {
            let d = |l: &Q, h: &Q, ri: usize| -> (Q2, Q2) {
                let p = c11.mul(&rq(ri));
                (lift(l).div(&p).sub(&one), lift(h).div(&p).sub(&one))
            };
            let (d27l, d27h) = d(l27, h27, 0);
            let (d9l, d9h) = d(l9, h9, 1);
            let z = qr(0, 1);
            let pos = z.lt(&d27l) && z.lt(&d9l);
            let neg = d27h.lt(&z) && d9h.lt(&z);
            let (rlo, rhi) = if pos {
                (d9l.div(&d27h), d9h.div(&d27l))
            } else if neg {
                (d9h.div(&d27l), d9l.div(&d27h))
            } else {
                (z.clone(), z.clone())
            };
            let inside = |a: i128, b: i128| !rlo.lt(&qr(a, 1)) && !qr(b, 1).lt(&rhi);
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
        "Z2O_RESULT|G1={}|G2={}|G3={}|G4={}|R1={}|S1={}|X1={}|W1={}|W2={}|W3={}|W4={w4}|W4_resolved_loops={w4_resolved}|W5={}|W6={}|inconclusive={}|json=0",
        p(g123.0),
        p(g123.1),
        p(g123.2),
        p(c.out.iter().filter(|l| l.starts_with("MATRIX|")).all(|l| l.contains("|w_self_adjoint=PASS|")) && c.out.iter().any(|l| l.starts_with("MATRIX|"))),
        p(r1),
        p(s1),
        match x1 {
            Some(b) => p(b),
            None => "INCONCLUSIVE",
        },
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

/// x (x - 4)^9 (x - 6)^20 (x - 8)^9 (x - 12) (x^2 - 12 x + 18)^4, coefficients low to high, as written by fgcheck
fn registered_charpoly() -> String {
    let mut p: Vec<i128> = vec![1];
    let mut mulby = |f: &[i128]| {
        let mut n = vec![0i128; p.len() + f.len() - 1];
        for (i, c) in p.iter().enumerate() {
            for (j, d) in f.iter().enumerate() {
                n[i + j] += c * d;
            }
        }
        p = n;
    };
    mulby(&[0, 1]);
    for (root, mult) in [(4i128, 9), (6, 20), (8, 9), (12, 1)] {
        for _ in 0..mult {
            mulby(&[-root, 1]);
        }
    }
    for _ in 0..4 {
        mulby(&[18, -12, 1]);
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

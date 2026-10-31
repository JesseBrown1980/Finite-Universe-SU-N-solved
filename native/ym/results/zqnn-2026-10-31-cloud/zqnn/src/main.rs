//! zqnn driver: the qutrit variational network on the finite-group lane. Upper bounds U from fgqnn jobs,
//! certified lower bounds L from the sealed fgwilsonmat/fgwilsoncount jobs (vacuum, LAM = 0), for S3, Q8 and
//! 2T on the open 2 x 3 strip; upper bounds on the open 2 x 4 ladder against the sealed zwilson and zq8
//! certificates and the sealed 2T series. Every computation is a job of the HEB-C HTTP-0 backend
//! (hebc_ym_http0 serve) over HTTP/1.1. See PREREG.md and ADDENDUM.md.
//! Rust 1.81, integer arithmetic only (i128 and big integers), json=0.
//! usage: zqnn PORT
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::net::TcpStream;

const JE_STRONG: i128 = 1_350_000_000_000_000;
const JE_WEAK: i128 = 135_000_000_000_000;
const RS2: [(i128, i128); 6] = [(1, 27), (1, 9), (1, 3), (1, 1), (3, 1), (9, 1)];
const RS3: [(i128, i128); 3] = [(1, 27), (1, 3), (1, 1)];
const GROUPS: [(i64, &str); 3] = [(1, "S3"), (2, "Q8"), (3, "2T")];
const ZWILSON_RESULT: &str = "../zwilson-2026-10-30-cloud/RESULT.hbp";
const ZWILSON_RESULT_SHA: &str = "8fe1c2a4d8fa0a7058d290f0ddc020a22b0422ca9bb888db1f32a6e1a794c818";
const ZWILSON_JE: i128 = 270_000_000_000_000;
const ZQ8_RESULT: &str = "../zq8-2026-10-30-cloud/RESULT.hbp";
const ZQ8_RESULT_SHA: &str = "b6c129eb5a9982e5014622fb0fead32e0012372b9a82b3af166c88cd07f1d50b";
const ZQ8_JE: i128 = 1_080_000_000_000_000;
const ZT2_RESULT: &str = "../zt2-2026-10-30-cloud/RESULT.hbp";
const ZT2_RESULT_SHA: &str = "45040a1708c015442c99f0eebd78240f0de2686884427afa11eda96667e11290";

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
    inconclusive: Vec<String>,
}

struct Qnn {
    start: Q,
    energy: Q,
    exact: bool,
    line: String,
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
    fn qnn(&mut self, g: i64, cols: usize, r: (i128, i128), seed: i64, sweeps: i64) -> Option<Qnn> {
        let spec = format!("fgqnn/{g}/2/{cols}/{}/{}/{seed}/{sweeps}", r.0, r.1);
        match self.job(&spec) {
            Ok((hs, shard)) => {
                let head = shard.lines().find(|l| l.starts_with("FG_QNN|")).unwrap_or("").to_string();
                let get = |p: &str| -> Option<Q> {
                    let l = shard.lines().find(|l| l.starts_with(p))?;
                    Some(Q { n: field(l, "numerator=")?.parse().ok()?, d: field(l, "denominator=")?.parse().ok()? })
                };
                let (Some(start), Some(energy)) = (get("FG_QNN_START|"), get("FG_QNN_ENERGY|")) else {
                    self.inconclusive.push(format!("qnn_{spec}"));
                    return None;
                };
                let exact = ["norm_exact=1", "inverse_exact=1", "recomputed_equal=1", "nonnegative=1"].iter().all(|k| head.contains(&format!("|{k}|")));
                let b = shard.lines().find(|l| l.starts_with("FG_QNN_B|")).unwrap_or("");
                self.out.push(format!("QNN|spec={spec}|{}|hbp_sha256={hs}|json=0", head.trim_start_matches("FG_QNN|").trim_end_matches("|json=0")));
                self.out.push(format!("QNN_B|spec={spec}|{}|json=0", b.trim_start_matches("FG_QNN_B|").trim_end_matches("|json=0")));
                self.out.push(format!("QNN_ENERGY|spec={spec}|start={}|U={}|U_exact={}/{}|json=0", start.dec(), energy.dec(), energy.n, energy.d));
                Some(Qnn { start, energy, exact, line: spec })
            }
            Err(e) => {
                self.out.push(format!("QNN|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0"));
                self.inconclusive.push(format!("qnn_{spec}"));
                None
            }
        }
    }
    fn matrix(&mut self, g: i64, jb: i128) -> bool {
        let spec = format!("fgwilsonmat/{g}/2/3/0/1/{jb}/0");
        if let Some(v) = self.cache.get(&spec) {
            return v.is_some();
        }
        let r = self.job(&spec);
        let ok = matches!(&r, Ok((_, s)) if s.contains("w_self_adjoint=PASS"));
        match &r {
            Ok((hs, shard)) => {
                let fm = shard.lines().find(|l| l.starts_with("FG_MATRIX|")).unwrap_or("FG_MATRIX|missing");
                self.out.push(format!("MATRIX|spec={spec}|ok={}|states={}|hbp_sha256={hs}|json=0", u8::from(ok), field(fm, "states=").unwrap_or("-")));
            }
            Err(e) => self.out.push(format!("MATRIX|spec={spec}|result=INCONCLUSIVE|reason={e}|json=0")),
        }
        self.cache.insert(spec, if ok { Some(1) } else { None });
        ok
    }
    fn count(&mut self, g: i64, jb: i128, je: i128, t: i128) -> Option<i64> {
        let spec = format!("fgwilsoncount/{g}/2/3/0/1/{jb}/0/{je}/{t}");
        if let Some(v) = self.cache.get(&spec) {
            return *v;
        }
        let v = match self.job(&spec) {
            Ok((_, shard)) => {
                let line = shard.lines().find(|l| l.starts_with("STURM-RS|")).unwrap_or("");
                match (field(line, "status="), field(line, "negative=").and_then(|v| v.parse::<i64>().ok())) {
                    (Some("COMPLETE"), Some(n)) => Some(n),
                    _ => None,
                }
            }
            Err(_) => None,
        };
        self.cache.insert(spec, v);
        v
    }
}

fn ceil_div(a: &BigInt, b: &BigInt) -> BigInt {
    -num_integer::Integer::div_floor(&(-a), b)
}

fn to_i128(x: &BigInt) -> Option<i128> {
    x.to_string().parse().ok()
}

/// the sealed ENERGY0 rows of a package: r -> T_lo
fn sealed_energy(path: &str, sha: &str) -> Option<BTreeMap<String, i128>> {
    let bytes = std::fs::read(path).ok()?;
    if sha256_hex(&bytes) != sha {
        return None;
    }
    let text = String::from_utf8(bytes).ok()?;
    Some(
        text.lines()
            .filter(|l| l.starts_with("ENERGY0|"))
            .filter_map(|l| Some((field(l, "r=")?.to_string(), field(l, "T_lo=")?.parse().ok()?)))
            .collect(),
    )
}

/// the sealed 2T series on the 2 x 4 ladder: E_k
fn zt2_series() -> Option<Vec<Q>> {
    let bytes = std::fs::read(ZT2_RESULT).ok()?;
    if sha256_hex(&bytes) != ZT2_RESULT_SHA {
        return None;
    }
    let text = String::from_utf8(bytes).ok()?;
    let l = text.lines().find(|l| l.starts_with("FGSERIES_ENERGY|group=2T|E_k=6/1"))?;
    field(l, "E_k=")?.split(',').map(Q::parse).collect()
}

#[allow(clippy::too_many_lines)]
fn main() {
    let port: u16 = std::env::args().nth(1).and_then(|p| p.parse().ok()).unwrap_or(0);
    let mut c = Client { port, out: Vec::new(), cache: BTreeMap::new(), inconclusive: Vec::new() };
    c.out.push(format!(
        "ZQNN|network=v7_qutrit_rotations_a=6561_eight_trits_per_b_two_layers|groups=S3,Q8,2T|S2=open_2x3_r=1/27,1/9,1/3,1,3,9_seeds=UNIFORM,IDENTITY_sweeps=4|S3=open_2x4_r=1/27,1/3,1_seed=UNIFORM_sweeps=2|JE_r_le_1={JE_STRONG}|JE_r_ge_3={JE_WEAK}|json=0"
    ));
    let (mut n1, mut n2, mut n3, mut n4, mut n5, mut n6) = (true, true, true, true, true, true);
    for &(g, gname) in &GROUPS {
        for &(rn, rd) in &RS2 {
            let je = if rn >= 3 * rd { JE_WEAK } else { JE_STRONG };
            let unit = BigInt::from(1000 * je);
            let jb = je * rn / rd;
            let a = c.qnn(g, 3, (rn, rd), 0, 4);
            let b = c.qnn(g, 3, (rn, rd), 1, 4);
            let (Some(a), Some(b)) = (a, b) else { continue };
            n1 &= a.exact && b.exact;
            let (la, lb) = (a.line.clone(), b.line.clone());
            n2 &= c.verify(&la) & c.verify(&lb);
            if rn >= 3 * rd {
                let ok = b.energy.lt(&b.start);
                n5 &= ok;
                c.out.push(format!("WEAK|group={gname}|r={rn}/{rd}|identity_start={}|identity_end={}|lowered={}|json=0", b.start.dec(), b.energy.dec(), u8::from(ok)));
            }
            let best = if a.energy.lt(&b.energy) || a.energy.eq(&b.energy) { &a } else { &b };
            let (u, s, seed) = (best.energy.clone(), best.start.clone(), if std::ptr::eq(best, &a) { "UNIFORM" } else { "IDENTITY" });
            if rn <= rd {
                let ok = u.lt(&q(4 * rn, rd));
                n4 &= ok;
                c.out.push(format!("STRONG|group={gname}|geometry=2x3|r={rn}/{rd}|U={}|uniform_2rNp={}|below={}|json=0", u.dec(), q(4 * rn, rd).dec(), u8::from(ok)));
            }
            // certified lower end
            if !c.matrix(g, jb) {
                c.inconclusive.push(format!("matrix_{gname}_{rn}_{rd}"));
                continue;
            }
            let hi_big = ceil_div(&(&u.n * &unit), &u.d);
            let Some(hi) = to_i128(&hi_big).filter(|&h| h < 9_000_000_000_000_000_000) else {
                c.inconclusive.push(format!("threshold_{gname}_{rn}_{rd}"));
                continue;
            };
            let (c0, ch) = (c.count(g, jb, je, 0), c.count(g, jb, je, hi));
            let (Some(0), Some(nh)) = (c0, ch) else {
                c.inconclusive.push(format!("count_{gname}_{rn}_{rd}"));
                continue;
            };
            let (mut lo, mut h2) = (0i128, hi);
            let mut ok_b = nh >= 1;
            while ok_b && h2 - lo > 1 {
                let mid = lo + (h2 - lo) / 2;
                match c.count(g, jb, je, mid) {
                    Some(0) => lo = mid,
                    Some(_) => h2 = mid,
                    None => ok_b = false,
                }
            }
            let v1 = c.verify(&format!("fgwilsoncount/{g}/2/3/0/1/{jb}/0/{je}/{lo}"));
            let v2 = c.verify(&format!("fgwilsoncount/{g}/2/3/0/1/{jb}/0/{je}/{h2}"));
            if !(ok_b && v1 && v2) {
                c.inconclusive.push(format!("bracket_{gname}_{rn}_{rd}"));
                continue;
            }
            let l = Q { n: BigInt::from(lo), d: unit.clone() };
            let encl = !u.lt(&l) && nh >= 1;
            n3 &= encl;
            let gap = u.sub(&l);
            let span = s.sub(&l);
            let closed = if span.n.is_zero() { q(0, 1) } else { span.sub(&gap).div(&span) };
            if rn == rd * 9 || (rn == 1 && rd == 27) {
                let ok = !q(2, 1).mul(&gap).sub(&span).n.is_positive();
                n6 &= ok;
                c.out.push(format!("HALF|group={gname}|r={rn}/{rd}|seed={seed}|U_minus_L={}|S_minus_L={}|closed_at_least_half={}|json=0", gap.dec(), span.dec(), u8::from(ok)));
            }
            c.out.push(format!(
                "ENCLOSURE|group={gname}|geometry=2x3|r={rn}/{rd}|JE={je}|T_lo={lo}|T_hi={h2}|L={}|U={}|U_minus_L={}|best_seed={seed}|start={}|fraction_closed={}|count_at_ceil_U={nh}|L_le_U={}|json=0",
                l.dec(),
                u.dec(),
                gap.dec(),
                s.dec(),
                closed.dec(),
                u8::from(encl)
            ));
        }
    }
    // the 2 x 4 ladder: upper bounds against the sealed certificates and series
    let zw = sealed_energy(ZWILSON_RESULT, ZWILSON_RESULT_SHA);
    let zq = sealed_energy(ZQ8_RESULT, ZQ8_RESULT_SHA);
    let ser = zt2_series();
    c.out.push(format!("SEALED|zwilson={}|zq8={}|zt2_series={}|json=0", u8::from(zw.is_some()), u8::from(zq.is_some()), u8::from(ser.is_some())));
    for &(g, gname) in &GROUPS {
        for &(rn, rd) in &RS3 {
            let Some(a) = c.qnn(g, 4, (rn, rd), 0, 2) else { continue };
            n1 &= a.exact;
            let la = a.line.clone();
            n2 &= c.verify(&la);
            let u = a.energy.clone();
            let uni = q(6 * rn, rd);
            let ok = u.lt(&uni);
            n4 &= ok;
            c.out.push(format!("STRONG|group={gname}|geometry=2x4|r={rn}/{rd}|U={}|uniform_2rNp={}|below={}|json=0", u.dec(), uni.dec(), u8::from(ok)));
            let key = format!("{rn}/{rd}");
            let sealed = match g {
                1 => zw.as_ref().and_then(|m| m.get(&key)).map(|&t| Q { n: BigInt::from(t), d: BigInt::from(1000 * ZWILSON_JE) }),
                2 => zq.as_ref().and_then(|m| m.get(&key)).map(|&t| Q { n: BigInt::from(t), d: BigInt::from(1000 * ZQ8_JE) }),
                _ => None,
            };
            match sealed {
                Some(l) => c.out.push(format!("LADDER|group={gname}|r={key}|U={}|sealed_L={}|U_minus_L={}|L_le_U={}|json=0", u.dec(), l.dec(), u.sub(&l).dec(), u8::from(!u.lt(&l)))),
                None => {
                    if let (3, Some(e)) = (g, ser.as_ref()) {
                        let rq = q(rn, rd);
                        let s5 = e.iter().enumerate().fold(q(0, 1), |acc, (k, ek)| acc.add(&ek.mul(&rq.pow(k + 1))));
                        c.out.push(format!("LADDER|group={gname}|r={key}|U={}|series_to_r5={}|U_minus_series={}|json=0", u.dec(), s5.dec(), u.sub(&s5).dec()));
                    }
                }
            }
        }
    }
    let inc = !c.inconclusive.is_empty();
    let v = |b: bool| if inc { "INCONCLUSIVE" } else if b { "PASS" } else { "FAIL" };
    c.out.push(format!(
        "ZQNN_RESULT|N1={}|N2={}|N3={}|N4={}|N5={}|N6={}|inconclusive={}|json=0",
        v(n1),
        v(n2),
        v(n3),
        v(n4),
        v(n5),
        v(n6),
        if inc { c.inconclusive.join(",") } else { "none".into() }
    ));
    for l in &c.out {
        println!("{l}");
    }
}

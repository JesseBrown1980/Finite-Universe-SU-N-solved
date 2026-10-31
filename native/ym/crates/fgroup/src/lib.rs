//! fgroup: finite-group lattice gauge lane for the SU(N) programme. Group data (S3 for regression against
//! the sealed S3 packages, Q8 from Gaussian-integer matrices, 2T from Hurwitz quaternions with Z[w]
//! characters, 2O from Z[sqrt2] quaternions with Z[sqrt2] characters), the open ROWS x COLS strip with a comb gauge
//! tree, gauge-orbit and configuration bases, the electric operator E = sum over links and g in S of
//! (1 - L_g), the magnetic diagonal M = sum over plaquettes (2 - chi_f(U_p)) and rectangular Wilson loops.
//! Integer arithmetic only.
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use std::collections::{BTreeMap, HashMap};

/// Quadratic field of the characters: theta^2 = p + q theta (Eisenstein w: p = -1, q = -1; sqrt2: p = 2,
/// q = 0; golden ratio: p = 1, q = 1). Groups with integer characters use the Eisenstein field with y = 0.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Field {
    pub p: i64,
    pub q: i64,
    /// name of theta in printed values
    pub theta: &'static str,
}

pub const EISENSTEIN: Field = Field { p: -1, q: -1, theta: "w" };
pub const SQRT2: Field = Field { p: 2, q: 0, theta: "r2" };

/// Algebraic integer x + y theta of the group's field (exact pairs; rational integers have y = 0).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Alg {
    pub x: i64,
    pub y: i64,
}

impl Alg {
    pub const fn int(x: i64) -> Alg {
        Alg { x, y: 0 }
    }
    pub fn plus(self, o: Alg) -> Alg {
        Alg { x: self.x + o.x, y: self.y + o.y }
    }
    pub fn minus(self, o: Alg) -> Alg {
        Alg { x: self.x - o.x, y: self.y - o.y }
    }
    pub fn scale(self, k: i64) -> Alg {
        Alg { x: self.x * k, y: self.y * k }
    }
    /// (a + b t)(c + d t) = (ac + bd p) + (ad + bc + bd q) t
    pub fn times(self, o: Alg, f: Field) -> Alg {
        Alg { x: self.x * o.x + self.y * o.y * f.p, y: self.x * o.y + self.y * o.x + self.y * o.y * f.q }
    }
    /// w^k in the Eisenstein field
    pub fn w_pow(k: i64) -> Alg {
        match k.rem_euclid(3) {
            0 => Alg::int(1),
            1 => Alg { x: 0, y: 1 },
            _ => Alg { x: -1, y: -1 },
        }
    }
    /// printed value: "x" when y = 0, otherwise "x+y<theta>" / "x-y<theta>"
    pub fn show(self, f: Field) -> String {
        if self.y == 0 {
            self.x.to_string()
        } else {
            format!("{}{}{}{}", self.x, if self.y < 0 { "-" } else { "+" }, self.y.abs(), f.theta)
        }
    }
}

/// A finite group with its electric generating set, fundamental character and full character table.
pub struct FGroup {
    pub name: &'static str,
    pub n: usize,
    pub mul: Vec<Vec<u8>>,
    pub inv: Vec<u8>,
    /// electric set S (closed under inverses, a union of conjugacy classes)
    pub electric: Vec<u8>,
    /// field of the character values
    pub field: Field,
    /// fundamental (Wilson-loop and plaquette) character, real, in the group's field
    pub chi_f: Vec<Alg>,
    /// irreducible characters (rows) in the group's field, and dimensions
    pub chars: Vec<Vec<Alg>>,
    pub dims: Vec<i64>,
    /// index of the fundamental irrep in `chars`
    pub fund: usize,
}

impl FGroup {
    pub fn m(&self, a: u8, b: u8) -> u8 {
        self.mul[a as usize][b as usize]
    }
    pub fn i(&self, a: u8) -> u8 {
        self.inv[a as usize]
    }
    /// electric eigenvalue of irrep k: |S| - sum_{g in S} chi_k(g) / d_k, an algebraic integer of the field
    /// (exact; both components of the sum must be divisible by d_k)
    pub fn e_irrep(&self, k: usize) -> Option<Alg> {
        let s = self.electric.iter().fold(Alg::int(0), |a, &g| a.plus(self.chars[k][g as usize]));
        let d = self.dims[k];
        (s.x % d == 0 && s.y % d == 0).then(|| Alg { x: self.electric.len() as i64 - s.x / d, y: -s.y / d })
    }
}

/// S3 exactly as in torusl (permutation order, composition t[x] = p[q[x]], transpositions as electric set).
pub fn s3() -> FGroup {
    const PERM: [[usize; 3]; 6] = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut mul = vec![vec![0u8; 6]; 6];
    for (a, p) in PERM.iter().enumerate() {
        for (b, q) in PERM.iter().enumerate() {
            let t = [p[q[0]], p[q[1]], p[q[2]]];
            mul[a][b] = PERM.iter().position(|x| *x == t).map_or(255, |v| v as u8);
        }
    }
    let inv = (0..6).map(|a| (0..6u8).find(|&b| mul[a][b as usize] == 0).unwrap_or(255)).collect();
    let sign = |g: usize| if [1usize, 2, 5].contains(&g) { -1 } else { 1 };
    let std = |g: usize| match g {
        0 => 2,
        3 | 4 => -1,
        _ => 0,
    };
    FGroup {
        name: "S3",
        n: 6,
        mul,
        inv,
        electric: vec![1, 2, 5],
        field: EISENSTEIN,
        chi_f: (0..6).map(|g| Alg::int(std(g))).collect(),
        chars: vec![vec![Alg::int(1); 6], (0..6).map(|g| Alg::int(sign(g))).collect(), (0..6).map(|g| Alg::int(std(g))).collect()],
        dims: vec![1, 1, 2],
        fund: 2,
    }
}

/// Gaussian integer a + b i
type Gi = (i64, i64);
type M2 = [[Gi; 2]; 2];

fn gmul(x: Gi, y: Gi) -> Gi {
    (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0)
}

fn mmul(a: &M2, b: &M2) -> M2 {
    let mut c = [[(0, 0); 2]; 2];
    for (r, row) in c.iter_mut().enumerate() {
        for (s, cell) in row.iter_mut().enumerate() {
            let p = gmul(a[r][0], b[0][s]);
            let q = gmul(a[r][1], b[1][s]);
            *cell = (p.0 + q.0, p.1 + q.1);
        }
    }
    c
}

fn mneg(a: &M2) -> M2 {
    let mut c = *a;
    for row in c.iter_mut() {
        for cell in row.iter_mut() {
            *cell = (-cell.0, -cell.1);
        }
    }
    c
}

/// The Q8 elements as 2x2 Gaussian-integer matrices, in the order 1, -1, i, -i, j, -j, k, -k (k = i j).
pub fn q8_matrices() -> Vec<M2> {
    let one: M2 = [[(1, 0), (0, 0)], [(0, 0), (1, 0)]];
    let qi: M2 = [[(0, 1), (0, 0)], [(0, 0), (0, -1)]];
    let qj: M2 = [[(0, 0), (1, 0)], [(-1, 0), (0, 0)]];
    let qk = mmul(&qi, &qj);
    vec![one, mneg(&one), qi, mneg(&qi), qj, mneg(&qj), qk, mneg(&qk)]
}

/// Q8 from its matrices: multiplication table from exact products (255 marks a product outside the set,
/// caught by the closure check), fundamental character = trace (real part; the imaginary part is checked 0).
pub fn q8() -> FGroup {
    let ms = q8_matrices();
    let n = ms.len();
    let mut mul = vec![vec![255u8; n]; n];
    for (a, x) in ms.iter().enumerate() {
        for (b, y) in ms.iter().enumerate() {
            let p = mmul(x, y);
            mul[a][b] = ms.iter().position(|z| *z == p).map_or(255, |v| v as u8);
        }
    }
    let inv = (0..n).map(|a| (0..n as u8).find(|&b| mul[a][b as usize] == 0).unwrap_or(255)).collect();
    let tr: Vec<i64> = ms.iter().map(|m| if m[0][0].1 + m[1][1].1 == 0 { m[0][0].0 + m[1][1].0 } else { i64::MIN }).collect();
    // 1-dimensional characters through Q8/{+-1} = {1, i, j, k}: +1 on the class pair named, -1 on the other two
    let one_dim = |keep: usize| -> Vec<Alg> { (0..n).map(|g| Alg::int(if g / 2 == 0 || g / 2 == keep { 1 } else { -1 })).collect() };
    FGroup {
        name: "Q8",
        n,
        mul,
        inv,
        electric: vec![2, 3, 4, 5, 6, 7],
        field: EISENSTEIN,
        chi_f: tr.iter().map(|&v| Alg::int(v)).collect(),
        chars: vec![vec![Alg::int(1); n], one_dim(1), one_dim(2), one_dim(3), tr.iter().map(|&v| Alg::int(v)).collect()],
        dims: vec![1, 1, 1, 1, 2],
        fund: 4,
    }
}

/// Hurwitz quaternion with doubled integer coordinates: q = (a + b i + c j + d k) / 2
type Hq = [i64; 4];

/// exact Hamilton product of doubled coordinates: 2 (p q) = (2p)(2q)/2, None if not divisible
fn hmul(p: &Hq, q: &Hq) -> Option<Hq> {
    let [a1, b1, c1, d1] = *p;
    let [a2, b2, c2, d2] = *q;
    let r = [
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    ];
    r.iter().all(|v| v % 2 == 0).then(|| [r[0] / 2, r[1] / 2, r[2] / 2, r[3] / 2])
}

/// The 24 elements of 2T as doubled Hurwitz coordinates: 1, -1, +-i, +-j, +-k, then (+-1 +-i +-j +-k)/2 in sign
/// order (a, b, c, d) with + before -.
pub fn t2_elements() -> Vec<Hq> {
    let mut v: Vec<Hq> = vec![[2, 0, 0, 0], [-2, 0, 0, 0], [0, 2, 0, 0], [0, -2, 0, 0], [0, 0, 2, 0], [0, 0, -2, 0], [0, 0, 0, 2], [0, 0, 0, -2]];
    for m in 0..16 {
        let sg = |bit: usize| if m & (8 >> bit) == 0 { 1 } else { -1 };
        v.push([sg(0), sg(1), sg(2), sg(3)]);
    }
    v
}

/// 2T: multiplication table from exact Hamilton products; chi_f = 2 Re(q) = a; characters in Z[w]:
/// trivial, 1' = w^k, 1'' = w^2k (k the coset of g in 2T/Q8 relative to a = (1 + i + j + k)/2), 3 = 4 Re^2 - 1,
/// fundamental, 2' = chi_f 1', 2'' = chi_f 1''. Electric set: the eight elements with Re(q) = 1/2.
pub fn t2() -> FGroup {
    let el = t2_elements();
    let n = el.len();
    let mut mul = vec![vec![255u8; n]; n];
    for (x, p) in el.iter().enumerate() {
        for (y, q) in el.iter().enumerate() {
            if let Some(r) = hmul(p, q) {
                mul[x][y] = el.iter().position(|z| *z == r).map_or(255, |v| v as u8);
            }
        }
    }
    let inv: Vec<u8> = (0..n).map(|x| (0..n as u8).find(|&y| mul[x][y as usize] == 0).unwrap_or(255)).collect();
    let chi_f: Vec<Alg> = el.iter().map(|q| Alg::int(q[0])).collect();
    // coset index k(g): g a^-k in Q8 (the first 8 elements), a = (1 + i + j + k)/2 = element 8
    let a = 8u8;
    let a_inv = inv[a as usize];
    let q8: Vec<u8> = (0..8).collect();
    let coset: Vec<i64> = (0..n as u8)
        .map(|g| {
            let mut h = g;
            for k in 0..3 {
                if q8.contains(&h) {
                    return k;
                }
                h = if (h as usize) < n && (a_inv as usize) < n { mul[h as usize][a_inv as usize] } else { 255 };
            }
            -1
        })
        .collect();
    let one1: Vec<Alg> = coset.iter().map(|&k| Alg::w_pow(k)).collect();
    let one2: Vec<Alg> = coset.iter().map(|&k| Alg::w_pow(2 * k)).collect();
    // chi_3 = 4 Re^2 - 1 = a^2 - 1 in doubled coordinates
    let three: Vec<Alg> = el.iter().map(|q| Alg::int(q[0] * q[0] - 1)).collect();
    let fund: Vec<Alg> = chi_f.clone();
    let two1: Vec<Alg> = (0..n).map(|g| fund[g].times(one1[g], EISENSTEIN)).collect();
    let two2: Vec<Alg> = (0..n).map(|g| fund[g].times(one2[g], EISENSTEIN)).collect();
    let electric: Vec<u8> = (0..n as u8).filter(|&g| el[g as usize][0] == 1).collect();
    FGroup {
        name: "2T",
        n,
        mul,
        inv,
        electric,
        field: EISENSTEIN,
        chi_f,
        chars: vec![vec![Alg::int(1); n], one1, one2, three, fund, two1, two2],
        dims: vec![1, 1, 1, 3, 2, 2, 2],
        fund: 4,
    }
}

/// Quaternion with doubled coordinates in Z[sqrt2]: q = (X + Y i + Z j + U k) / 2
type Hq2 = [Alg; 4];

/// exact Hamilton product of doubled Z[sqrt2] coordinates: 2 (p q) = (2p)(2q)/2, None if not divisible
fn hmul2(p: &Hq2, q: &Hq2) -> Option<Hq2> {
    let f = SQRT2;
    let m = |a: Alg, b: Alg| a.times(b, f);
    let [a1, b1, c1, d1] = *p;
    let [a2, b2, c2, d2] = *q;
    let r = [
        m(a1, a2).minus(m(b1, b2)).minus(m(c1, c2)).minus(m(d1, d2)),
        m(a1, b2).plus(m(b1, a2)).plus(m(c1, d2)).minus(m(d1, c2)),
        m(a1, c2).minus(m(b1, d2)).plus(m(c1, a2)).plus(m(d1, b2)),
        m(a1, d2).plus(m(b1, c2)).minus(m(c1, b2)).plus(m(d1, a2)),
    ];
    r.iter().all(|v| v.x % 2 == 0 && v.y % 2 == 0).then(|| [0, 1, 2, 3].map(|k| Alg { x: r[k].x / 2, y: r[k].y / 2 }))
}

/// The 48 elements of 2O as doubled Z[sqrt2] coordinates: the 24 of 2T in the t2_elements order, then
/// (+-u +-v)/sqrt2 for the pairs (u, v) = (1,i), (1,j), (1,k), (i,j), (i,k), (j,k) with signs ++, +-, -+, --.
pub fn o2_elements() -> Vec<Hq2> {
    let mut v: Vec<Hq2> = t2_elements().iter().map(|q| q.map(Alg::int)).collect();
    for (a, b) in [(0usize, 1usize), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
        for (sa, sb) in [(1i64, 1i64), (1, -1), (-1, 1), (-1, -1)] {
            let mut e = [Alg::int(0); 4];
            e[a] = Alg { x: 0, y: sa };
            e[b] = Alg { x: 0, y: sb };
            v.push(e);
        }
    }
    v
}

/// 2O: multiplication table from exact Z[sqrt2] Hamilton products; chi_f = 2 Re(q) in Z[sqrt2]; characters
/// trivial, sign (+1 on 2T, -1 off it), chi_2 (2 on Q8, -1 on 2T \ Q8, 0 off 2T), chi_3 = 4 Re^2 - 1,
/// chi_3' = chi_3 sign, chi_f, chi_f' = chi_f sign, chi_4 = chi_f chi_2. Electric set: the six elements with
/// Re(q) = 1/sqrt2.
pub fn o2() -> FGroup {
    let el = o2_elements();
    let n = el.len();
    let mut mul = vec![vec![255u8; n]; n];
    for (x, p) in el.iter().enumerate() {
        for (y, q) in el.iter().enumerate() {
            if let Some(r) = hmul2(p, q) {
                mul[x][y] = el.iter().position(|z| *z == r).map_or(255, |v| v as u8);
            }
        }
    }
    let inv: Vec<u8> = (0..n).map(|x| (0..n as u8).find(|&y| mul[x][y as usize] == 0).unwrap_or(255)).collect();
    let f = SQRT2;
    let chi_f: Vec<Alg> = el.iter().map(|q| q[0]).collect();
    let sign: Vec<Alg> = (0..n).map(|g| Alg::int(if g < 24 { 1 } else { -1 })).collect();
    let chi2: Vec<Alg> = (0..n).map(|g| Alg::int(if g < 8 { 2 } else if g < 24 { -1 } else { 0 })).collect();
    let chi3: Vec<Alg> = el.iter().map(|q| q[0].times(q[0], f).minus(Alg::int(1))).collect();
    let chi3p: Vec<Alg> = (0..n).map(|g| chi3[g].times(sign[g], f)).collect();
    let chifp: Vec<Alg> = (0..n).map(|g| chi_f[g].times(sign[g], f)).collect();
    let chi4: Vec<Alg> = (0..n).map(|g| chi_f[g].times(chi2[g], f)).collect();
    let electric: Vec<u8> = (0..n as u8).filter(|&g| el[g as usize][0] == Alg { x: 0, y: 1 }).collect();
    FGroup {
        name: "2O",
        n,
        mul,
        inv,
        electric,
        field: f,
        chi_f: chi_f.clone(),
        chars: vec![vec![Alg::int(1); n], sign, chi2, chi3, chi3p, chi_f, chifp, chi4],
        dims: vec![1, 1, 2, 3, 3, 2, 2, 4],
        fund: 5,
    }
}

pub fn group(id: i64) -> Option<FGroup> {
    match id {
        1 => Some(s3()),
        2 => Some(q8()),
        3 => Some(t2()),
        4 => Some(o2()),
        _ => None,
    }
}

/// Open ROWS x COLS strip: comb spanning tree (column-0 verticals, every horizontal), free links = the other
/// verticals (one per plaquette).
pub struct Strip {
    pub rows: usize,
    pub cols: usize,
    pub edges: Vec<(usize, usize)>,
    pub hid: Vec<usize>,
    pub vid: Vec<usize>,
    pub tree: Vec<usize>,
    pub free: Vec<usize>,
    pub plaq: Vec<[usize; 4]>,
    pub sites: usize,
}

pub fn strip(rows: usize, cols: usize) -> Strip {
    let site = |x: usize, y: usize| y * cols + x;
    let (mut edges, mut hid, mut vid) = (Vec::new(), vec![usize::MAX; rows * cols], vec![usize::MAX; rows * cols]);
    for y in 0..rows {
        for x in 0..cols - 1 {
            hid[site(x, y)] = edges.len();
            edges.push((site(x, y), site(x + 1, y)));
        }
    }
    for y in 0..rows - 1 {
        for x in 0..cols {
            vid[site(x, y)] = edges.len();
            edges.push((site(x, y), site(x, y + 1)));
        }
    }
    let mut tree: Vec<usize> = (0..rows - 1).map(|y| vid[site(0, y)]).collect();
    for y in 0..rows {
        for x in 0..cols - 1 {
            tree.push(hid[site(x, y)]);
        }
    }
    let mut free = Vec::new();
    for y in 0..rows - 1 {
        for x in 1..cols {
            free.push(vid[site(x, y)]);
        }
    }
    let mut plaq = Vec::new();
    for y in 0..rows - 1 {
        for x in 0..cols - 1 {
            plaq.push([hid[site(x, y)], vid[site(x + 1, y)], hid[site(x, y + 1)], vid[site(x, y)]]);
        }
    }
    Strip { rows, cols, edges, hid, vid, tree, free, plaq, sites: rows * cols }
}

pub struct Sys<'a> {
    pub g: &'a FGroup,
    pub s: &'a Strip,
}

impl Sys<'_> {
    pub fn raw(&self, q: &[u8]) -> Vec<u8> {
        let mut l = vec![0u8; self.s.edges.len()];
        for (i, &e) in self.s.free.iter().enumerate() {
            l[e] = q[i];
        }
        l
    }
    /// gauge-fix a raw configuration on the tree (root (0,0)): free-link values h_t U h_h^-1
    pub fn fix(&self, l: &[u8]) -> Vec<u8> {
        let mut h = vec![0u8; self.s.sites];
        for &e in &self.s.tree {
            let (t, hd) = self.s.edges[e];
            h[hd] = self.g.m(h[t], l[e]);
        }
        self.s.free.iter().map(|&e| self.g.m(self.g.m(h[self.s.edges[e].0], l[e]), self.g.i(h[self.s.edges[e].1]))).collect()
    }
    pub fn code(&self, q: &[u8]) -> u64 {
        let n = self.g.n as u64;
        q.iter().rev().fold(0u64, |a, &v| a * n + u64::from(v))
    }
    pub fn conjugates(&self, q: &[u8]) -> Vec<u64> {
        (0..self.g.n as u8).map(|c| self.code(&q.iter().map(|&x| self.g.m(self.g.m(c, x), self.g.i(c))).collect::<Vec<u8>>())).collect()
    }
    pub fn key(&self, q: &[u8]) -> u64 {
        self.conjugates(q).into_iter().min().unwrap_or(0)
    }
    pub fn magnetic(&self, l: &[u8]) -> Alg {
        self.s.plaq.iter().fold(Alg::int(0), |acc, p| acc.plus(Alg::int(2).minus(self.g.chi_f[self.g.m(self.g.m(self.g.m(l[p[0]], l[p[1]]), self.g.i(l[p[2]])), self.g.i(l[p[3]])) as usize])))
    }
    /// chi_f of the rectangle on rows 0..1, columns x0..x1, counter-clockwise from (x0, 0)
    pub fn wilson(&self, l: &[u8], x0: usize, x1: usize) -> Alg {
        let s = self.s;
        let g = self.g;
        let mut u = 0u8;
        for x in x0..x1 {
            u = g.m(u, l[s.hid[x]]);
        }
        u = g.m(u, l[s.vid[x1]]);
        for x in (x0..x1).rev() {
            u = g.m(u, g.i(l[s.hid[s.cols + x]]));
        }
        u = g.m(u, g.i(l[s.vid[x0]]));
        g.chi_f[u as usize]
    }
    /// local gauge transformation by h at vertex v: U_e -> h U_e on links leaving v, U_e -> U_e h^-1 on
    /// links entering v
    pub fn gauge(&self, l: &[u8], v: usize, h: u8) -> Vec<u8> {
        let mut out = l.to_vec();
        for (e, &(t, hd)) in self.s.edges.iter().enumerate() {
            if t == v {
                out[e] = self.g.m(h, out[e]);
            }
            if hd == v {
                out[e] = self.g.m(out[e], self.g.i(h));
            }
        }
        out
    }
}

/// Basis: representatives, weights, index by code (orbit keys or configuration codes), electric rows.
pub struct Basis {
    pub reps: Vec<Vec<u8>>,
    pub w: Vec<i64>,
    pub index: HashMap<u64, usize>,
    pub e: Vec<BTreeMap<usize, i64>>,
    pub orbit: bool,
}

pub fn basis(sys: &Sys, orbit: bool, electric_links: &[usize]) -> Basis {
    let g = sys.g;
    let nf = sys.s.free.len();
    let n = g.n as u64;
    let total = n.pow(nf as u32);
    let mut reps = Vec::new();
    let mut w = Vec::new();
    for c in 0..total {
        let q: Vec<u8> = (0..nf).map(|i| ((c / n.pow(i as u32)) % n) as u8).collect();
        if orbit {
            let cs = sys.conjugates(&q);
            if cs.iter().min() == Some(&c) {
                let mut s = cs.clone();
                s.sort_unstable();
                s.dedup();
                w.push(s.len() as i64);
                reps.push(q);
            }
        } else {
            w.push(1);
            reps.push(q);
        }
    }
    let index: HashMap<u64, usize> = reps.iter().enumerate().map(|(i, q)| (sys.code(q), i)).collect();
    let e = move_matrix(sys, &reps, &index, orbit, electric_links);
    Basis { reps, w, index, e, orbit }
}

/// E restricted to the given links: rows {j: E_ij} with E_ij = |S| |links| delta_ij - #moves i -> j
pub fn move_matrix(sys: &Sys, reps: &[Vec<u8>], index: &HashMap<u64, usize>, orbit: bool, links: &[usize]) -> Vec<BTreeMap<usize, i64>> {
    let g = sys.g;
    reps.iter()
        .map(|q| {
            let raw = sys.raw(q);
            let mut row: BTreeMap<usize, i64> = BTreeMap::new();
            *row.entry(index[&sys.code(q)]).or_insert(0) += (g.electric.len() * links.len()) as i64;
            for &e in links {
                for &t in &g.electric {
                    let mut l = raw.clone();
                    l[e] = g.m(t, l[e]);
                    let f = sys.fix(&l);
                    let k = if orbit { sys.key(&f) } else { sys.code(&f) };
                    *row.entry(index[&k]).or_insert(0) -= 1;
                }
            }
            row.retain(|_, v| *v != 0);
            row
        })
        .collect()
}

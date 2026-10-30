//! fgroup: finite-group lattice gauge lane for the SU(N) programme. Group data (S3 for regression against
//! the sealed S3 packages, Q8 from Gaussian-integer matrices, 2T from Hurwitz quaternions with Z[w]
//! characters), the open ROWS x COLS strip with a comb gauge
//! tree, gauge-orbit and configuration bases, the electric operator E = sum over links and g in S of
//! (1 - L_g), the magnetic diagonal M = sum over plaquettes (2 - chi_f(U_p)) and rectangular Wilson loops.
//! Integer arithmetic only.
#![forbid(unsafe_code)]
#![deny(clippy::float_arithmetic, clippy::float_cmp, clippy::cast_precision_loss)]

use std::collections::{BTreeMap, HashMap};

/// Eisenstein integer x + y w, w^2 + w + 1 = 0 (exact pairs; integers have y = 0).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Zw {
    pub x: i64,
    pub y: i64,
}

impl Zw {
    pub const fn int(x: i64) -> Zw {
        Zw { x, y: 0 }
    }
    pub fn plus(self, o: Zw) -> Zw {
        Zw { x: self.x + o.x, y: self.y + o.y }
    }
    /// (a + b w)(c + d w) = (ac - bd) + (ad + bc - bd) w
    pub fn times(self, o: Zw) -> Zw {
        Zw { x: self.x * o.x - self.y * o.y, y: self.x * o.y + self.y * o.x - self.y * o.y }
    }
    /// complex conjugate: w -> w^2 = -1 - w
    pub fn conj(self) -> Zw {
        Zw { x: self.x - self.y, y: -self.y }
    }
    /// w^k
    pub fn w_pow(k: i64) -> Zw {
        match k.rem_euclid(3) {
            0 => Zw::int(1),
            1 => Zw { x: 0, y: 1 },
            _ => Zw { x: -1, y: -1 },
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
    /// fundamental (Wilson-loop and plaquette) character, real and integer valued
    pub chi_f: Vec<i64>,
    /// irreducible characters (rows), Eisenstein-integer valued, and dimensions
    pub chars: Vec<Vec<Zw>>,
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
    /// electric eigenvalue of irrep k: |S| - sum_{g in S} chi_k(g) / d_k (exact; the sum must be a rational
    /// integer, which holds because S is inverse-closed, and divisible by d_k)
    pub fn e_irrep(&self, k: usize) -> Option<i64> {
        let s = self.electric.iter().fold(Zw::int(0), |a, &g| a.plus(self.chars[k][g as usize]));
        let d = self.dims[k];
        (s.y == 0 && s.x % d == 0).then(|| self.electric.len() as i64 - s.x / d)
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
        chi_f: (0..6).map(std).collect(),
        chars: vec![vec![Zw::int(1); 6], (0..6).map(|g| Zw::int(sign(g))).collect(), (0..6).map(|g| Zw::int(std(g))).collect()],
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
    let one_dim = |keep: usize| -> Vec<Zw> { (0..n).map(|g| Zw::int(if g / 2 == 0 || g / 2 == keep { 1 } else { -1 })).collect() };
    FGroup {
        name: "Q8",
        n,
        mul,
        inv,
        electric: vec![2, 3, 4, 5, 6, 7],
        chi_f: tr.clone(),
        chars: vec![vec![Zw::int(1); n], one_dim(1), one_dim(2), one_dim(3), tr.iter().map(|&v| Zw::int(v)).collect()],
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
    let chi_f: Vec<i64> = el.iter().map(|q| q[0]).collect();
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
    let one1: Vec<Zw> = coset.iter().map(|&k| Zw::w_pow(k)).collect();
    let one2: Vec<Zw> = coset.iter().map(|&k| Zw::w_pow(2 * k)).collect();
    // chi_3 = 4 Re^2 - 1 = a^2 - 1 in doubled coordinates
    let three: Vec<Zw> = el.iter().map(|q| Zw::int(q[0] * q[0] - 1)).collect();
    let fund: Vec<Zw> = chi_f.iter().map(|&v| Zw::int(v)).collect();
    let two1: Vec<Zw> = (0..n).map(|g| fund[g].times(one1[g])).collect();
    let two2: Vec<Zw> = (0..n).map(|g| fund[g].times(one2[g])).collect();
    let electric: Vec<u8> = (0..n as u8).filter(|&g| el[g as usize][0] == 1).collect();
    FGroup {
        name: "2T",
        n,
        mul,
        inv,
        electric,
        chi_f,
        chars: vec![vec![Zw::int(1); n], one1, one2, three, fund, two1, two2],
        dims: vec![1, 1, 1, 3, 2, 2, 2],
        fund: 4,
    }
}

pub fn group(id: i64) -> Option<FGroup> {
    match id {
        1 => Some(s3()),
        2 => Some(q8()),
        3 => Some(t2()),
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
    pub fn magnetic(&self, l: &[u8]) -> i64 {
        self.s.plaq.iter().map(|p| 2 - self.g.chi_f[self.g.m(self.g.m(self.g.m(l[p[0]], l[p[1]]), self.g.i(l[p[2]])), self.g.i(l[p[3]])) as usize]).sum()
    }
    /// chi_f of the rectangle on rows 0..1, columns x0..x1, counter-clockwise from (x0, 0)
    pub fn wilson(&self, l: &[u8], x0: usize, x1: usize) -> i64 {
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

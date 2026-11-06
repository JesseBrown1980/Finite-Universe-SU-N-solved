#!/usr/bin/env python3
"""Exact integer 2x2x2 x Nt 4D slice-transfer TEST for Z3 and S3 <= SU(3).

Only the unprojected entry kernel is exact for the complete gauge model.
The small four-state boundary catalogue is a non-closed *diagnostic projection*;
its trace, gap or spectrum is not the full finite-universe gauge spectrum.
No floating-point, continuum assumption, JSON, or external dependency.
"""
import hashlib
from itertools import permutations, product
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
VERTICES = tuple(range(8))
DIRECTIONS = (0, 1, 2)
LINKS = tuple((x, mu) for x in VERTICES for mu in DIRECTIONS)
LINK_POS = {v: i for i, v in enumerate(LINKS)}
PAIR_EDGES = tuple((u, u ^ (1 << axis)) for u in VERTICES for axis in DIRECTIONS if u < (u ^ (1 << axis)))
assert len(LINKS) == 24 and len(PAIR_EDGES) == 12


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


class Group:
    def __init__(self, name):
        self.name = name
        if name == 'Z3':
            self.elements = (0, 1, 2)
            self.identity = 0
            self.mul = tuple(tuple((a + b) % 3 for b in range(3)) for a in range(3))
            self.inv = (0, 2, 1)
            self.weight = (4, 1, 1)
        elif name == 'S3':
            ps = tuple(permutations(range(3)))
            ident = (0, 1, 2)
            # The sign-twisted 3x3 permutation matrices form a faithful subgroup of SU(3).
            where = {p: i for i, p in enumerate(ps)}
            self.elements = tuple(range(6))
            self.identity = where[ident]
            self.mul = tuple(tuple(where[tuple(p[q[i]] for i in range(3))] for q in ps) for p in ps)
            self.inv = tuple(where[tuple(p.index(i) for i in range(3))] for p in ps)
            weights = []
            signs = []
            for p in ps:
                inv_count = sum(int(p[i] > p[j]) for i in range(3) for j in range(i + 1, 3))
                sign = 1 if inv_count % 2 == 0 else -1
                signs.append(sign)
                trace = sign * sum(int(p[i] == i) for i in range(3))
                weights.append(2 + trace)
            self.weight = tuple(weights)
            self.permutation = ps
            self.signs = tuple(signs)
            assert sorted(self.weight) == [1,1,1,2,2,5]
        else:
            raise ValueError('unknown finite group')
        self.n = len(self.elements)
        for g in self.elements:
            assert self.mul[self.identity][g] == g == self.mul[g][self.identity]
            assert self.mul[g][self.inv[g]] == self.identity
            assert self.mul[self.inv[g]][g] == self.identity
            for h in self.elements:
                assert self.weight[self.mul[h][self.mul[g][self.inv[h]]]] == self.weight[g]

    def word(self, *vals):
        out = self.identity
        for v in vals:
            out = self.mul[out][v]
        return out


def cold(group):
    return tuple(group.identity for unused in LINKS)


def shifted(x, mu):
    return x ^ (1 << mu)


def magnetic_holonomies(group, spatial):
    assert len(spatial) == 24
    for x in VERTICES:
        for mu in DIRECTIONS:
            for nu in range(mu + 1, 3):
                yield group.word(spatial[LINK_POS[x, mu]],
                                 spatial[LINK_POS[shifted(x, mu), nu]],
                                 group.inv[spatial[LINK_POS[shifted(x, nu), mu]]],
                                 group.inv[spatial[LINK_POS[x, nu]]])


def magnetic_weight(group, spatial):
    return product_int(group.weight[g] for g in magnetic_holonomies(group, spatial))


def product_int(vals):
    a = 1
    for b in vals:
        a *= b
    return a


def temporal_holonomy(group, before, after, tlinks, x, mu):
    y = shifted(x, mu)
    return group.word(before[LINK_POS[x, mu]],tlinks[y],group.inv[after[LINK_POS[x, mu]]],group.inv[tlinks[x]])


def temporal_plaquette_weight(group, before, after, tlinks):
    return product_int(group.weight[temporal_holonomy(group,before,after,tlinks,x,mu)] for x, mu in LINKS)


def independent_gauge_transform(group, spatial, h):
    assert len(h) == 8
    return tuple(group.word(h[x], spatial[LINK_POS[x, mu]],group.inv[h[shifted(x, mu)]]) for x,mu in LINKS)


def pair_tables(group, before, after):
    # Time-like plaquettes on the two directed links of every undirected cube edge.
    q = group.n
    out = []
    for u, v in PAIR_EDGES:
        mu = ((u ^ v).bit_length() - 1)
        tab = []
        for gu in range(q):
            row = []
            for gv in range(q):
                a = group.word(before[LINK_POS[u, mu]],gv,group.inv[after[LINK_POS[u, mu]]],group.inv[gu])
                b = group.word(before[LINK_POS[v, mu]],gu,group.inv[after[LINK_POS[v, mu]]],group.inv[gv])
                row.append(group.weight[a]*group.weight[b])
            tab.append(tuple(row))
        out.append(((u,v),tuple(tab)))
    return tuple(out)


def factor_contract(group, before, after):
    """Variable elimination of 8 temporal links; exact sum over G^8."""
    q = group.n
    factors = [(scope,{(i,j): tab[i][j] for i in range(q) for j in range(q)})
               for scope,tab in pair_tables(group,before,after)]
    candidates = set(VERTICES)
    max_scope = 2
    while candidates:
        # Choose a variable with the smallest next cluster; deterministic tie-break.
        choices = []
        for var in sorted(candidates):
            touched = [f for f in factors if var in f[0]]
            scope = tuple(sorted(set(a for terms, _ in touched for a in terms) - {var}))
            choices.append((len(scope),var,scope,touched))
        arity,var,scope,touched = min(choices,key=lambda x:(x[0],x[1]))
        max_scope = max(max_scope,arity+1)
        residual = [f for f in factors if var not in f[0]]
        entries = {}
        for other in product(range(q),repeat=len(scope)):
            assign = dict(zip(scope,other))
            total = 0
            for g in range(q):
                assign[var] = g
                term = 1
                for fs,tab in touched:
                    term *= tab[tuple(assign[x] for x in fs)]
                total += term
            entries[other] = total
        residual.append((scope,entries))
        factors = residual
        candidates.remove(var)
    assert all(not vars for vars,unused in factors)
    return product_int(tab[()] for unused,tab in factors),max_scope


def reference_bruteforce(group, before, after):
    """Independent sum over all 3^8 or 6^8 temporal-link assignments.

    Uses the exact plaquette pair table but *not* variable elimination. Each pair
    entry is separately checked against the direct four-dimensional holonomy.
    """
    tables = pair_tables(group,before,after)
    total = 0
    for g in product(group.elements,repeat=8):
        term = 1
        for (u,v),tab in tables:
            term *= tab[g[u]][g[v]]
        total += term
    return total


def raw_time_entry(group, before, after, tlinks):
    tabprod = 1
    for (u,v),tab in pair_tables(group,before,after):
        tabprod *= tab[tlinks[u]][tlinks[v]]
    return tabprod,temporal_plaquette_weight(group,before,after,tlinks)


def catalogue(group):
    e=group.identity
    vals=[cold(group)]
    a=next(z for z in group.elements if z!=e)
    if group.name=='S3':
        a=next(z for z in group.elements if group.weight[z]==1)
        b=next(z for z in group.elements if group.weight[z]==2)
        assert group.mul[a][b]!=group.mul[b][a]
    else:
        b=2
    def changes(*edits):
        s=list(cold(group))
        for x,mu,value in edits:
            s[LINK_POS[x,mu]]=value
        return tuple(s)
    vals.append(changes((0,0,a)))
    vals.append(changes((0,0,a),(0,1,b)))
    vals.append(changes((0,0,a),(0,1,b),(1,2,a)))
    return tuple(vals)


def dense_mul(a,b):
    n=len(a)
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)) for i in range(n))


def eye(n):
    return tuple(tuple(int(i==j) for j in range(n)) for i in range(n))


def power(a,n):
    assert n>=0
    z=eye(len(a));cur=a
    while n:
        if n&1:
            z=dense_mul(z,cur)
        n//=2
        if n:
            cur=dense_mul(cur,cur)
    return z


def matvec(a,vec):
    return tuple(sum(row[j]*vec[j] for j in range(len(vec))) for row in a)


def vec_bytes(v):
    return (','.join(str(x) for x in v)+'\n').encode('ascii')


def trace(a):
    return sum(row[i] for i,row in enumerate(a))


def chain(group,tm,steps=27):
    vec=(1,0,0,0)
    parent='0'*64
    chainlines=[]
    leaves=[]
    for t in range(steps+1):
        leaf=sha(vec_bytes(vec))
        row=f'SLICE|group={group.name}|t={t}|dim=4|state_sha256={leaf}|previous={parent}|json=0'
        parent=sha((row+'\n').encode('ascii'))
        chainlines.append(row)
        leaves.append(leaf)
        if t<steps:
            vec=matvec(tm,vec)
    assert vec==matvec(power(tm,steps),(1,0,0,0))
    for t in [3,6,9,27]:
        if t<=steps:
            a=(1,0,0,0)
            for unused in range(t):
                a=matvec(tm,a)
            assert matvec(power(tm,t),(1,0,0,0))==a
    return chainlines,leaves,parent,vec


def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--out',default=str(ROOT))
    args=p.parse_args()
    dest=Path(args.out)
    dest.mkdir(parents=True,exist_ok=True)
    rows=[];chainrows=[];rootrows=[]; tm_bytes=[]
    for group_name in ('Z3','S3'):
        g=Group(group_name)
        cold_boundary=cold(g)
        C=catalogue(g)
        assert len(tuple(magnetic_holonomies(g,cold_boundary)))==24
        expected_cold=(4 if group_name=='Z3' else 5)**24
        assert magnetic_weight(g,cold_boundary)==expected_cold
        if group_name=='S3':
            a=next(z for z in g.elements if g.weight[z]==1)
            b=next(z for z in g.elements if g.weight[z]==2)
            ab=g.mul[a][b];ba=g.mul[b][a]
            assert ab!=ba
            assert g.word(a,b,g.inv[a],g.inv[b])!=g.identity
        rows.append(f'MODEL|group={group_name}|sites=8|spatial_links=24|spatial_plaquettes=24|temporal_plaquettes=24|weight_cold={expected_cold}|json=0')
        magnets=[magnetic_weight(g,c) for c in C]
        kernels={};scopes=[]
        for i,u in enumerate(C):
            for j,v in enumerate(C):
                value,arity=factor_contract(g,u,v)
                kernels[i,j]=value
                scopes.append(arity)
        for i in range(4):
            for j in range(4):
                assert kernels[i,j]==kernels[j,i]
        tm=tuple(tuple(magnets[i]*kernels[i,j] for j in range(4)) for i in range(4))
        rows.append(f'CONTRACTION|group={group_name}|boundary_catalog=4|entries=16|temporal_variables=8|max_scope={max(scopes)}|cold_K={kernels[0,0]}|cold_T={tm[0][0]}|json=0')
        for idx in ((0,0),(0,2)):
            i,j=idx
            brute=reference_bruteforce(g,C[i],C[j])
            assert brute==kernels[i,j]
            rows.append(f'BRUTE|group={group_name}|pair={i},{j}|assignment_count={g.n**8}|match=1|json=0')
        for i,j in ((0,0),(0,2),(2,1)):
            for k in range(32):
                tlinks=tuple((k*3+x*(k+1)+x*x)%g.n for x in VERTICES)
                a,b=raw_time_entry(g,C[i],C[j],tlinks)
                assert a==b
        rows.append(f'PLAQUETTE|group={group_name}|direct_vs_edge_pairs=96|match=1|json=0')
        # Independent gauge changes for before/after, with no commuting assumption.
        h=tuple((x*3+1)%g.n for x in VERTICES)
        k=tuple((x*5+2)%g.n for x in VERTICES)
        for i,j in ((0,0),(2,3)):
            before=independent_gauge_transform(g,C[i],h)
            after=independent_gauge_transform(g,C[j],k)
            assert magnetic_weight(g,before)==magnets[i]
            assert factor_contract(g,before,after)[0]==kernels[i,j]
        rows.append(f'GAUGE|group={group_name}|independent_before_after_transform=2|pass=1|json=0')
        for size in (2,3,6):
            dense=power(tm,size)
            direct=eye(4)
            for unused in range(size):
                direct=dense_mul(tm,direct)
            assert dense==direct
            rows.append(f'CATALOG_TRACE|group={group_name}|Nt={size}|trace={trace(dense)}|projected_not_full=1|json=0')
        ch,leaves,root,vec=chain(g,tm)
        chainrows.extend(ch)
        rootrows.append(f'CHAIN_ROOT|group={group_name}|steps=27|root={root}|end_state_sha256={leaves[-1]}|json=0')
        rows.append(f'JUMP|group={group_name}|steps=27|all_sequential_hashes=28|direct_vs_power=1|root={root}|json=0')
        tm_bytes.extend(f'MATRIX|group={group_name}|i={i}|j={j}|T={tm[i][j]}|json=0' for i in range(4) for j in range(4))
        for test_name,c in [('cold',C[0]),('twisted',C[2])]:
            rows.append(f'WILSON|group={group_name}|configuration={test_name}|magnetic_weight={magnetic_weight(g,c)}|json=0')
    files={'RESULT.hbp':'\n'.join(rows)+'\n','SLICE-CHAIN.hbp':'\n'.join(chainrows)+'\n','ROOT.hbp':'\n'.join(rootrows)+'\n','MATRIX.hbp':'\n'.join(tm_bytes)+'\n'}
    for f,s in files.items():
        (dest/f).write_text(s,encoding='ascii')
    print('PASS: '+str(len(rows))+' result rows, '+str(len(chainrows))+' linked steps')
    for s in rows:
        if s.startswith(('CONTRACTION','BRUTE','GAUGE','JUMP')):
            print(s)


if __name__=='__main__':
    main()

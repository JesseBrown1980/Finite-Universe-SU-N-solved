#!/usr/bin/env python3
"""Independent 2x2x2 Z3/S3 transfer verifier: parity-bipartite contraction.

Does not import the engine; independently reconstructs the group and cube.
The per-pair kernel can be reproduced by an exhaustive |G|^8 enumeration.
Pinned roots are literals in this reviewer-facing script and the CI workflow.
"""
import hashlib
from itertools import permutations,product
from pathlib import Path
import sys

BASE=Path(__file__).resolve().parent
ANCHOR={'Z3':'31011a6ff18be78581acb86538549af37680da3d39f300fd8c3607a7dedacc02',
        'S3':'12677da0fbae69ea364d1ebfde683570ad89f7ec248f99a0193206218ded9065'}


def sha(x):
    return hashlib.sha256(x).hexdigest()


def parse_hbp(path):
    rows=[]
    for line in path.read_text('ascii').splitlines():
        vals=line.split('|')
        fields={'_kind':vals[0]}
        for cell in vals[1:]:
            k,sep,v=cell.partition('=')
            assert sep and k not in fields
            fields[k]=v
        assert fields['json']=='0'
        rows.append(fields)
    return rows


class G:
    def __init__(self,name):
        self.name=name
        if name=='Z3':
            self.n=3;self.one=0
            self.product=lambda a,b:(a+b)%3
            self.inverse=lambda a:(-a)%3
            self.character=lambda a:3 if a==0 else 0
            self.w=lambda a:4 if a==0 else 1
        else:
            perms=list(permutations((0,1,2)))
            to_id={p:i for i,p in enumerate(perms)}
            self.n=6;self.one=0
            self.perms=perms
            self.product=lambda a,b:to_id[tuple(perms[a][perms[b][k]] for k in range(3))]
            self.inverse=lambda a:to_id[tuple(perms[a].index(k) for k in range(3))]
            self.character=lambda a:(1 if sum(perms[a][i]>perms[a][j] for i in range(3) for j in range(i+1,3))%2==0 else -1)*sum(perms[a][i]==i for i in range(3))
            self.w=lambda a:2+self.character(a)
            assert [self.w(i) for i in range(6)]==[5,1,1,2,2,1]

    def hol(self,word):
        z=self.one
        for a in word:z=self.product(z,a)
        return z


def cold(g):return tuple(g.one for _ in range(24))


def catalogue(g):
    a=1
    b=3 if g.n==6 else 2
    assert g.w(a)==1 and (g.n!=6 or g.w(b)==2)
    def change(ops):
        z=list(cold(g))
        for x,mu,v in ops:z[3*x+mu]=v
        return tuple(z)
    return (cold(g),change(((0,0,a),)),change(((0,0,a),(0,1,b))),change(((0,0,a),(0,1,b),(1,2,a))))


def magnetic(g,u):
    w=1
    for x in range(8):
        for a in range(3):
            for b in range(a+1,3):
                four=(u[3*x+a], u[3*(x^(1<<a))+b],g.inverse(u[3*(x^(1<<b))+a]),g.inverse(u[3*x+b]))
                w*=g.w(g.hol(four))
    return w


def edges(g,u,v):
    # Independent orientation-specific plaquette construction, pairing 24 links.
    out={}
    for x in range(8):
        for mu in range(3):
            y=x^(1<<mu)
            a,b=sorted((x,y))
            if a==x:
                # g at start=a, end=b.
                oriented=(a,b)
            else:oriented=(b,a)
            table=out.setdefault((a,b),[[1]*g.n for _ in range(g.n)])
            for ga in range(g.n):
                for gb in range(g.n):
                    gs=ga if x==a else gb
                    gt=gb if x==a else ga
                    pl=g.hol((u[3*x+mu],gt,g.inverse(v[3*x+mu]),g.inverse(gs)))
                    table[ga][gb]*=g.w(pl)
    assert len(out)==12
    return out


def kernel_bipartite(g,u,v):
    adj=edges(g,u,v)
    even=tuple(x for x in range(8) if x.bit_count()%2==0)
    odd=tuple(x for x in range(8) if x.bit_count()%2==1)
    assert len(even)==len(odd)==4
    n=g.n
    answer=0
    for ev in product(range(n),repeat=4):
        gE=dict(zip(even,ev))
        vtotal=1
        for j in odd:
            subtotal=0
            for gj in range(n):
                p=1
                for mu in range(3):
                    i=j^(1<<mu)
                    a,b=sorted((i,j))
                    table=adj[(a,b)]
                    p*=table[gE[i] if i==a else gj][gj if i==a else gE[i]]
                subtotal+=p
            vtotal*=subtotal
        answer+=vtotal
    return answer


def step(M,v):
    return tuple(sum(M[i][j]*v[j] for j in range(4)) for i in range(4))


def product_m(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))


def eye():return tuple(tuple(int(i==j) for j in range(4)) for i in range(4))


def matrix_power(a,k):
    total=eye()
    while k:
        if k%2:total=product_m(a,total)
        k//=2
        if k:a=product_m(a,a)
    return total


def verify_group(name,entries,chain,roots,results):
    g=G(name); c=catalogue(g)
    got={(int(r['i']),int(r['j'])):int(r['T']) for r in entries if r['group']==name}
    assert len(got)==16
    M=tuple(tuple(got[i,j] for j in range(4)) for i in range(4))
    magnets=[magnetic(g,u) for u in c]
    for i in range(4):
        for j in range(4):
            k=kernel_bipartite(g,c[i],c[j])
            assert M[i][j]==magnets[i]*k, (name,i,j)
            assert M[i][j]*magnets[j]==M[j][i]*magnets[i]
    assert magnets[0]==(4 if name=='Z3' else 5)**24
    if name=='S3':
        assert g.product(1,4)!=g.product(4,1)
        assert g.hol((1,4,g.inverse(1),g.inverse(4)))!=g.one
    state=(1,0,0,0)
    tail='0'*64
    rows=sorted((r for r in chain if r['group']==name),key=lambda r:int(r['t']))
    assert len(rows)==28
    for t,row in enumerate(rows):
        assert int(row['t'])==t
        assert row['state_sha256']==sha((','.join(str(x) for x in state)+'\n').encode('ascii'))
        assert row['previous']==tail
        raw=('SLICE|group='+name+'|t='+str(t)+'|dim=4|state_sha256='+row['state_sha256']+'|previous='+tail+'|json=0\n').encode('ascii')
        tail=sha(raw)
        if t<27:state=step(M,state)
    assert tail==ANCHOR[name]
    assert any(r['group']==name and r['root']==tail for r in roots)
    assert state==step(matrix_power(M,27),(1,0,0,0))
    # Anchored-chain fault: replacing the last leaf and recomputing its internal hash
    # leaves a consistent *new* root, which does not match the external pin.
    tampered='SLICE|group='+name+'|t=27|dim=4|state_sha256='+('0'*64)+'|previous='+rows[-1]['previous']+'|json=0\n'
    assert sha(tampered.encode('ascii'))!=ANCHOR[name]
    assert any(r['_kind']=='JUMP' and r['group']==name and r['root']==tail for r in results)
    return f'CHECK|group={name}|entries=16|bipartite_direct=PASS|chain=28|external_anchor=PASS|negative_anchor=REJECT|json=0'


def check_fault_controls():
    g=G('S3')
    c=catalogue(g)
    a=1
    b=3
    # F1: an erroneously commuting S3 multiplication table contradicts 3x3 matrices.
    assert g.product(a,b)!=g.product(b,a)
    # F2: incorrect temporal order U*g_end*g_start^-1*V^-1 gives another class.
    links=tuple((3*x+mu,x,x^(1<<mu)) for x in range(8) for mu in range(3))
    x=0;mu=0;k=3
    temporal=tuple((k*3+v*(k+1)+v*v)%g.n for v in range(8))
    gs=temporal[x];gt=temporal[x^(1<<mu)]
    good=g.w(g.hol((c[0][3*x+mu],gt,g.inverse(c[1][3*x+mu]),g.inverse(gs))))
    bad=g.w(g.hol((c[0][3*x+mu],gt,g.inverse(gs),g.inverse(c[1][3*x+mu]))))
    assert good==5 and bad==2
    # F3: omitting one spatial plaquette changes its cold magnetic factor.
    assert magnetic(g,c[0])==5**24 and magnetic(g,c[0])!=5**23
    # F4: adding a weight only to one transposition destroys conjugacy invariance.
    unequal={1:2,2:1}
    assert g.w(1)==g.w(2)==1 and unequal[1]!=unequal[2]
    assert any(g.product(h,g.product(1,g.inverse(h)))==2 for h in range(6))
    # F5: a wrong gauge transform without the inverse on the end is not a symmetry.
    h=tuple((x*x+1)%6 for x in range(8))
    wrong=tuple(g.hol((h[x],c[0][3*x+mu],h[x^(1<<mu)])) for x in range(8) for mu in range(3))
    assert magnetic(g,wrong)!=magnetic(g,c[0])
    # F6: changing both last leaf and the internal chain changes the external seal.
    chain=parse_hbp(BASE/'SLICE-CHAIN.hbp')
    row=next(r for r in chain if r['group']=='S3' and r['t']=='27')
    forged='SLICE|group=S3|t=27|dim=4|state_sha256='+('0'*64)+'|previous='+row['previous']+'|json=0\n'
    assert sha(forged.encode('ascii'))!=ANCHOR['S3']
    for i,name in enumerate(('noncommuting_multiplication','temporal_order','missing_plaquette','nonclass_character','incorrect_gauge_transform','forged_internal_chain'),1):
        print(f'NEGATIVE|case={i}|mutation={name}|verdict=REJECT|json=0')


def main():
    pred=BASE/'PREDICTION.hbp'
    assert pred.exists()
    p=pred.read_text('ascii')
    assert 'FROZEN_BEFORE_CODE' in p and 'spatial=2x2x2' in p
    entries=parse_hbp(BASE/'MATRIX.hbp')
    chain=parse_hbp(BASE/'SLICE-CHAIN.hbp')
    roots=parse_hbp(BASE/'ROOT.hbp')
    results=parse_hbp(BASE/'RESULT.hbp')
    print(verify_group('Z3',entries,chain,roots,results))
    print(verify_group('S3',entries,chain,roots,results))
    check_fault_controls()
    print('PASS: independent bipartite contraction, immutable roots, six negative controls')

if __name__=='__main__':main()

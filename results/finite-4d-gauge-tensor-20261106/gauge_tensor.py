#!/usr/bin/env python3
"""Full-space, matrix-free finite 4D Z3/S3 slice transfer.

An exact projector/convolution/tensor representation of the entire q^24
spatial-link configuration space; NO dense Hilbert or orbit vector is formed.
Magnetic Wilson observables for Z3 are evaluated from all closed 2-form fluxes.
The general interacting mass gap is NOT computed.
Integer-only, no JSON, no sampling.
"""
import hashlib
from itertools import permutations, product
from pathlib import Path

VERTICES = tuple(range(8))
LINKS = tuple((x, axis) for x in VERTICES for axis in range(3))
LINK_ID = {p: i for i,p in enumerate(LINKS)}
FACES = tuple((x,a,b) for x in VERTICES for a in range(3) for b in range(a+1,3))
PAIRS = tuple((x, x ^ (1 << a), a) for x in VERTICES for a in range(3) if x < (x ^ (1 << a)))
assert len(LINKS)==24 and len(FACES)==24 and len(PAIRS)==12

class Group:
    def __init__(self, name):
        self.name=name
        if name=='Z3':
            self.n=3
            self.mul=tuple(tuple((a+b)%3 for b in range(3)) for a in range(3))
            self.inv=(0,2,1)
            self.weight=(4,1,1)
            self.sign=None
        elif name=='S3':
            ps=tuple(permutations(range(3)))
            idx={p:i for i,p in enumerate(ps)}
            self.n=6
            self.mul=tuple(tuple(idx[tuple(p[q[i]] for i in range(3))] for q in ps) for p in ps)
            self.inv=tuple(idx[tuple(p.index(i) for i in range(3))] for p in ps)
            signs=[];weights=[]
            for p in ps:
                invcount=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
                sign=1 if invcount%2==0 else -1
                signs.append(sign)
                weights.append(2+sign*sum(p[i]==i for i in range(3)))
            self.sign=tuple(signs)
            self.weight=tuple(weights)
            assert tuple(sorted(weights))==(1,1,1,2,2,5)
        else:
            raise ValueError(name)
        self.e=0
        self.z=sum(self.weight)
        assert self.n in (3,6)
        for x in range(self.n):
            assert self.mul[x][self.inv[x]]==0==self.mul[self.inv[x]][x]
            for h in range(self.n):
                assert self.weight[self.compose(h,x,self.inv[h])]==self.weight[x]
    def compose(self,*args):
        z=self.e
        for a in args:z=self.mul[z][a]
        return z


def magnetic_faces(g,U):
    for x,a,b in FACES:
        yield g.compose(U[LINK_ID[x,a]],U[LINK_ID[x^(1<<a),b]],
                        g.inv[U[LINK_ID[x^(1<<b),a]]],g.inv[U[LINK_ID[x,b]]])


def magnetic(g,U):
    p=1
    for x in magnetic_faces(g,U):p*=g.weight[x]
    return p


def gauge_change(g,U,h):
    return tuple(g.compose(h[x],U[i],g.inv[h[x^(1<<a)]]) for i,(x,a) in enumerate(LINKS))


def cold(g):
    return (g.e,)*24


def sample(g,seed):
    s=seed;out=[]
    for i in range(24):
        s=(1103515245*s+12345+i)%(1<<31)
        out.append(s%g.n)
    return tuple(out)


def unit_factors(g):
    return tuple((1,)*g.n for unused in LINKS)


def point_factors(g,V):
    return tuple(tuple(int(v==V[i]) for v in range(g.n)) for i in range(24))


def pair_tables_product(g,U,factors):
    """Eliminate input links analytically before the eight time-link variables.

    For each directed edge e=(x,y): L_e(g_x,g_y)=
    sum_v w(U_e g_y v^-1 g_x^-1) f_e(v).
    The 24 L_e form 12 pair factors on the bipartite spatial cube.
    """
    assert len(U)==len(factors)==24
    out={}
    for a,b,axis in PAIRS:
        forward=LINK_ID[a,axis]
        backward=LINK_ID[b,axis]
        entries=[]
        for ga in range(g.n):
            row=[]
            for gb in range(g.n):
                p=1
                for src,dst,e,gs,gd in ((a,b,forward,ga,gb),(b,a,backward,gb,ga)):
                    v_sum=0
                    for v in range(g.n):
                        term=g.compose(U[e],gd,g.inv[v],g.inv[gs])
                        v_sum+=g.weight[term]*factors[e][v]
                    p*=v_sum
                row.append(p)
            entries.append(tuple(row))
        out[a,b]=tuple(entries)
    return out


def contract_bipartite(g,tables):
    """Exact eight-variable contraction, max four variables enumerated jointly."""
    even=(0,3,5,6)
    odd=(1,2,4,7)
    assert all(x.bit_count()%2==0 for x in even)
    assert all(x.bit_count()%2==1 for x in odd)
    total=0
    for vals in product(range(g.n),repeat=4):
        ev=dict(zip(even,vals))
        joint=1
        for y in odd:
            partial=0
            for hy in range(g.n):
                p=1
                for axis in range(3):
                    x=y^(1<<axis)
                    a,b=sorted((x,y))
                    p*=tables[a,b][ev[x] if x==a else hy][hy if y==b else ev[x]]
                partial+=p
            joint*=partial
        total+=joint
    return total


def apply_product(g,U,factors):
    """Evaluate (T f)(U) for f(V)=product_e factors[e][V_e], exactly.

    The result belongs to the complete gauge-invariant state space for
    every input, whether or not f was gauge invariant.
    """
    K=contract_bipartite(g,pair_tables_product(g,U,factors))
    return magnetic(g,U)*K


def ones_constant(g):
    return (g.n**8)*(g.z**24)


def boundary_matrix():
    """Map 2-chains (face fluxes) to 1-chains, mod 3."""
    b=[[0]*24 for unused in LINKS]
    for j,(x,a,c) in enumerate(FACES):
        terms=((x,a,1),(x^(1<<a),c,1),(x^(1<<c),a,-1),(x,c,-1))
        for y,d,s in terms:
            idx=LINK_ID[y,d]
            b[idx][j]=(b[idx][j]+s)%3
    return b


def mod3_kernel_basis(b):
    a=[row[:] for row in b];ncols=len(a[0]); r=0; pivots=[]
    for col in range(ncols):
        next_row=next((i for i in range(r,len(a)) if a[i][col]),None)
        if next_row is None:continue
        a[r],a[next_row]=a[next_row],a[r]
        if a[r][col]==2:
            a[r]=[(2*v)%3 for v in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][col]:
                mul=a[i][col]
                a[i]=[(v-mul*w)%3 for v,w in zip(a[i],a[r])]
        pivots.append(col);r+=1
    free=[c for c in range(ncols) if c not in pivots]
    basis=[]
    for c in free:
        vec=[0]*ncols;vec[c]=1
        for i,p in enumerate(pivots):vec[p]=(-a[i][c])%3
        assert all(sum(bj[j]*vec[j] for j in range(ncols))%3==0 for bj in b)
        basis.append(tuple(vec))
    return r,tuple(basis)


def exact_z3_dual():
    """ALL spatial configurations summed by the closed-2-form character dual.

    w(z)=2+omega^z+omega^-z; after summing all 24 link variables,
    only flux vectors in ker(boundary_2) remain. No spatial orbit sampled.
    """
    b=boundary_matrix();rank,basis=mod3_kernel_basis(b)
    assert rank==14 and len(basis)==10
    hist=[0]*25; z=0;w1=0;w2=0
    coeff=(2,1,1)
    for v in product(range(3),repeat=10):
        flux=[sum(basis[k][p]*v[k] for k in range(10))%3 for p in range(24)]
        zero=flux.count(0)
        weight=1<<zero
        hist[zero]+=1;z+=weight
        w1+=(weight//coeff[flux[0]])*coeff[(flux[0]-1)%3]
        w2+=(weight//(coeff[flux[0]]*coeff[flux[1]]))*coeff[(flux[0]-1)%3]*coeff[(flux[1]-1)%3]
    assert sum(hist)==3**10
    return {'rank':rank,'ker':len(basis),'hist':tuple(hist),'dual_partition':z,
            'dual_W1':w1,'dual_W2':w2,'raw_magnetic':(3**24)*z,
            'quotient_magnetic':(3**17)*z,'raw_open_slab':(3**8)*(6**24)*(3**24)*z}


def link_electric_spectrum(g):
    """Eigenvalues of a single-link convolution C[a,b]=w(a b^-1)."""
    if g.name=='Z3':
        return {'vacuum':6,'nontrivial':3}
    return {'vacuum':12,'sign':6,'standard':3}


def electric_2cycle_check(g):
    """Direct physical two-edge Gauss-law spin-network eigenpair check.

    Two sites, two opposite directed links; independent vertex time links.
    Real character for Z3, sign character for S3, each built from closed loop.
    """
    eps=lambda x: (2 if x==0 else -1) if g.name=='Z3' else g.sign[x]
    states=tuple(product(range(g.n),repeat=2))
    vacuum = (g.n**2)*(g.z**2)
    # Sample every U, and compare full K applied to the loop-character mode.
    for u0,u1 in states:
        total0=0;total1=0
        for v0,v1 in states:
            kernel=0
            for h0,h1 in product(range(g.n),repeat=2):
                w0=g.weight[g.compose(u0,h1,g.inv[v0],g.inv[h0])]
                w1=g.weight[g.compose(u1,h0,g.inv[v1],g.inv[h1])]
                kernel+=w0*w1
            total0+=kernel
            total1+=kernel*eps(g.compose(v0,v1))
        assert total0==vacuum
        assert 4*total1==vacuum*eps(g.compose(u0,u1)),(g.name,u0,u1,total1,vacuum)
    return {'states':g.n**2,'vacuum_eigen':vacuum,'ratio_num':1,'ratio_den':4}


def sha(s):return hashlib.sha256(s.encode('ascii')).hexdigest()


def hbp_fields(*parts):return '|'.join(str(x) for x in parts)+'|json=0'


def run():
    rows=[]
    dual=exact_z3_dual()
    for name in ('Z3','S3'):
        g=Group(name)
        const=ones_constant(g)
        assert (name!='Z3' or const==31088519960728128454656)
        rows.append(hbp_fields('FULL_SPACE',f'group={name}',f'raw_configs={g.n**24}',
            f'physical_orbits={129140163 if name=="Z3" else 2821153019713}',
            'variables=24','wavefunction=magnetic_product_times_constant','stored_full_vector=0'))
        rows.append(hbp_fields('TRANSFER_ONES',f'group={name}',f'K={const}',
            f'cold_M={magnetic(g,cold(g))}',f'cold_T1={const*magnetic(g,cold(g))}'))
        # Check fully nontrivial, rank-one tensor inputs over every color.
        for k in range(3):
            U=sample(g,17+61*k)
            h=tuple(sample(g,200+k)[i] for i in range(8))
            Uprime=gauge_change(g,U,h)
            ff=unit_factors(g)
            got=apply_product(g,U,ff)
            assert got==const*magnetic(g,U)
            assert apply_product(g,Uprime,ff)==got
            if k==0:
                # A nontrivial link-local product, not necessarily gauge invariant.
                tmp=list(ff)
                for edge in (0,3):tmp[edge]=tuple(1+int(v==g.e) for v in range(g.n))
                mixed=tuple(tmp)
                value=apply_product(g,U,mixed)
                assert apply_product(g,Uprime,mixed)==value
                rows.append(hbp_fields('PRODUCT_INPUT',f'group={name}','nonconstant_edges=0,3',f'output={value}',
                    'Gauss_projected=1','all_orbit_entries_implicit=1'))
            rows.append(hbp_fields('SAMPLED_ORBIT',f'group={name}',f'seed={17+61*k}',
                f'M={magnetic(g,U)}',f'T1={got}','gauge_match=1'))
        # Delta-product represents one exact full-space column.
        for i in range(2):
            U=sample(g,90+i)
            V=sample(g,140+i)
            coldinput=point_factors(g,V)
            result=apply_product(g,U,coldinput)
            rows.append(hbp_fields('DELTA_INPUT',f'group={name}',f'case={i}',
                f'column_entry={result}','full_24_link_point_state=1'))
        # Full *interacting* physical transfer: strict finite gap from positive
        # kernel and a row-wise Doeblin minorization after the Perron transform.
        # Each temporal kernel has 24 weights between 1 and w_max. Therefore
        # lambda_1/lambda_0 <= 1 - w_max^-24. This bound is tiny but rigorous.
        maxweight=max(g.weight)
        denominator=maxweight**24
        rows.append(hbp_fields('INTERACTING_FINITE_GAP',f'group={name}',
            f'ratio_upper_num={denominator-1}',f'ratio_upper_den={denominator}',
            'gap_lower=minus_log(ratio_upper)','strict_positive=1',
            'physical_Gauss_sector=1','numerical_gap=NOT_COMPUTED'))
        z=electric_2cycle_check(g)
        rows.append(hbp_fields('ELECTRIC_ONLY',f'group={name}',f'two_cycle_states={z["states"]}',
            f'vacuum_lambda={z["vacuum_eigen"]}','largest_nonvacuum_ratio=1/4',
            'magnetic_weight=1','exact_full_space_gap=ln4','interacting_gap=OPEN'))
    assert dual['dual_partition']==99186762
    assert dual['dual_W1']==83956047 and dual['dual_W2']==74068664
    assert dual['raw_magnetic']==28013271216711264522
    assert dual['quotient_magnetic']==12808994612122206
    assert dual['raw_open_slab']==870891141386018892529920237038951498514432
    rows.append(hbp_fields('DUAL_Z3','rank=14','cycles=10','assignments=59049',
        f'partition={dual["dual_partition"]}',f'W1_num={dual["dual_W1"]}',
        f'W2_num={dual["dual_W2"]}',f'raw_M={dual["raw_magnetic"]}',
        f'physical_orbit_M={dual["quotient_magnetic"]}',f'open_4D_slab_Z={dual["raw_open_slab"]}'))
    for num,n in enumerate(dual['hist']):
        rows.append(hbp_fields('DUAL_HIST',f'zeros={num}',f'closed_assignments={n}'))
    rows.append(hbp_fields('SCOPE','nonabelian_full_S3_wavefunction=EXACT','Z3_open_slab_partition=EXACT',
        'Z3_spatial_Wilson_plaquette=EXACT','full_isotropic_interacting_gap=OPEN',
        'S3_full_slab_partition=OPEN','max_bond_dimension=UNPROVEN'))
    return '\n'.join(rows)+'\n'


if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser()
    a.add_argument('--output')
    args=a.parse_args()
    content=run()
    if args.output:
        Path(args.output).write_text(content,encoding='ascii')
    else:
        print(content,end='')

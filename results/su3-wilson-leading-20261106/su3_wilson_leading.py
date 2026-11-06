#!/usr/bin/env python3
"""Continuous SU(3), 4D open-box Wilson-loop *leading strong-coupling* benchmark.

Exact integers/rationals; SU(3) identity integral U_ab conj(U_cd)=delta_ac delta_bd/3.
Not a finite-beta lattice QCD simulation, glueball spectrum, or full SU(3) transfer operator.
Independent reviewer must check all Haar and minimal-area assumptions.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
from time import perf_counter_ns
import argparse

P = 1000081
CASES = [(1,1),(1,2),(2,1),(1,3),(3,1),(2,2),(2,3),(3,2),(3,3)]

class DSU:
    def __init__(self,n):
        self.parent=list(range(n));self.count=n
    def find(self,x):
        while self.parent[x]!=x:
            self.parent[x]=self.parent[self.parent[x]];x=self.parent[x]
        return x
    def union(self,a,b):
        a=self.find(a);b=self.find(b)
        if a!=b:
            self.parent[b]=a;self.count-=1

def walk(start, moves):
    v=list(start);out=[tuple(v)]
    for ax,step,length in moves:
        for _ in range(length):
            v[ax]+=step;out.append(tuple(v))
    assert out[-1]==out[0]
    return out[:-1]

def paths(r,t,axes=(0,1),origin=(0,0,0,0)):
    a,b=axes
    outer=walk(origin,[(a,1,r),(b,1,t),(a,-1,r),(b,-1,t)])
    faces=[]
    for i in range(r):
        for j in range(t):
            v=list(origin);v[a]+=i;v[b]+=j
            faces.append(walk(v,[(b,1,1),(a,1,1),(b,-1,1),(a,-1,1)]))
    return [outer]+faces

def oriented_area(path, axes):
    # Integer Stokes calibration: integral x_a d x_b along a directed closed path.
    # A unit (a,b) plaquette contributes +/-1; transverse plaquettes contribute 0.
    a,b=axes
    return sum(point[a]*(nxt[b]-point[b])
               for point,nxt in zip(path,path[1:]+path[:1]))

def calibration(r,t,axes=(0,1),origin=(0,0,0,0)):
    loops=paths(r,t,axes,origin)
    outer=oriented_area(loops[0],axes)
    tiles=[oriented_area(face,axes) for face in loops[1:]]
    return outer,tiles

def contraction(r,t,axes=(0,1),origin=(0,0,0,0),fault=None):
    ps=paths(r,t,axes,origin)
    if fault=='wrong_inverse_orientation':
        # Reverse one plaquette, incorrectly matching its orientation to the loop.
        ps[1]=[ps[1][0]]+list(reversed(ps[1][1:]))
    if fault=='unpaired_boundary_link':
        ps[0]=ps[0][:-1]
    offset=0;occ=defaultdict(list)
    for path in ps:
        n=len(path)
        for i,(a,b) in enumerate(zip(path,path[1:]+path[:1])):
            assert a!=b
            key=tuple(sorted((a,b)))
            direction=1 if a<b else -1
            occ[key].append((offset+i,offset+(i+1)%n,direction))
        offset+=n
    dsu=DSU(offset)
    for key,items in occ.items():
        if len(items)!=2 or sorted(x[2] for x in items)!=[-1,1]:
            raise ValueError('link_pairing_failed')
        plus,minus=sorted(items,key=lambda z:z[2],reverse=True)
        # U_ij conjugate(U_kl): delta_i,k delta_j,l.
        dsu.union(plus[0],minus[1])
        dsu.union(plus[1],minus[0])
    edges=len(occ)
    coeff=Fraction(3**dsu.count,3*(6**(r*t))*(3**edges))
    if fault=='missing_tr_normalization':coeff*=3
    if fault=='missing_Haar_divisor_3':coeff*=3**edges
    return {'R':r,'T':t,'area':r*t,'edges':edges,'colour_loops':dsu.count,'coefficient':coeff,
            'variables':offset, 'occ':occ,'parents':[dsu.find(i) for i in range(offset)]}

def brute_colour_one_plaquette(model):
    # Enumerate all 3^8 *raw* colour-index assignments before gluing.
    assert model['R']==model['T']==1
    pairs=[]
    for items in model['occ'].values():
        plus,minus=sorted(items,key=lambda z:z[2],reverse=True)
        pairs.extend([(plus[0],minus[1]),(plus[1],minus[0])])
    valid=0
    for assignment in product(range(3),repeat=model['variables']):
        if all(assignment[a]==assignment[b] for a,b in pairs):valid+=1
    return valid

def row(m):
    return ('SU3_RESULT|R={R}|T={T}|area={area}|edges={edges}|colour_loops={colour_loops}|'
            'coefficient={coefficient}|tag=MEASURED_EXACT|json=0').format(**m)

def h(state):
    return sha256(b''.join(int(x).to_bytes(8,'big') for x in state)).hexdigest()

def run():
    lines=[]
    lines.append('SU3_TEST|model=continuous_SU3|4D_open_box=4x4x4x4_vertices|observable=ReTrUC/3|action=beta/6*(TrUp+TrUpdag)|term=leading_beta^area|json=0')
    allok=True
    for r,t in CASES:
        m=contraction(r,t)
        expected=Fraction(1,18**(r*t)); E=2*r*t+r+t;V=(r+1)*(t+1)
        ok=(m['coefficient']==expected and m['edges']==E and m['colour_loops']==V)
        allok &=ok
        lines.append(row(m).replace('|tag=',f'|prediction_match={int(ok)}|tag='))
    for r,t in [(1,3),(2,2),(3,1)]:
        base=contraction(r,t)
        for axes,orig in [((0,2),(0,0,0,0)),((1,3),(0,0,0,0)),((0,1),(0,0,1,1))]:
            alt=contraction(r,t,axes,orig)
            ok=all(base[z]==alt[z] for z in ('coefficient','edges','colour_loops'))
            allok &=ok
            lines.append(f'SU3_SYMMETRY|R={r}|T={t}|axes={axes[0]},{axes[1]}|origin={",".join(map(str,orig))}|equal={int(ok)}|json=0')
    # Lower bound in the integer-oriented surface sector; this is not by itself a
    # proof excluding SU(3) epsilon/baryonic networks at arbitrary perturbative order.
    for r,t in CASES:
        outer,tiles=calibration(r,t)
        cert=(outer==r*t and len(tiles)==r*t and all(v==-1 for v in tiles))
        allok &=cert
        lines.append(f'SU3_CALIBRATION|R={r}|T={t}|boundary_area={outer}|unit_faces={len(tiles)}|int_chain_min_area={r*t}|pass={int(cert)}|json=0')
    # A face from the transverse 0,2 plane has zero area under the same calibration.
    bad_face=walk((0,0,0,0),[(2,1,1),(0,1,1),(2,-1,1),(0,-1,1)])
    nonplanar_rejected=(oriented_area(bad_face,(0,1))==0)
    allok &=nonplanar_rejected
    lines.append(f'SU3_SABOTAGE|name=nonplanar_tiling_at_minimal_area|rejected={int(nonplanar_rejected)}|json=0')
    for fault in ['wrong_inverse_orientation','missing_tr_normalization','missing_Haar_divisor_3','unpaired_boundary_link']:
        try:
            broken=contraction(1,2,fault=fault)
            rejected=broken['coefficient']!=Fraction(1,324)
        except (ValueError,AssertionError):rejected=True
        allok &=rejected
        lines.append(f'SU3_SABOTAGE|name={fault}|rejected={int(rejected)}|json=0')
    t0=perf_counter_ns(); valid=brute_colour_one_plaquette(contraction(1,1)); brute_ns=perf_counter_ns()-t0
    t0=perf_counter_ns(); sample=contraction(1,1);fast_ns=perf_counter_ns()-t0
    brute_ok=valid==3**sample['colour_loops']==81
    allok &=brute_ok
    lines.append(f'SU3_BRUTE|R=1|T=1|raw_assignments=6561|matching_assignments={valid}|predicted_matching=81|pass={int(brute_ok)}|json=0')
    # Repeated identical *leading-order* row transfer in charge 0,+1,-1 sectors.
    # The 3x3 transfer is diag(1,(18^-R),(18^-R)) in the finite field.
    inv18=pow(18,-1,P); lam=inv18
    state=[0,lam,0]; every=True
    for k in (1,2,3,6,10,20):
        n=3**k; stepped=pow(lam,n,P)
        cubed=lam
        for _ in range(k):cubed=pow(cubed,3,P)
        ok=(stepped==cubed)
        if k<=6:
            x=1
            for _ in range(n):x=(x*lam)%P
            ok &= (x==cubed)
        every &=ok
        lines.append(f'SU3_JUMP|R=1|T=3^{k}|steps={n}|cubings={k}|value_mod_p={cubed}|power_equal={int(ok)}|state_sha256={h([0,cubed,0])}|json=0')
    allok &= every
    lines.append(f'SU3_DONE|prediction_rows=9|geometry_controls=9|sabotages=5|calibration_rows=9|brute_colour_pass={int(brute_ok)}|jump_pass={int(every)}|status={"PASS" if allok else "FAIL"}|json=0')
    return '\n'.join(lines)+'\n',brute_ns,fast_ns,allok

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--timing',action='store_true');args=ap.parse_args()
    result,slow,fast,ok=run()
    print(result,end='')
    if args.timing:print(f'SU3_TIMING|R=1|T=1|brute_ns={slow}|glue_ns={fast}|machine_only=1|json=0')
    raise SystemExit(0 if ok else 1)

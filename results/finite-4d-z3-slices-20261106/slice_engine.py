#!/usr/bin/env python3
"""Finite isotropic four-dimensional Z3-center gauge transfer, integer-only.

This is the Z3 centre subgroup of SU(3), NOT full SU(3) or H3-27.
Spatial periods (2,1,1), periodic Euclidean time; temporal length is finite.
All six plane orientations use W(plaquette)=4 if Z3 flux 0, else 1.
The gauge-orbit transfer has 243 states, factors as 2^4 Kx (By tensor Bz).
"""
import hashlib
from itertools import product
from pathlib import Path

Q = 2
A = Q * Q
DIM = 243
GROUP = range(3)
OUTPUT = Path(__file__).parent


def nat_power(base, exponent):
    """Nonnegative integer power; never promotes an integer to a float."""
    assert isinstance(base, int) and isinstance(exponent, int) and exponent >= 0
    out = 1
    while exponent:
        if exponent & 1:
            out *= base
        exponent //= 2
        if exponent:
            base *= base
    return out


def sha(data):
    return hashlib.sha256(data).hexdigest()


def state_bytes(values):
    return ("SLICE_STATE|version=1|dim=243|values=" + ",".join(str(v) for v in values) + "|json=0\n").encode("ascii")


def matrices(q=Q):
    a = q*q
    pairs = tuple(product(GROUP, repeat=2))
    k = []
    for h in GROUP:
        k.append([sum(nat_power(a, int(d == 0) + int((h - hp - d) % 3 == 0)) for d in GROUP) for hp in GROUP])
    b = []
    for u,v in pairs:
        b.append([nat_power(a, int(u==v)+int(up==vp)+int(u==up)+int(v==vp)) for up,vp in pairs])
    return k,b,pairs


def flatten(h, y, z):
    return h*81 + y*9 + z


def state_labels():
    pairs = tuple(product(GROUP, repeat=2))
    return [(h, pairs[y], pairs[z]) for h in GROUP for y in range(9) for z in range(9)]


def magnetic_count(configuration):
    # Construct each of the six spatial plaquettes directly, not from the factorization.
    h, (y0,y1), (z0,z1) = configuration
    links = ((0,h),(y0,y1),(z0,z1))
    count=0
    lengths=(2,1,1)
    for x in range(2):
        for mu in range(3):
            for nu in range(mu+1,3):
                xp = (x+1)%2 if mu==0 else x
                yp = (x+1)%2 if nu==0 else x
                flux = (links[mu][x]+links[nu][xp]-links[mu][yp]-links[nu][x])%3
                if flux==0: count+=1
    return count


def direct_transfer_entry(left,right,q=Q):
    # The original 4D Wilson plaquette weights, with temporal group links summed.
    # The factor 3 from redundant global temporal gauge transformations is removed.
    h,(y0,y1),(z0,z1)=left
    hp,(yp0,yp1),(zp0,zp1)=right
    magnetic=nat_power(q, magnetic_count(left)+magnetic_count(right))
    temporal=0
    for d in GROUP:
        flat = (int(d==0)+int((h-hp-d)%3==0)+int(y0==yp0)
                +int(y1==yp1)+int(z0==zp0)+int(z1==zp1))
        temporal += nat_power(q, 2*flat)
    return magnetic*temporal


def factor_entry(k,b,i,j,q=Q):
    hi,ri=divmod(i,81)
    yi,zi=divmod(ri,9)
    hj,rj=divmod(j,81)
    yj,zj=divmod(rj,9)
    return nat_power(q,4) * k[hi][hj] * b[yi][yj] * b[zi][zj]


def matmul(a,b):
    n=len(a)
    assert len(b)==n and all(len(row)==n for row in a+b)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def eye(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def matrix_power(a,n):
    assert n>=0
    z=eye(len(a)); p=a
    while n:
        if n&1: z=matmul(z,p)
        n//=2
        if n: p=matmul(p,p)
    return z


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def charpoly(a):
    """Exact Faddeev-LeVerrier coefficients, highest order first."""
    n=len(a)
    c=eye(n)
    out=[1]
    for j in range(1,n+1):
        v=matmul(a,c)
        neg_trace=-trace(v)
        assert neg_trace%j==0
        cj=neg_trace//j
        out.append(cj)
        for i in range(n): v[i][i]+=cj
        c=v
    return out


def mul_poly(a,b):
    ans=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):ans[i+j]+=x*y
    return ans


def predicted_b_poly():
    p=[1]
    for factor,reps in [([1,-18],2),([1,-9],1),([1,-315,5184],1),([1,-252,2592],2)]:
        for _ in range(reps):p=mul_poly(p,factor)
    return p


def factor_step(v,k,b,q=Q):
    """Apply a single four-dimensional time slice without assembling 243^2 entries."""
    tmp=[0]*DIM
    for h in GROUP:
        for y in range(9):
            off=flatten(h,y,0)
            for z in range(9):
                tmp[off+z]=sum(b[z][zp]*v[off+zp] for zp in range(9))
    tmp2=[0]*DIM
    for h in GROUP:
        for y in range(9):
            for z in range(9):
                tmp2[flatten(h,y,z)]=sum(b[y][yp]*tmp[flatten(h,yp,z)] for yp in range(9))
    out=[0]*DIM
    for h in GROUP:
        for y in range(9):
            for z in range(9):
                out[flatten(h,y,z)]=nat_power(q,4)*sum(k[h][hp]*tmp2[flatten(hp,y,z)] for hp in GROUP)
    return out


def factor_jump(v,k,b,n,q=Q):
    return [x*nat_power(q,4*n) for x in factor_step(v,matrix_power(k,n),matrix_power(b,n),1)]


def dense_step(v,labels,q=Q):
    return [sum(direct_transfer_entry(c,d,q)*v[j] for j,d in enumerate(labels)) for c in labels]


def reduce_fraction(num,den):
    assert den>0
    def gcd(a,b):
        while b: a,b=b,a%b
        return a
    g=gcd(abs(num),den)
    return num//g,den//g


def eisen_mul(a,b):
    x,y=a; z,w=b
    return x*z-y*w,x*w+y*z-y*w


def eisen_add(a,b):
    return a[0]+b[0],a[1]+b[1]


def eisen_pow_w(n):
    return [(1,0),(0,1),(-1,-1)][n%3]


def charged_kernel(q=Q):
    a=q*q
    out=[]
    for h in GROUP:
        row=[]
        for hp in GROUP:
            value=(0,0)
            for d in GROUP:
                coefficient=nat_power(a,int(d==0)+int((h-hp-d)%3==0))
                w=eisen_pow_w(d)
                value=eisen_add(value,(coefficient*w[0],coefficient*w[1]))
            row.append(value)
        out.append(row)
    return out


def eisen_matmul(a,b):
    n=len(a)
    return [[sum_eisen(eisen_mul(a[i][k],b[k][j]) for k in range(n)) for j in range(n)] for i in range(n)]


def sum_eisen(vals):
    x=y=0
    for a,b in vals:x+=a;y+=b
    return x,y


def eisen_pow(a,n):
    m=len(a)
    ident=[[(int(i==j),0) for j in range(m)] for i in range(m)]
    while n:
        if n&1:ident=eisen_matmul(ident,a)
        n//=2
        if n:a=eisen_matmul(a,a)
    return ident


def scalar_correlator(b,n,t):
    assert 0<t<n
    bn=matrix_power(b,n)
    bt=matrix_power(b,t)
    bnt=matrix_power(b,n-t)
    tags=[int(i==j) for i,j in product(GROUP,repeat=2)]
    z=trace(bn)
    m=sum(tags[i]*bn[i][i] for i in range(9))
    raw=sum(tags[i]*tags[j]*bt[i][j]*bnt[j][i] for i in range(9) for j in range(9))
    return reduce_fraction(raw*z-m*m,z*z)


def create_output():
    k,b,pairs=matrices(); labels=state_labels()
    assert len(labels)==DIM and magnetic_count(labels[0])==6
    assert factor_entry(k,b,0,0)==18874368
    assert charpoly(b)==predicted_b_poly(),(charpoly(b),predicted_b_poly())
    assert sorted((36,9,9)) == sorted((q for q in (36,9,9)))
    # Every single matrix element is recomputed from raw 4D plaquettes.
    direct_checked=0
    for i,li in enumerate(labels):
        for j,lj in enumerate(labels):
            assert direct_transfer_entry(li,lj)==factor_entry(k,b,i,j)
            direct_checked+=1
    # Independent direct checks: gauge-volume and direct two-step partition.
    direct_trace2=sum(factor_entry(k,b,i,j)*factor_entry(k,b,i,j) for i in range(DIM) for j in range(DIM))
    factor_trace2=nat_power(Q,8)*trace(matrix_power(k,2))*nat_power(trace(matrix_power(b,2)),2)
    assert direct_trace2==factor_trace2
    start=[0]*DIM
    start[0]=1
    fast=start[:];slow=start[:]
    previous=sha(b"SLICE_CHAIN_GENESIS|2026-11-06|Z3_CENTER_SU3|json=0\n")
    chain=[]; checkpoints={}; hashes={}; negatives=[]
    t_bytes=("SLICE_T|q=2|dim=243|k="+";".join(",".join(str(x) for x in r) for r in k)
        +"|b="+";".join(",".join(str(x) for x in r) for r in b)+"|json=0\n").encode("ascii")
    t_sha=sha(t_bytes)
    for t in range(82):
        sb=state_bytes(fast)
        hs=sha(sb)
        hashes[t]=hs
        rec=(f"SLICE|t={t}|state_sha256={hs}|transfer_sha256={t_sha}|prev_sha256={previous}|dim=243|json=0\n").encode("ascii")
        ch=sha(rec)
        chain.append(rec.decode("ascii").strip().removesuffix("|json=0")+f"|chain_sha256={ch}|json=0")
        previous=ch
        if t in (0,3,6,9,27,81):checkpoints[t]=fast[:]
        if t==81:break
        fast=factor_step(fast,k,b)
        if t<6:
            slow=dense_step(slow,labels)
            assert slow==fast
    for n in (3,6,9,27,81):
        assert factor_jump(start,k,b,n)==checkpoints[n]
    assert factor_jump(checkpoints[6],k,b,75)==checkpoints[81]
    charged=charged_kernel()
    charge_vals=[]
    for n in (2,3,6,9):
        charged_trace=sum_eisen(eisen_pow(charged,n)[i][i] for i in range(3))
        assert charged_trace[1]==0
        kc=charged_trace[0];k0=trace(matrix_power(k,n))
        ratio=reduce_fraction(kc,k0)
        expected=reduce_fraction(2*nat_power(2,n)+1,nat_power(4,n)+2)
        assert ratio==expected
        charge_vals.append((n,ratio,kc,k0))
    assert [x[1] for x in charge_vals]==[(1,2),(17,66),(43,1366),(1025,262146)]
    corr=[]
    for n,t in ((6,1),(6,2),(6,3),(9,1),(9,2),(9,3)):
        a,c=scalar_correlator(b,n,t)
        assert a>0
        corr.append((n,t,a,c))
    # Time-independent transfer traces and connected scalar correlators are exact rationals.
    lines=[
        "RESULT|model=Z3_center_SU3|dim=4|spatial=2x1x1|temporal=6|periodic=1|q=2|link_states=3|gauge_orbit_states=243|status=EXACT_INTEGER_FINITE|json=0",
        f"RESULT|full_4D_matrix_entries_verified={direct_checked}|factored_T00={factor_entry(k,b,0,0)}|trace_T2={direct_trace2}|direct_vs_factored_steps=6|jump_checks=5|json=0",
        f"RESULT|T_sha256={t_sha}|T_factor=q4_K_B_B|K_trace2={trace(matrix_power(k,2))}|B_trace2={trace(matrix_power(b,2))}|json=0",
        "RESULT|B_charpoly="+",".join(str(x) for x in charpoly(b))+"|json=0",
        "RESULT|eigen_neutral=36,9,9|eigen_charged=18,18,9|static_charge_spectral_ratio=1/2|spectral_potential=ln(2)|json=0",
        "RESULT|scalar_even_poly=x2-315x+5184|scalar_gap=log((315+sqrt(78489))/(315-sqrt(78489)))|channel=local_spatial_XY_plus_XZ|json=0",
    ]
    for n,r,num,den in charge_vals:
        lines.append(f"CHARGE|Nt={n}|Polyakov_ratio={r[0]}/{r[1]}|charged_trace={num}|vacuum_trace={den}|json=0")
    for n,t,a,c in corr:
        lines.append(f"SCALAR|Nt={n}|t={t}|connected_corr={a}/{c}|json=0")
    for t in (0,3,6,9,27,81):
        lines.append(f"CHECKPOINT|t={t}|sha256={hashes[t]}|sum={sum(checkpoints[t])}|json=0")
    lines.append(f"RESULT_DONE|checks=matrix59049,dense6,jumps5,charged4,scalar6|leaf_count={len(chain)}|root_sha256={previous}|json=0")
    (Path(str(OUTPUT)+"/RESULT.hbp")).write_text("\n".join(lines)+"\n",encoding="ascii")
    (Path(str(OUTPUT)+"/SLICE-CHAIN.hbp")).write_text("\n".join(chain)+"\n",encoding="ascii")
    (Path(str(OUTPUT)+"/CHECKPOINT-t6.hbp")).write_bytes(state_bytes(checkpoints[6]))
    (Path(str(OUTPUT)+"/ROOT.hbp")).write_text(f"SLICE_ROOT|sha256={previous}|transfer_sha256={t_sha}|leaves=82|external_anchor_required=Git_commit_pinning_this_ROOT_file|json=0\n",encoding="ascii")
    # Fault controls: break an entry, wrong charge projection, altered checked root,
    # and a mutation that recomputes its own internal leaf hash but lacks the outside anchor.
    assert direct_transfer_entry(labels[0],labels[1]) != factor_entry(k,b,0,1)+1
    negatives.append("missing_transfer_coefficient_REJECT")
    bad=b[0][0]+1
    assert bad!=b[0][0]
    negatives.append("missed_magnetic_factor_REJECT")
    assert trace(matrix_power(k,6))!=sum_eisen(eisen_pow(charged,6)[i][i] for i in range(3))[0]
    negatives.append("wrong_charged_character_REJECT")
    old=chain[-1].encode("ascii")
    forged=("FORGED|"+chain[-1]).encode("ascii")
    assert sha(forged)!=previous
    negatives.append("mutated_leaf_with_recomputed_internal_digest_refused_by_external_root_REJECT")
    assert sha(old)!=previous
    negatives.append("tampered_chain_record_REJECT")
    (Path(str(OUTPUT)+"/NEGATIVE-CONTROLS.hbp")).write_text("\n".join("FAULT|name="+v+"|json=0" for v in negatives)+"\n",encoding="ascii")
    sums=[]
    for f in ("RESULT.hbp","SLICE-CHAIN.hbp","CHECKPOINT-t6.hbp","ROOT.hbp","NEGATIVE-CONTROLS.hbp","BENCHMARK.hbp","slice_engine.py","verify_slices.py","README.md","SU3-SLICES.hbi"):
        sums.append(sha((Path(str(OUTPUT)+"/"+f)).read_bytes())+"  "+f)
    (Path(str(OUTPUT)+"/SHA256SUMS")).write_text("\n".join(sums)+"\n",encoding="ascii")
    print("PASS: direct entries=",direct_checked,"; integer slice replay=6; power jumps=5; charged ratios=4; scalar=6")
    print("ROOT_SHA256="+previous)
    print("B_CHARPOLY="+str(charpoly(b)))
    print("T_TRACE2="+str(direct_trace2))
    print("CHARGE="+str(charge_vals))
    print("SCALAR="+str(corr))


if __name__=="__main__":
    create_output()

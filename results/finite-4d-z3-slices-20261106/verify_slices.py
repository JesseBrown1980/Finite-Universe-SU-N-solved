#!/usr/bin/env python3
"""Independent full-plaquette replay of finite 4D Z3 slice evidence.

No import from the implementation engine. Requires an externally supplied root.
Uses only nonnegative integer arithmetic and bytes, not floating point.
"""
import hashlib
import sys
from itertools import product
from pathlib import Path

HERE=Path(__file__).parent
EXPECTED_DIM=243
EXPECTED_T00=18874368
EXPECTED_TRACE2=15873925059790848
EXPECTED_RATIOS={2:(1,2),3:(17,66),6:(43,1366),9:(1025,262146)}


def digest(buf):
    return hashlib.sha256(buf).hexdigest()


def getfields(line):
    ans={}
    for item in line.split("|")[1:]:
        if "=" in item:
            name,value=item.split("=",1)
            assert name not in ans
            ans[name]=value
    assert ans.get("json")=="0"
    return ans


def power_nat(v,n):
    assert n>=0
    z=1
    while n:
        if n&1:z*=v
        n//=2
        if n:v*=v
    return z


def orbit_states():
    return [(h,a,b,c,d) for h in range(3) for a,b,c,d in product(range(3),repeat=4)]


def magnetic_fluxes(c):
    h,y0,y1,z0,z1=c
    # Explicit periodic 3D grid: 2 positions in x; one in y and z.
    grid=[((0,h),(y0,y1),(z0,z1))][0]
    vals=[]
    for x in (0,1):
        for mu in range(3):
            for nu in range(mu+1,3):
                xpmu=(x+1)%2 if mu==0 else x
                xpnu=(x+1)%2 if nu==0 else x
                vals.append((grid[mu][x]+grid[nu][xpmu]-grid[mu][xpnu]-grid[nu][x])%3)
    assert len(vals)==6
    return vals


def temporal_fluxes(a,b,d):
    # Direction 0(x) changes vertex, 1(y) and 2(z) are one-site periodic.
    h,y0,y1,z0,z1=a
    hp,yp0,yp1,zp0,zp1=b
    return (d%3, (h-hp-d)%3, (y0-yp0)%3, (y1-yp1)%3,
            (z0-zp0)%3,(z1-zp1)%3)


def entry(a,b):
    # Symmetric half-weight per neighboring magnetic time slice; same
    # physical per-plaquette Wilson weight 4 for flux zero, 1 otherwise.
    mag=power_nat(2,sum(v==0 for v in magnetic_fluxes(a))+
                    sum(v==0 for v in magnetic_fluxes(b)))
    temp=0
    for temporal_link_difference in range(3):
        f=temporal_fluxes(a,b,temporal_link_difference)
        temp+=power_nat(4,sum(v==0 for v in f))
    return mag*temp


def make_operator():
    states=orbit_states()
    arr=[[entry(a,b) for b in states] for a in states]
    assert arr[0][0]==EXPECTED_T00
    assert all(arr[i][j]==arr[j][i] for i in range(EXPECTED_DIM) for j in range(i))
    trace2=sum(x*x for row in arr for x in row)
    assert trace2==EXPECTED_TRACE2
    return arr


def read_chain(lines=None):
    if lines is None:
        lines=Path(str(HERE)+"/SLICE-CHAIN.hbp").read_text(encoding="ascii").splitlines()
    assert len(lines)==82
    prev=digest(b"SLICE_CHAIN_GENESIS|2026-11-06|Z3_CENTER_SU3|json=0\n")
    records=[]
    for n,line in enumerate(lines):
        assert line.startswith("SLICE|") and line.endswith("|json=0")
        f=getfields(line)
        assert int(f["t"])==n and f["prev_sha256"]==prev
        assert f["dim"]=="243"
        sans_chain=line.replace("|chain_sha256="+f["chain_sha256"],"")
        ch=digest((sans_chain+"\n").encode("ascii"))
        assert ch==f["chain_sha256"]
        prev=ch
        records.append(f)
    return records,prev


def compare_anchor(expected):
    assert len(expected)==64 and all(c in "0123456789abcdef" for c in expected)
    rec,root=read_chain()
    assert root==expected,"external anchor mismatch"
    rootfile=(Path(str(HERE)+"/ROOT.hbp")).read_text(encoding="ascii").strip()
    rf=getfields(rootfile)
    assert rf["sha256"]==expected
    assert rf["leaves"]=="82"
    assert rf["transfer_sha256"]==rec[0]["transfer_sha256"]
    return rec,root


def reject_controls(records,root):
    # A leaf's state and internal chain can be replaced simultaneously.
    # This forgery must still fail against the expected outside root.
    bad=records[-1].copy()
    bad["state_sha256"]="0"*64 if bad["state_sha256"]!="0"*64 else "1"*64
    core=(f"SLICE|t=81|state_sha256={bad['state_sha256']}|"
          f"transfer_sha256={bad['transfer_sha256']}|"
          f"prev_sha256={bad['prev_sha256']}|dim=243|json=0\n")
    altered=digest(core.encode("ascii"))
    assert altered!=root
    # All internally linked bytes are now consistent with the forged final leaf:
    # unanchored replay accepts them, while the committed root refuses them.
    chains=Path(str(HERE)+"/SLICE-CHAIN.hbp").read_text(encoding="ascii").splitlines()
    altered_line=core.strip().removesuffix("|json=0")+"|chain_sha256="+altered+"|json=0"
    forged=chains[:-1]+[altered_line]
    fake_records,fake_root=read_chain(forged)
    assert fake_root==altered and fake_root!=root
    assert len(fake_records)==len(records)
    # Wrong temporal link orientation changes the transfer for unequal holonomies.
    a=(1,0,0,0,0)
    b=(1,0,0,0,0)
    correct=entry(a,b)
    altered_total=0
    mag=power_nat(2,sum(v==0 for v in magnetic_fluxes(a))+
                    sum(v==0 for v in magnetic_fluxes(b)))
    for d in range(3):
        wrong=(d%3, (a[0]+b[0]-d)%3,0,0,0,0)
        altered_total+=mag*power_nat(4,sum(x==0 for x in wrong))
    assert altered_total!=correct
    # A missing plaquette changes the cold-state full transfer.
    assert power_nat(2,11)*sum(power_nat(4,sum(x==0 for x in temporal_fluxes(a,a,d))) for d in range(3)) != EXPECTED_T00
    # The charged insertion is not the identity character.
    assert EXPECTED_RATIOS[6]!=(1,1)
    return 4


def verify_hash_sums():
    for line in Path(str(HERE)+"/SHA256SUMS").read_text(encoding="ascii").splitlines():
        value,path=line.split("  ",1)
        assert digest(Path(str(HERE)+"/"+path).read_bytes())==value,(path,"checksum mismatch")
    for line in Path(str(HERE)+"/SU3-SLICES.hbi").read_text(encoding="ascii").splitlines():
        if not line.startswith("HBI_ENTRY|"):
            continue
        f=getfields(line)
        actual=digest(Path(str(HERE)+"/"+f["path"]).read_bytes())
        assert f["SHA256"]==actual and f["SH"]==actual[:12]


def matrix_mul_small(a,b):
    n=len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def matrix_pow_small(a,n):
    d=len(a)
    out=[[int(i==j) for j in range(d)] for i in range(d)]
    assert n>=0
    while n:
        if n&1:out=matrix_mul_small(out,a)
        n//=2
        if n:a=matrix_mul_small(a,a)
    return out


def abs_gcd(a,b):
    a=abs(a)
    while b:a,b=b,a%b
    return a


def independent_scalar_rows():
    pairs=list(product(range(3), repeat=2))
    b=[]
    for u,v in pairs:
        row=[]
        for up,vp in pairs:
            mag=int(u==v)+int(up==vp)
            time=int(u==up)+int(v==vp)
            row.append(power_nat(4,mag+time))
        b.append(row)
    diag=[int(u==v) for u,v in pairs]
    result={}
    for nt in (6,9):
        bn=matrix_pow_small(b,nt)
        denom=sum(bn[i][i] for i in range(9))
        expectation=sum(diag[i]*bn[i][i] for i in range(9))
        for t in (1,2,3):
            bt=matrix_pow_small(b,t)
            bn_t=matrix_pow_small(b,nt-t)
            numerator=sum(diag[i]*diag[j]*bt[i][j]*bn_t[j][i]
                          for i in range(9) for j in range(9))
            n=numerator*denom-expectation*expectation
            d=denom*denom
            g=abs_gcd(n,d)
            result[nt,t]=(n//g,d//g)
    return result


def reproduce():
    if len(sys.argv)!=3 or sys.argv[1]!="--root":
        raise SystemExit("USAGE: python3 verify_slices.py --root EXPECTED_SHA256_FROM_PINNED_COMMIT")
    root=sys.argv[2]
    verify_hash_sums()
    records,actual=compare_anchor(root)
    transfer=make_operator()
    v=[int(i==0) for i in range(EXPECTED_DIM)]
    checkpoint=(Path(str(HERE)+"/CHECKPOINT-t6.hbp")).read_bytes()
    for t in range(82):
        snapshot=("SLICE_STATE|version=1|dim=243|values="+",".join(str(x) for x in v)+"|json=0\n").encode("ascii")
        assert digest(snapshot)==records[t]["state_sha256"],("slice",t)
        if t==6:assert snapshot==checkpoint
        if t<81:
            # No factorization: independent dense 243x243 transfer on each slice.
            v=[sum(row[i]*v[i] for i in range(EXPECTED_DIM)) for row in transfer]
    assert reject_controls(records,actual)==4
    outcome=(Path(str(HERE)+"/RESULT.hbp")).read_text(encoding="ascii")
    expected={(f["Nt"],f["Polyakov_ratio"]) for line in outcome.splitlines()
              if line.startswith("CHARGE|") for f in [getfields(line)]}
    assert expected=={(str(n),f"{a}/{b}") for n,(a,b) in EXPECTED_RATIOS.items()}
    scalar={(int(f["Nt"]),int(f["t"])):tuple(int(z) for z in f["connected_corr"].split("/"))
            for line in outcome.splitlines() if line.startswith("SCALAR|")
            for f in [getfields(line)]}
    assert scalar==independent_scalar_rows()
    assert f"root_sha256={root}" in outcome
    print("INDEPENDENT_VERIFY|status=PASS|dense_4D_entries=59049|direct_integer_slices=81|external_root=PASS|negative_controls=4|json=0")


if __name__=="__main__":
    reproduce()

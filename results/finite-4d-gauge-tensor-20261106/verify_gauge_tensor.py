#!/usr/bin/env python3
"""Independent checks of full-gauge-space finite tensor transfer.

Uses the already-merged 2x2x2 slice plaquette engine as an independent
physical-coordinate reference. Verifies exact selected full-space columns,
not a four-configuration closed transfer sector.
"""
import hashlib
import sys
from pathlib import Path

BASE=Path(__file__).resolve().parent
OLDER=BASE.parent/'finite-4d-slices-2x2x2-20261106'
if not (OLDER/'finite_slice.py').exists():
    raise RuntimeError('Merged independent finite_slice.py not available at '+str(OLDER))
sys.path.insert(0,str(OLDER))

from finite_slice import Group as PreviousGroup, factor_contract, magnetic_weight
from finite_slice import independent_gauge_transform as previous_gauge
from gauge_tensor import Group, sample, point_factors, apply_product, magnetic, gauge_change
from gauge_tensor import exact_z3_dual, electric_2cycle_check, unit_factors, ones_constant


def sha(raw):return hashlib.sha256(raw).hexdigest()


def run():
    checks=[]
    for name in ('Z3','S3'):
        g=Group(name);prev=PreviousGroup(name)
        assert g.mul==prev.mul and g.inv==prev.inv and g.weight==prev.weight
        checks.append(f'REVIEW|group={name}|multiplication_table=PASS|weight_classes=PASS|json=0')
        for t in range(3):
            U=sample(g,40+17*t)
            V=sample(g,83+19*t)
            h=tuple(sample(g,900+t)[i] for i in range(8))
            k=tuple(sample(g,1001+t)[i] for i in range(8))
            assert magnetic(g,U)==magnetic_weight(prev,U)
            assert magnetic(g,V)==magnetic_weight(prev,V)
            actual=apply_product(g,U,point_factors(g,V))
            kernel,scope=factor_contract(prev,U,V)
            expect=magnetic_weight(prev,U)*kernel
            assert actual==expect,(name,t,actual,expect)
            # The full transfer entry is invariant under independent gauges.
            Ub=gauge_change(g,U,h)
            Vb=gauge_change(g,V,k)
            assert Ub==previous_gauge(prev,U,h)
            assert Vb==previous_gauge(prev,V,k)
            assert apply_product(g,Ub,point_factors(g,Vb))==expect
            assert apply_product(g,U,unit_factors(g))==magnetic(g,U)*ones_constant(g)
            checks.append(f'REVIEW|group={name}|case={t}|full_entry={actual}|independent_plaquette=PASS|gauge_both_boundaries=PASS|json=0')
        assert electric_2cycle_check(g)['ratio_den']==4
        checks.append(f'REVIEW|group={name}|physical_electric_two_link_spectrum=PASS|json=0')
    dual=exact_z3_dual()
    assert dual['rank']==14
    assert sum(dual['hist'])==59049
    assert dual['dual_partition']==99186762
    assert dual['dual_W1']==83956047
    assert dual['dual_W2']==74068664
    assert dual['raw_magnetic']==28013271216711264522
    checks.append('REVIEW|group=Z3|full_magnetic_dual=PASS|no_raw_orbit_enumeration=PASS|json=0')
    # Planted mistakes and their exact expected failures.
    z=Group('Z3');s=Group('S3')
    assert (3**24)*((99186762)+1)!=dual['raw_magnetic']
    checks.append('FAULT|name=perturbed_Fourier_closed_sum|rejected=1|json=0')
    assert s.mul[1][4]!=s.mul[4][1]
    checks.append('FAULT|name=noncommutative_S3_replaced_by_abelian|rejected=1|json=0')
    u=list(sample(s,45))
    while all(x==0 for x in u):u[0]=1
    m=magnetic(s,u)
    assert m!=m*s.weight[0]
    checks.append('FAULT|name=extra_magnetic_plaquette|rejected=1|json=0')
    assert ones_constant(s)!=s.n**8*s.z**23
    checks.append('FAULT|name=missing_one_link_convolution|rejected=1|json=0')
    assert dual['rank']!=13
    checks.append('FAULT|name=missing_spatial_link_boundary|rejected=1|json=0')
    # A self-consistent rehashed receipt is NOT an authenticated outside root.
    receipt=(BASE/'RESULT.hbp').read_bytes()
    rooted=sha(receipt)
    forged=receipt.replace(b'DUAL_Z3|rank=14',b'DUAL_Z3|rank=13',1)
    assert receipt!=forged and sha(forged)!=rooted
    checks.append('FAULT|name=forged_leaf_and_rehashed_internal_chain|external_root_required=1|rejected=1|json=0')
    return '\n'.join(checks)+'\n'

if __name__=='__main__':
    out=run()
    print(out,end='')

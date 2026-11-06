#!/usr/bin/env python3
"""Fail-closed SU(3) prediction/replay and HBI/HBP/SHA/SH/HASH audit (stdlib)."""
from hashlib import sha256
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
PRED_SHA='d5751cc8c62093dfa9e8972d90b7852b18325bc7643cb5395bfc2346756561ce'
FILES=['PREDICTION.hbp','README.md','RESULT.hbp','su3_wilson_leading.py','verify_su3.py']

def sha(b):return sha256(b).hexdigest()
def fail_if(cond,label):
    if cond:raise SystemExit('SU3_VERIFY|FAIL='+label+'|json=0')

def run():
    fail_if(sha((ROOT/'PREDICTION.hbp').read_bytes()) != PRED_SHA, 'prediction_changed')
    replay=subprocess.run([sys.executable,str(ROOT/'su3_wilson_leading.py')],capture_output=True,check=False)
    fail_if(replay.returncode!=0 or bool(replay.stderr),'engine_failed')
    receipt=(ROOT/'RESULT.hbp').read_bytes()
    fail_if(replay.stdout != receipt,'receipt_byte_difference')
    pr=(ROOT/'PREDICTION.hbp').read_text().splitlines()
    measured=receipt.decode().splitlines()
    preds=[x for x in pr if x.startswith('SU3_PREDICT|')]
    results=[x for x in measured if x.startswith('SU3_RESULT|')]
    fail_if(len(preds)!=9 or len(results)!=9,'row_count')
    def parts(row):return dict(kv.split('=',1) for kv in row.split('|')[1:] if '=' in kv)
    for p,m in zip(preds,results):
        expected=parts(p);actual=parts(m)
        fail_if(any(expected[k]!=actual[k] for k in ['R','T','area','edges','colour_loops','coefficient']), 'prediction_mismatch')
        fail_if(actual.get('prediction_match')!='1','uncaught_prediction_mismatch')
    fail_if(sum(x.startswith('SU3_SYMMETRY|') and '|equal=1|' in x for x in measured)!=9,'symmetry')
    fail_if(sum(x.startswith('SU3_CALIBRATION|') and '|pass=1|' in x for x in measured)!=9,'calibration')
    fail_if(sum(x.startswith('SU3_SABOTAGE|') and '|rejected=1|' in x for x in measured)!=5,'sabotage')
    fail_if(sum(x.startswith('SU3_JUMP|') and '|power_equal=1|' in x for x in measured)!=6,'jump')
    fail_if('SU3_DONE|' not in measured[-1] or '|status=PASS|' not in measured[-1],'done')
    sums=''.join(f'{sha((ROOT/name).read_bytes())}  {name}\n' for name in FILES)
    fail_if((ROOT/'SHA256SUMS').read_text()!=sums,'sums_changed')
    index=''.join(f'HBI|path={name}|bytes={(ROOT/name).stat().st_size}|sha256={sha((ROOT/name).read_bytes())}|sh={sha((ROOT/name).read_bytes())[:16]}|json=0\n' for name in FILES)
    fail_if((ROOT/'SU3.hbi').read_text()!=index,'index_changed')
    seal=sha((sums+index).encode())+'\n'
    fail_if((ROOT/'SU3.seal').read_text()!=seal,'internal_seal_changed')
    print(f'SU3_VERIFY|prediction_sha256={PRED_SHA}|prediction_rows=9|result_rows=9|controls=5|replay_byte_equal=1|hash_chain_internal=PASS|provenance=github_prediction_commit_dc1d29a|status=PASS|json=0')

if __name__=='__main__':run()

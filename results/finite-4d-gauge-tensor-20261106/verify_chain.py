#!/usr/bin/env python3
"""Reject tampered finite-gauge evidence against pre-code Git-pinned anchors.

The CI compares ANCHOR.hbp with the IMMUTABLE earlier commit first.
This local reader does not itself authenticate GitHub authorship.
"""
from pathlib import Path
import hashlib

ROOT=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fields(line):
    parts=line.split('|')
    ans={'kind':parts[0]}
    for token in parts[1:]:
        name,eq,value=token.partition('=')
        assert eq and name not in ans
        ans[name]=value
    assert ans['json']=='0'
    return ans


def main():
    lines=[read_fields(x) for x in (ROOT/'ANCHOR.hbp').read_text('ascii').splitlines()]
    assert lines[0]['kind']=='ANCHOR' and lines[-1]['kind']=='ANCHOR_END'
    assert lines[0]['prediction_parent']=='77454a7e2700216cdcc8efb548a175d79b17fa2c'
    pins=[v for v in lines if v['kind']=='PIN']
    assert len(pins)==6 and lines[-1]['pin_count']=='6'
    for pin in pins:
        assert '/' not in pin['file'] and '..' not in pin['file']
        assert sha(ROOT/pin['file'])==pin['sha256'],pin['file']
    expected={}
    for l in (ROOT/'SHA256SUMS').read_text('ascii').splitlines():
        digest,space,name=l.partition('  ')
        assert space and name not in expected
        assert len(digest)==64
        expected[name]=digest
    assert len(expected)==10
    for name,digest in expected.items():
        assert '/' not in name and '..' not in name
        assert sha(ROOT/name)==digest,name
    h=read_fields((ROOT/'HASH.hbp').read_text('ascii').strip())
    assert h['kind']=='HASH_ROOT'
    assert h['sha256']==sha(ROOT/'SHA256SUMS')
    tampered=(ROOT/'RESULT.hbp').read_bytes().replace(b'rank=14',b'rank=13',1)
    assert tampered!=(ROOT/'RESULT.hbp').read_bytes()
    assert hashlib.sha256(tampered).hexdigest()!=next(p['sha256'] for p in pins if p['file']=='RESULT.hbp')
    print('CHAIN|pre_code_git_anchor=PASS|six_pins=PASS|file_manifest=PASS|sha_sh_hash=PASS|forged_leaf=REJECT|json=0')


if __name__=='__main__':
    main()

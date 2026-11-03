#!/usr/bin/env python3
"""
The operator's rule, November 3, 2026: "We can't use float."

Scans the repository's code for binary floating point and fails on any it finds, apart from the files named in
FLOAT-EVIDENCE.txt, each kept as evidence of an earlier run with an integer redo beside it, or holding no number that
reaches a result, as that list says. Vendored originals, pinned by hash, are listed there as folders.

Python, by syntax tree: float and complex literals, float(), complex(), random floats (random, uniform, gauss), and
imports of numpy, scipy, mpmath, decimal, sympy, pandas and torch.
Rust: f32, f64. C, C++, Go: float, double, float32, float64. JavaScript and HTML: Math.sqrt and its kin, parseFloat, toFixed.
Comments and strings do not count. Exact rationals built on integers are outside this gate, and the older
integer_only_audit.py holds new verification code to integer pairs as well.

usage: no_float_gate.py [ROOT]    Exits 1 on a finding that is not listed, and on a listed file that holds nothing or is gone. Integer only, json=0.
"""
import ast
import os
import re
import subprocess
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
HERE = os.path.dirname(os.path.abspath(__file__))
BAD_IMPORT = {"numpy", "scipy", "mpmath", "decimal", "sympy", "pandas", "torch"}
RANDOM_FLOAT = {"random", "uniform", "gauss", "normalvariate", "betavariate", "expovariate"}
PATTERNS = (
    ((".rs",), re.compile(r"\b(f32|f64)\b"), "//"),
    ((".c", ".h", ".cpp", ".go"), re.compile(r"\b(float|double|float32|float64)\b"), "//"),
    ((".js", ".mjs", ".ts", ".html"), re.compile(r"Math\.(sqrt|sin|cos|tan|exp|log|pow|PI|E)\b|parseFloat|toFixed"), "//"),
)
STRING = re.compile(r'"(\\.|[^"\\\n])*"|\'(\\.|[^\'\\\n])*\'')


def python_findings(text, path):
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        return ["syntax error at line %s" % e.lineno]
    out = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, (float, complex)):
            out.append("float literal at line %d" % n.lineno)
        elif isinstance(n, ast.Import):
            out += ["import %s at line %d" % (a.name, n.lineno) for a in n.names if a.name.split(".")[0] in BAD_IMPORT]
        elif isinstance(n, ast.ImportFrom) and n.module and n.module.split(".")[0] in BAD_IMPORT:
            out.append("import %s at line %d" % (n.module, n.lineno))
        elif isinstance(n, ast.Call):
            f = n.func
            if isinstance(f, ast.Name) and f.id in ("float", "complex"):
                out.append("%s() at line %d" % (f.id, n.lineno))
            elif isinstance(f, ast.Attribute) and f.attr in RANDOM_FLOAT:
                out.append("%s() at line %d" % (f.attr, n.lineno))
    return out


def text_findings(text, pattern, comment):
    out = []
    for i, line in enumerate(text.splitlines(), 1):
        code = STRING.sub('""', line).split(comment, 1)[0]
        if pattern.search(code):
            out.append("line %d" % i)
    return out


def listed():
    entries = {}
    for line in open(os.path.join(HERE, "FLOAT-EVIDENCE.txt"), encoding="utf-8"):
        line = line.rstrip("\r\n")
        if line and not line.startswith("#"):
            path, _, why = line.partition("|")
            entries[path.strip()] = why.strip()
    return entries


def main():
    tracked = subprocess.run(["git", "-C", ROOT, "ls-files"], capture_output=True, text=True, check=True).stdout.splitlines()
    entries = listed()
    folders = [e for e in entries if e.endswith("/")]
    files, fails, kept, used = 0, 0, 0, set()
    for path in tracked:
        name = path.lower()
        if name.endswith(".py"):
            check = lambda text: python_findings(text, path)
        else:
            match = [(pat, com) for exts, pat, com in PATTERNS if name.endswith(exts)]
            if not match:
                continue
            check = lambda text, m=match[0]: text_findings(text, m[0], m[1])
        full = os.path.join(ROOT, path)
        if not os.path.isfile(full):
            continue
        files += 1
        found = check(open(full, encoding="utf-8", errors="replace").read())
        folder = [f for f in folders if path.startswith(f)]
        if folder:
            used.add(folder[0])
            kept += 1 if found else 0
            continue
        if path in entries:
            used.add(path)
            if found:
                kept += 1
            else:
                print("NO_FLOAT|listed_with_nothing_found=%s|json=0" % path)
                fails += 1
            continue
        if found:
            fails += 1
            print("NO_FLOAT|file=%s|first=%s|count=%d|json=0" % (path, found[0], len(found)))
    for e in sorted(set(entries) - used):
        print("NO_FLOAT|listed_and_gone=%s|json=0" % e)
        fails += 1
    print("NO_FLOAT_DONE|code_files=%d|kept_as_listed=%d|failed=%d|json=0" % (files, kept, fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

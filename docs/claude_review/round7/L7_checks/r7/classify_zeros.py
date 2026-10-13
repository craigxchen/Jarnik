"""Exact classification of the rational zeros of Res found by fastsearch.py (Data)."""
import subprocess, collections
reasons = collections.Counter(); n = 0; genuine = 0
for line in open('fs6_zeros.txt'):
    _, a, b, c, d = line.split()
    out = subprocess.run(['python3', 'verify_cub13.py', a, b, c, d], capture_output=True, text=True).stdout
    n += 1
    if '39 divisors 0): True' in out: genuine += 1
    for l in out.splitlines():
        if l.strip().startswith(('shared', 'class mismatch', 'non-transversal')):
            reasons[l.strip().split(' at ')[0] if 'mismatch' in l else l.strip()] += 1
print('zeros', n, 'genuine', genuine)
for k, v in reasons.most_common(): print(v, k)

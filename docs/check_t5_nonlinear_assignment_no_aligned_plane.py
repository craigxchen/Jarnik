#!/usr/bin/env python3
from itertools import combinations

def dot(a,b): return bin(a & b).count('1') & 1

def main():
    assignment = [28,21,11,30,24,23,13,3,12,7,25,4,31,1,22,3,
                  19,18,10,9,2,16,17,27,26,6,15,5,14,8,20,29]
    assert sorted(assignment) == sorted(list(range(1,32)) + [3])
    seen, hits = set(), 0
    for u in range(1,32):
        for v in range(u+1,32):
            U = tuple(sorted((0,u,v,u^v)))
            for x0 in range(32):
                P = tuple(sorted(x0^w for w in U))
                if P in seen: continue
                seen.add(P)
                us = U[1:]
                for a,b in ((a,b) for a in us for b in us if a != b):
                    vals = {(dot(assignment[x],a),dot(assignment[x],b)) for x in P}
                    if len(vals) == 1 and next(iter(vals)) != (0,0):
                        hits += 1; break
    assert len(seen) == 1240 and hits == 0
    print('PASS: 1240 affine planes checked; no aligned-plane witness')

if __name__ == '__main__': main()

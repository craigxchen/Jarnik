#!/usr/bin/env python3
from itertools import combinations
from math import isqrt

def main():
    cases = squares = 0
    for k in range(1, 5):
        for D in range(2*k, 18):
            for cs in combinations(range(1, D + 1), 2*k):
                for Z in range(D + 1, 4_000):
                    p = 1
                    for c in cs: p *= Z-c
                    if isqrt(p)**2 == p:
                        squares += 1
                        assert Z <= (64*k*D)**(6*k)
                cases += 1
    print(f"PASS: {cases} products and {squares} square values")

if __name__ == '__main__': main()

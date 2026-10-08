#!/usr/bin/env python3
"""Exact rational countermodel: strict obtuseness does not force flip weight."""
from fractions import Fraction
from itertools import combinations

def parity(x): return bin(x).count('1') & 1

def check(t):
    M, b, eps = 1 << t, 2 * (1 << t), Fraction(1, 100)
    rows = range(M); labels = range(1, M)
    # One designated physical copy per row, with the symplectic label Jx.
    def J(x):
        y = 0
        for k in range(0, t, 2):
            y |= ((x >> k) & 1) << (k + 1)
            y |= ((x >> (k + 1)) & 1) << k
        return y
    designated = {(J(x) if x else 1, x) for x in rows}
    columns = []
    for a in labels:
        for k in range(b):
            dr = next((x for aa, x in designated if aa == a and k == x), None)
            columns.append((a, dr, eps if dr is not None else Fraction(1)))
    for x, y in combinations(rows, 2):
        gram = sum(w * (-1) ** (parity(a & x) + parity(a & y)
                              + int(dr == x) + int(dr == y))
                   for a, dr, w in columns)
        assert gram < 0
    F = (M - 1) * eps
    W = sum(w for _, _, w in columns)
    assert F / (W / b) == Fraction(b * (M - 1) * eps, W)
    return F, W

def main():
    for t in (2, 4): check(t)
    print('PASS: exact strict-obtuse weighted symplectic flip models for M=4,16')

if __name__ == '__main__': main()

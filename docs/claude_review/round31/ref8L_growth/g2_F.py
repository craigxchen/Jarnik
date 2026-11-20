#!/usr/bin/env python3
"""Exact count of F_delta(P) truncated at h <= B (VIII.6 integer form).

F counts u = eps*b/conj(b) with b primitive (gcd(b, conj b)=1), norm h supported on P,
h >= (prod P)^gamma (H_0 = 1), 0 < |eps b - conj b|^2 <= h^(1-delta).
Equivalently signed exponent vectors e != 0, b = prod pi_p^{e_p} (conj if e_p<0), h = prod p^|e_p|,
r = arg b reduced mod pi/4 into (-pi/8, pi/8], condition 4 sin^2 r <= h^-delta.

Meet in the middle on the angle; float filter with explicit safety margin, then exact integer
verification of every candidate whose float decision is within the margin.  Every accepted
solution of the final list is ALSO re-verified exactly (cheap enough at these sizes) when
--verify-all is given.

Usage: g2_F.py s logB delta_num/delta_den [P-spec]
"""
import math, sys, itertools
import numpy as np
from fractions import Fraction
from sympy import primerange

def gprime(p):
    for y in range(1, math.isqrt(p) + 1):
        x2 = p - y*y; x = math.isqrt(x2)
        if x*x == x2:
            return (max(x, y), min(x, y))

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gpow(a, e):
    r = (1, 0); base = a
    while e:
        if e & 1: r = gmul(r, base)
        base = gmul(base, base); e >>= 1
    return r

Q4 = math.pi/4

def half_enum(ps, logB):
    """all signed exponent vectors on ps with sum |e| log p <= logB. returns E (int16 matrix), logs, angs"""
    E = np.zeros((1, 0), dtype=np.int16)
    logs = np.zeros(1); angs = np.zeros(1)
    for p in ps:
        x, y = gprime(p)
        phi = math.atan2(y, x); lp = math.log(p)
        amax = int(logB / lp)
        newE = []; newl = []; newa = []
        for a in range(-amax, amax+1):
            nl = logs + abs(a)*lp
            keep = nl <= logB + 1e-12
            if not keep.any(): continue
            col = np.full((int(keep.sum()), 1), a, dtype=np.int16)
            newE.append(np.hstack([E[keep], col])); newl.append(nl[keep]); newa.append(angs[keep] + a*phi)
        E = np.vstack(newE); logs = np.concatenate(newl); angs = np.concatenate(newa)
    return E, logs, np.mod(angs, Q4)

def exact_check(ps, pis, e, dnum, dden, gnum, gden, prodP, B):
    b = (1, 0); h = 1
    for p, pi, a in zip(ps, pis, e):
        if a > 0: b = gmul(b, gpow(pi, a))
        elif a < 0: b = gmul(b, gpow((pi[0], -pi[1]), -a))
        h *= p**abs(int(a))
    x, y = b
    L = min(4*y*y, 4*x*x, 2*(x+y)**2, 2*(x-y)**2)
    if L == 0: return False, h
    # L <= h^(1-delta)  <=>  L^dden <= h^(dden-dnum)
    ok = L**dden <= h**(dden-dnum)
    ok = ok and h**gden >= prodP**gnum and h <= B
    return ok, h

def count(ps, logB, delta, gammas, verify_all=False, tol=1e-9):
    s = len(ps)
    pis = [gprime(p) for p in ps]
    prodP = 1
    for p in ps: prodP *= p
    B = int(math.exp(logB))
    gmin = min(gammas)
    logY0 = gmin*math.log(prodP)
    s1 = (s+1)//2
    E1, l1, a1 = half_enum(ps[:s1], logB)
    E2, l2, a2 = half_enum(ps[s1:], logB)
    sig = float(delta)/2
    # bucket half 2 by floor(log h2)
    bk = np.floor(l2).astype(np.int64)
    sols = []   # list of (e tuple, logh)
    nexact = 0
    for bval in np.unique(bk):
        idx = np.where(bk == bval)[0]
        order = np.argsort(a2[idx]); idx = idx[order]
        A = a2[idx]
        A3 = np.concatenate([A - Q4, A, A + Q4]); I3 = np.concatenate([idx, idx, idx])
        lmin2 = float(bval)
        # admissible e1: l1 + l2 <= logB
        sel = np.where(l1 + lmin2 <= logB + 1e-9)[0]
        if len(sel) == 0: continue
        lh = np.maximum(l1[sel] + lmin2, logY0)
        w = np.arcsin(np.minimum(1.0, 0.5*np.exp(-sig*lh))) + tol
        t = np.mod(-a1[sel], Q4)
        lo = np.searchsorted(A3, t - w, 'left'); hi = np.searchsorted(A3, t + w, 'right')
        cnt = hi - lo
        tot = int(cnt.sum())
        if tot == 0: continue
        if tot > 5e7: raise RuntimeError('too many candidates %d' % tot)
        rep = np.repeat(np.arange(len(sel)), cnt)
        offs = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        j3 = np.repeat(lo, cnt) + offs
        i1 = sel[rep]; i2 = I3[j3]
        pair = np.unique(i1.astype(np.int64) * len(l2) + i2)
        i1 = pair // len(l2); i2 = pair % len(l2)
        lh = l1[i1] + l2[i2]
        r = np.mod(a1[i1] + a2[i2] + Q4/2, Q4) - Q4/2
        thr = np.arcsin(np.minimum(1.0, 0.5*np.exp(-sig*lh)))
        nonzero = ~((np.abs(E1[i1]).sum(1) == 0) & (np.abs(E2[i2]).sum(1) == 0))
        inb = (lh <= logB + 1e-9) & (lh >= logY0 - 1e-9) & nonzero
        sure_in = inb & (np.abs(r) <= thr - tol) & (lh <= logB - 1e-9) & (lh >= logY0 + 1e-9)
        maybe = inb & ~sure_in & (np.abs(r) <= thr + tol)
        for k in np.where(sure_in)[0]:
            e = tuple(int(v) for v in E1[i1[k]]) + tuple(int(v) for v in E2[i2[k]])
            sols.append((e, float(lh[k])))
        for k in np.where(maybe)[0]:
            e = tuple(int(v) for v in E1[i1[k]]) + tuple(int(v) for v in E2[i2[k]])
            nexact += 1
            ok, h = exact_check(ps, pis, e, delta.numerator, delta.denominator,
                                Fraction(gmin).limit_denominator(1000).numerator,
                                Fraction(gmin).limit_denominator(1000).denominator, prodP, B)
            if ok: sols.append((e, float(lh[k])))
    if verify_all:
        bad = 0
        for e, lh in sols:
            ok, h = exact_check(ps, pis, e, delta.numerator, delta.denominator,
                                Fraction(gmin).limit_denominator(1000).numerator,
                                Fraction(gmin).limit_denominator(1000).denominator, prodP, B)
            if not ok: bad += 1
        if bad: print('EXACT VERIFICATION FAILURES:', bad)
    res = {}
    lp = math.log(prodP)
    for g in gammas:
        res[g] = sum(1 for e, lh in sols if lh >= g*lp - 1e-12)
    return res, sols, nexact, len(l1), len(l2)

if __name__ == '__main__':
    s = int(sys.argv[1]); logB = float(sys.argv[2]); delta = Fraction(sys.argv[3])
    ps = [p for p in primerange(5, 10**6) if p % 4 == 1][:s]
    gammas = [0.1, 0.25, 0.4, 0.49]
    res, sols, nex, n1, n2 = count(ps, logB, delta, gammas, verify_all=('-v' in sys.argv))
    print('s=%d P=%s logB=%g delta=%s halves=%d,%d exact_checks=%d' % (s, ps, logB, delta, n1, n2, nex))
    print('  counts by gamma:', res)

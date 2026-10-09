"""Relaxed version used for the restriction argument (Remark 3.4 of construct.md):
is there an injective map pi : F_2^4 minus 0 -> F_2^4 (value 0 allowed) with (x+y).(pi(x)+pi(y)) = 1 for all
distinct x, y?  If not, no Walsh core of order 2^t, t >= 4, admits a skew pairing (restrict to a 4-dim subspace)."""
t = 4; M = 1 << t
els = list(range(1, M)); vals = list(range(M))
dot = [[bin(a & b).count("1") & 1 for b in range(M)] for a in range(M)]
pi = {}; used = set(); nodes = [0]
def bt(i):
    nodes[0] += 1
    if i == len(els):
        return True
    x = els[i]
    for v in vals:
        if v in used:
            continue
        if all(dot[x ^ y][v ^ pi[y]] == 1 for y in els[:i]):
            pi[x] = v; used.add(v)
            if bt(i + 1):
                return True
            used.discard(v); del pi[x]
    return False
print("t=4 relaxed (value 0 allowed): %s, search nodes %d" % ("FOUND " + str(pi) if bt(0) else "no such injection", nodes[0]))

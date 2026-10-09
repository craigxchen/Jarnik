"""Is the Sylvester matrix of order 2^t equivalent to a skew-Hadamard matrix?
Equivalent (derivation in construct.md, Remark 3.4): there is a bijection pi of F_2^t with
(x+y).(pi(x)+pi(y)) = 1 (mod 2) for all x != y.  Backtracking search; exhaustive for t <= 4."""
import sys, itertools
t = int(sys.argv[1]); M = 1 << t
dot = [[bin(a & b).count("1") & 1 for b in range(M)] for a in range(M)]
pi = [-1] * M; used = [False] * M
sols = 0
def bt(x):
    global sols
    if x == M:
        sols += 1
        return True
    for v in range(M):
        if used[v]:
            continue
        ok = True
        for y in range(x):
            if dot[x ^ y][v ^ pi[y]] != 1:
                ok = False; break
        if ok:
            pi[x] = v; used[v] = True
            if bt(x + 1):
                return True
            used[v] = False
    pi[x] = -1
    return False
found = bt(0)
print("t=%d M=%d: skew-equivalent bijection %s" % (t, M, "FOUND " + str(pi) if found else "does NOT exist (exhaustive)"))

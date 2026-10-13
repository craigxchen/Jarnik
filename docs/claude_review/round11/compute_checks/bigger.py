"""rho vs zeta for larger n (d = 5..7) and small T: looking for the first faithful slab-boxes (zeta < rho).
Profiles: Sylvester-Hadamard (n = 8), all balanced patterns (n = 6), random sign profiles."""
import sys, itertools, random, time
from core import build_A, width
from zono import zeta
from scan import rho_p
from models import fullrank, canon, random_profile

def hadamard8():
    H = [[1]]
    for _ in range(3):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    pats = []
    for r in H[1:]:
        S = tuple(x for x in range(8) if r[x] == 1)
        pats.append(canon(8, S))
    return tuple(sorted(set(pats)))

def balanced(n):
    return tuple(sorted(set(canon(n, S) for S in itertools.combinations(range(n), n // 2))))

if __name__ == "__main__":
    which = sys.argv[1]
    T = int(sys.argv[2]); L = int(sys.argv[3])
    if which == "had8":
        n, prof = 8, hadamard8()
    elif which.startswith("pairs"):
        n = int(which[5:]); prof = tuple(sorted(set(canon(n, S) for S in itertools.combinations(range(n), 2))))
    elif which.startswith("bal"):
        n = int(which[3:]); prof = balanced(n)
    else:
        n = int(which[4:]); rng = random.Random(int(sys.argv[4]))
        while True:
            prof = random_profile(n, rng)
            if fullrank(n, prof):
                break
    t0 = time.time()
    A = build_A(n, T, [(S, L) for S in prof])
    rmax = max(width(A, S) for S in prof)
    r = rho_p(A)
    t1 = time.time()
    print("  |A|=%d rmax=%d rho_p=%d (%.1fs)" % (len(A), rmax, r, t1 - t0), flush=True)
    if len(sys.argv) > 5 and sys.argv[5] == "norho":
        sys.exit(0)
    z, Z = zeta(n, T, prof, L, r)
    t2 = time.time()
    print(which, "n=%d T=%d L=%d prof=%s |A|=%d box=%d rmax=%d rho_p=%d zeta=%d  (%.1fs, %.1fs)" % (
        n, T, L, prof, len(A), (2*T+1)**(n-1), rmax, r, z, t1 - t0, t2 - t1), flush=True)
    print("  zonotope:", Z)

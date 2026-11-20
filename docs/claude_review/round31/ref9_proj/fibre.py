import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from gz import nu_formula, nu_brute
from sympy import primerange
tot = 0; bad = 0; mx = 0
for p in primerange(3, 252):
    for D in range(p):
        tot += 1
        f, b = nu_formula(p, D), nu_brute(p, D)
        if f != b: bad += 1
        if D % p:
            mx = max(mx, b - (p+3)/2)
# also prime powers? not needed. p=2:
print("pairs (p,D) checked", tot, "mismatches", bad, "max nu_p-(p+3)/2 over p not dividing D", mx)
# also check the point count p - chi4(p) on the fibre
for p in primerange(3, 60):
    chi4 = 1 if p % 4 == 1 else -1
    for D in range(1, p):
        cnt = sum((x*x+y*y-D) % p == 0 for x in range(p) for y in range(p))
        assert cnt == p - chi4
print("fibre point counts p - chi4(p) verified for p<60")

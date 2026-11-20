# Judge: the 10-point record cluster and K_C values (imports the judge's own viii1_test helpers only).
import sys, math
sys.argv = ['x', '10']
exec(open('viii1_test.py').read().split("for N in range(5, NMAX+1):")[0])
def KC(C, sharp):
    return sum(Sset(C, sharp)) if C > (math.sqrt(2) if sharp else 1) else 0
for C in (1, 1.3, math.sqrt(2)+1e-9, 1.7, 2, 3, 4, 5.9, 8, 8.2394, 10):
    print('C=%.4f  K_C Astra %d  sharpened %d' % (C, KC(C, False), KC(C, True)))
N = 1176852625
pts = lattice_points(N); n = len(pts); R4 = N**0.25
ang = [math.atan2(z[1], z[0]) % (2*math.pi) for z in pts]
best = min(range(n), key=lambda s: (ang[(s+9) % n]-ang[s]) % (2*math.pi))
Z = [pts[(best+k) % n] for k in range(10)]
g = Z[0]
for z in Z[1:]: g = ggcd(g, z)
Zp = [gdivexact(z, g) for z in Z]; f = factor(N // gnorm(g)); primes = sorted(f)
rows = [[gval(z, pi_of(p)) for p in primes] for z in Zp]
print('points', n, 'arc C', (ang[(best+9) % n]-ang[best]) % (2*math.pi)*R4, 'gcd norm', gnorm(g), 'N\'=', f)
Cp = max(math.sqrt(gnorm(gsub(a, b))) for a in Zp for b in Zp)/(N//gnorm(g))**0.25
print('chord C\'', Cp, 'affine rank', affine_rank(rows))
enl = [list(r) for r in rows]
for l in (3, 7, 11):
    cl = sorted(set((z[0] % l, z[1] % l) for z in Zp)); idx = {c: t for t, c in enumerate(cl)}
    for i, z in enumerate(Zp):
        oh = [0]*len(cl); oh[idx[(z[0] % l, z[1] % l)]] = 1; enl[i] += oh
    print(' + block mod', l, ': affine rank', affine_rank(enl), '(classes hit %d)' % len(cl))
# z mod 2 block
e2 = [list(r) + [ (z[0] % 2) ] for r, z in zip(rows, Zp)]
print(' rows + z mod 2 indicator: affine rank', affine_rank(e2))

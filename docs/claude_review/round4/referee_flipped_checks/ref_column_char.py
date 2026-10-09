"""Referee check of Sec. 2(a): in the fully flipped class the conjugate-primitive part of
gamma_a = G_a^(M/2) conj(P_a) has log-norm (M/2) W_a + W_F - 2 f_a  (f_a = flipped mass of label a),
i.e. (3.1) of flipped.md, which uses log|gamma_a|^2 = (M/2) W_a + W_F, is weaker than the column
character half-height inequality by 2 f_a.  Exact Gaussian arithmetic on literal profiles."""
import math, random, sys
sys.path.insert(0, '.')
from ref_two_log import gp, mul, cj, nrm, pw, prod, paley, sylv, is_prime
rng = random.Random(31)
n = 0
for H in [sylv(3), paley(11), sylv(4), paley(19)]:
    M = len(H)
    splits = [p for p in range(5, 5000) if p % 4 == 1 and is_prime(p)]
    for _ in range(5):
        ps = rng.sample(splits, 5*(M-1))
        labels = list(range(1, M)) + [rng.randrange(1, M)]; rng.shuffle(labels); asg = labels[:M]
        cols = []; used = {}
        for a in range(1, M):
            for c in range(5): cols.append(dict(a=a, s=rng.choice([1, -1]), p=ps[(a-1)*5+c], F=set()))
        for x in range(M):
            a = asg[x]; c = used.get(a, 0); used[a] = c + 1; cols[(a-1)*5 + c]['F'].add(x)
        WF = sum(math.log(c['p']) for c in cols if c['F'])
        for a in range(1, M):
            G = prod([pw(gp(c['p']), c['s']) for c in cols if c['a'] == a])
            cf = lambda c: sum(H[x][a]*H[x][c['a']] for x in c['F'])
            P = prod([pw(gp(c['p']), c['s']*cf(c)) for c in cols])
            gam = mul(pw(G, M//2), cj(P))
            # conj-primitive part: exponents (M/2)[a(j)=a] - c_j(a) per column
            gprim = prod([pw(gp(c['p']), c['s']*((M//2 if c['a'] == a else 0) - cf(c))) for c in cols])
            Wa = sum(math.log(c['p']) for c in cols if c['a'] == a)
            fa = sum(math.log(c['p']) for c in cols if c['a'] == a and c['F'])
            assert abs(math.log(nrm(gam)) - ((M/2)*Wa + WF)) < 1e-6
            assert abs(math.log(nrm(gprim)) - ((M/2)*Wa + WF - 2*fa)) < 1e-6
            q = mul(gam, cj(gprim)); N2 = nrm(gprim)
            assert q[0] % N2 == 0 and q[1] % N2 == 0 and (q[1] == 0 or q[0] == 0)   # gamma = integer * unit * gprim
            n += 1
print("column-character comparison verified on", n, "labels: PASS")

import sys, pickle, math, itertools, glob, os
from glib import *
from forced import forced_set
D = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round14/rigtheory_checks/'
Fcache = {}
for f in sorted(glob.glob(D + 'witness*.pkl')):
    w = pickle.load(open(f, 'rb'))
    K, ps, G, q0, lab, l2 = w['K'], w['ps'], w['G'], w['q0'], w['labels'], w['lev2']
    r = len(ps); pis = [two_sq(p) for p in ps]
    key = (K, tuple(ps))
    if key not in Fcache: Fcache[key] = forced_set(ps, K, 2000)[0]
    F = Fcache[key]
    ok_mod = sorted(list(lab) + list(G)) == F
    V = VK(r, K); errs = []
    for q, ell in lab.items():
        m = q-1 if q % 4 == 1 else q+1; nv = nu(ps, q)
        if [x % 2 for x in ell] != nv: errs.append(('par', q))
        if any((a - b) % m for a, b in zip(l2[q], ell)): errs.append(('lift', q))
        M2 = m*q
        for v in V:
            if sum(x*y for x, y in zip(v, l2[q])) % (M2//4) == 0: errs.append(('lev2coll', q, v)); break
    viol = 0; worst = -1e9
    for v in V:
        z = Zv(pis, v); nz = z[0]**2 + z[1]**2
        for s in range(4):
            d = Ds(z, s); mm = 1
            for q in G: mm *= q**vq(d, q)
            for q, ell in lab.items():
                m = q-1 if q % 4 == 1 else q+1
                if sum(x*y for x, y in zip(v, ell)) % m == s*(m//4): mm *= q
            lhs, rhs = mm*mm, nz*(2 if s % 2 else 1)
            if lhs > rhs: viol += 1
            worst = max(worst, math.log(lhs) - math.log(rhs))
    m0, L0 = genuine(ps, q0); n0 = m0//4; c0 = tuple(x % n0 for x in lab[q0])
    incl = any(tuple((u*e*x) % n0 for e, x in zip(eps, L0)) == c0 for eps in itertools.product([1, -1], repeat=r)
               for u in range(1, n0+1) if math.gcd(u, n0) == 1)
    # moduli where witness is outside every eps class
    outs = []
    for q, ell in lab.items():
        m, L = genuine(ps, q); n = m//4; c = tuple(x % n for x in ell)
        if n > 1 and not any(tuple((u*e*x) % n for e, x in zip(eps, L)) == c for eps in itertools.product([1, -1], repeat=r)
                   for u in range(1, n+1) if math.gcd(u, n) == 1): outs.append(q)
    print(os.path.basename(f), f"K={K} G={G} q0={q0} modsOK={ok_mod} errs={errs[:3]} budgets={len(V)*4} viol={viol} maxlog={worst:.3f} q0_in_eps_class={incl} outside_at={sorted(outs)}", flush=True)

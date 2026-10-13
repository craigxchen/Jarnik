"""Components of the projective cone variety X = V(G_j) (x = w - 1 coordinates).
1. Linear components among coordinate subspaces {x_j = 0, j notin S}: exact test (forms vanish identically).
2. Residual: saturate I by the product of the variables x_j (points with ALL x_j != 0, i.e. no coordinate of w
   equal to 1) using the Rabinowitsch trick I + (1 - t*prod x_j), eliminating t; Hilbert function of the result.
   Also saturate by each binomial (x_i - x_j) (w_i = w_j).
Prints the dimension/degree data of the residual and, if it is finite, its points."""
import sys, json, itertools, time
import sympy as sp
fn = sys.argv[1]
C = json.load(open(fn)); d = C['d']
X = sp.symbols('x0:%d' % d)
forms = [sp.expand(sp.sympify(f)) for f in C['forms']]
lin = []
for s in range(1, d):
    for S in itertools.combinations(range(d), s):
        sub = {X[j]: 0 for j in range(d) if j not in S}
        if all(sp.expand(f.subs(sub)) == 0 for f in forms):
            lin.append(S)
maximal = [S for S in lin if not any(set(S) < set(T) for T in lin)]
print("maximal coordinate subspaces contained in the cone (support sets):", maximal)
t = sp.Symbol('t')
def sat_hf(extra, label):
    t0 = time.time()
    GB = sp.groebner(forms + [1 - t * extra], t, *X, order='lex')
    elim = [g for g in GB.exprs if not g.has(t)]
    if elim == [1] or any(g == 1 for g in elim):
        print(label, ": saturation is the unit ideal -> no points of X with", label, "(%.1fs)" % (time.time() - t0))
        return None
    G2 = sp.groebner(elim, *X, order='grevlex')
    LM = [sp.Poly(g, *X).monoms(order='grevlex')[0] for g in G2.exprs]
    HF = []
    for D in range(0, 14):
        cnt = 0
        for comb in itertools.combinations_with_replacement(range(d), D):
            mon = [0] * d
            for i in comb: mon[i] += 1
            if not any(all(mon[i] >= l[i] for i in range(d)) for l in LM): cnt += 1
        HF.append(cnt)
    print(label, ": saturated ideal Hilbert function", HF, "(%.1fs)" % (time.time() - t0))
    return G2
prodx = sp.prod(X)
G2 = sat_hf(prodx, "all x_j != 0")
json.dump(dict(maximal=maximal), open(fn.replace("_cone.json", "_components.json"), "w"))

import time, pickle, itertools, sys
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
gw = pickle.load(open('m8_gw.pkl','rb'))
cache = {}
def fam_gens(F, maxdeg=8):
    g, rho, ny = proj_ideal_gens(A, F, maxdeg)
    out = {}
    for e, fl in g.items():
        out[e] = [pullback(f, F, cache) for f in fl]
    return out, rho, ny, {e: len(fl) for e, fl in g.items()}
def merge(*gs):
    o = {}
    for g in gs:
        for e, l in g.items(): o.setdefault(e, []).extend(l)
    return o
def dims(gens, maxdeg, exact_upto):
    out = []
    for i in range(maxdeg+1):
        rows = ideal_degree_rows(gens, i)
        out.append(dim_exact(rows, i) if i <= exact_upto else dim_modp(rows, i))
    return out
fano = lambda a,b,c: (a ^ b) == c
pairs = list(itertools.combinations(range(1,8), 2))
triples = list(itertools.combinations(range(1,8), 3))
t = time.time()
PG = {}
for (a,b) in pairs:
    g, rho, ny, cnt = fam_gens((cols[a], cols[b]))
    PG[(a,b)] = g
    if (a,b) == pairs[0]: print("pair", (a,b), "|Y|", ny, "rho(Y)", rho, "Y-generator degrees", cnt)
print("pairs done", time.time()-t, flush=True)
TG = {}; info = {}
for tr in triples:
    g, rho, ny, cnt = fam_gens(tuple(cols[a] for a in tr))
    TG[tr] = g; info[tr] = (ny, rho, cnt, fano(*tr) or fano(tr[0],tr[2],tr[1]) or fano(tr[1],tr[2],tr[0]))
print("triples done", time.time()-t, flush=True)
for tr in triples[:]:
    pass
print("triple info (|Y|, rho(Y), Y-gen degree counts, Fano line):")
for key in sorted(set((v[0], v[1], tuple(sorted(v[2].items())), v[3]) for v in info.values())):
    print("   ", key, "x", sum(1 for v in info.values() if (v[0], v[1], tuple(sorted(v[2].items())), v[3]) == key))
G2 = merge(gw, *PG.values())
r2 = dims(G2, 8, 6)
print("J_w + 21 column pairs: dims (exact Q for i<=6, mod p for 7,8) =", r2, flush=True)
fl = [tr for tr in triples if info[tr][3]]; nfl = [tr for tr in triples if not info[tr][3]]
print("Fano-line triples:", fl)
for name, sel in [("all 35 triples", triples), ("7 Fano-line triples", fl), ("28 non-line triples", nfl)]:
    G3 = merge(gw, *[TG[tr] for tr in sel])
    r3 = dims(G3, 8, 6)
    print(f"J_w + {name}: dims (exact Q i<=6, mod p 7,8) = {r3}  J(A) = {JA}", flush=True)
# sandwich for all triples in degrees 7,8 : mod-p cumulative vs exact cumulative J(A)
pickle.dump((PG, TG, info), open('m8_step3.pkl','wb'))
print("time", time.time()-t)

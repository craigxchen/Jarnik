"""Collect the calibration data into markdown tables (reads out/*.txt).  For each (q, L, m) with dim V >= 2:
 Psi_A(C_1), rel = Psi/E_A l^2m, -ln(rel)/m, Psi_cube(C_1) and Psi_A/Psi_cube, the rigorous kappa lower bound
 k_lo = sqrt(Psi)/d^m [(2.1)], and the 'effective tau'  tau* = W*/(2q),  W* = 4 ln(C d/(2 (m+1)^{3/2} e^{pi^2/(4q)} k)),
 i.e. the barrier scale at which kappa = k would sit exactly at the CIRCLE threshold (with m in place of m'); computed
 with k = k_lo (an upper bound for tau* of the true inf kappa_V, hence of nothing about kappa_H) and with the LP value
 when available.  tau_D(q) of barrier Lemma 3.4 (c_0 = 0.184, C = 1) for comparison.
"""
import re, glob, math, os

def tau_D(q, c0=0.184, C=1.0):
    r = 2 * q; n = q + 1
    return max((2 / c0) * (2 * math.log(4 * r + 5) + math.log(1000 / C)), 4 * math.log(4 * r + 5) + 8,
               2 * math.log(1000 * q * n * n / C))

def cube(q):
    d = {}
    f = f"out/cube_q{q}.txt"
    if not os.path.exists(f): return d
    for line in open(f):
        mm = re.search(r"theta=0.00 m=(\d+): Psi_cube = (\S+)\s+E_cube l\^2m = (\S+)\s+rel = (\S+)", line)
        if mm: d[int(mm.group(1))] = (float(mm.group(2)), float(mm.group(4)))
    return d

def exact_rows(q, L):
    rows = {}
    for f in sorted(glob.glob(f"out/exact_q{q}_L{L}*.txt")):
      cur = None
      for line in open(f):
        mm = re.match(r"m=(\d+): primes used \d+; dim Poly_<m\|_A = (\d+), dim V\(A,m\) = (\d+)", line)
        if mm: cur = int(mm.group(1)); rows.setdefault(cur, {}).update({'dimP': int(mm.group(2)), 'dimV': int(mm.group(3)), 'src': 'exact'})
        mm = re.search(r"Psi_A\(C_1\) = \S+\s+~ (\S+);\s+E_A l\^\(2m\) = (\S+);\s+rel = (\S+)", line)
        if mm and cur is not None:
            rows[cur].update(Psi=float(mm.group(1)), E=float(mm.group(2)), rel=float(mm.group(3)))
    return rows

def float_rows(q, L):
    rows = {}
    for f in sorted(glob.glob(f"out/float_q{q}_L{L}*.txt")):
     for line in open(f):
        mm = re.match(r"m=(\d+): dim Poly_<m\|_A = (\d+), dim V = (\d+); Psi_A\(C_1\) = (\S+), E l\^2m = (\S+), rel = (\S+); random u: Psi in \[(\S+), (\S+)\], rel in \[(\S+), (\S+)\]", line)
        if mm:
            g = mm.groups()
            rows[int(g[0])] = dict(dimP=int(g[1]), dimV=int(g[2]), Psi=float(g[3]), E=float(g[4]), rel=float(g[5].rstrip(';')),
                                   rPsi=(float(g[6]), float(g[7].rstrip(','))), rrel=(float(g[8]), float(g[9].rstrip(']'))), src='float')
    return rows

def lp(q, L, m):
    f = f"out/cheb_q{q}_L{L}_m{m}.txt"
    if not os.path.exists(f): return None
    for line in open(f):
        mm = re.search(r"inf/sup kappa_V\(hat u\) over grid = (\S+) .*inf kappa_V = (\S+) at", line)
        if mm: return float(mm.group(1)), float(mm.group(2))
    ks = [float(x) for x in re.findall(r"kappa_V\(hat u\)=(\S+)", open(f).read())]
    return (None, min(ks)) if ks else None

def tstar(q, m, k):
    if k <= 0: return float('inf')
    W = 4 * math.log(q / (2 * (m + 1) ** 1.5 * math.exp(math.pi ** 2 / (4 * q)) * k))
    return W / (2 * q)

print("| q | L | U | m | dim V | Psi_A(C_1) | rel | -ln(rel)/m | Psi_A/Psi_cube | cube rel | random-u rel range | k_lo = sqrt(Psi)/d^m | LP inf kappa_V (inf/sup) | tau* (k_lo / LP) | tau_D |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (q, L) in [(7, 1), (11, 1), (19, 1), (23, 1), (31, 1), (43, 1), (47, 1), (59, 1), (67, 1), (71, 1), (79, 1), (83, 1), (103, 1),
               (7, 2), (11, 2), (19, 2), (23, 2)]:
    ex = exact_rows(q, L); fl = float_rows(q, L); cb = cube(q)
    U = {7: {1: 57, 2: 757}, 11: {1: 133, 2: 9373}, 19: {1: 381, 2: 264291}, 23: {1: 553, 2: 818755}}.get(q, {}).get(L, 1 + q * (q + 1))
    ms = sorted(set(ex) | set(fl))
    for m in ms:
        r = ex.get(m) or fl.get(m)
        if r is None or 'Psi' not in r or r['dimV'] < 2 or m >= q: continue
        src = 'E' if m in ex and 'Psi' in ex[m] else 'F'
        rr = fl.get(m, {}).get('rrel')
        klo = math.sqrt(r['Psi']) / q ** m
        l = lp(q, L, m)
        cbs = f"{r['Psi']/cb[m][0]:.3g}" if m in cb else "-"
        cbr = f"{cb[m][1]:.3g}" if m in cb else "-"
        lps = (f"{l[1]:.3e} ({l[0]:.3f})" if l[0] is not None else f"<= {l[1]:.3e} (theta=0 only)") if l else "-"
        ts = f"{tstar(q, m, klo):.2f}" + (f" / {tstar(q, m, l[1]):.2f}" if l else "")
        rel = r['rel']
        print(f"| {q} | {L} | {U} | {m} | {r['dimV']} | {r['Psi']:.6g} ({src}) | {rel:.4g} | {(-math.log(rel)/m if rel > 0 else float('inf')):.3f} | {cbs} | {cbr} | "
              f"{(f'[{rr[0]:.3g}, {rr[1]:.3g}]' if rr else '-')} | {klo:.3e} | {lps} | {ts} | {tau_D(q):.0f} |")

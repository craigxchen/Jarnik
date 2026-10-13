import sympy as sp, pickle, sys
from genus import monodromy_genus
d = pickle.load(open(sys.argv[1], 'rb'))
G = d['G']; a, b = [sp.Symbol(s) for s in d['free']]
for (x, y) in ((a, b), (b, a)):
    g, n, m = monodromy_genus(G, x, y, verbose=False)
    print('project to', x, ': sheets', n, 'critical pts', m, 'genus', g, flush=True)

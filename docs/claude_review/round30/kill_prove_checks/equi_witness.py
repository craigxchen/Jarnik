# exact witnesses: for each 6-point equidistant class, a weight-2 integer vector in span_Q(differences)
import itertools, numpy as np, sympy as sp
exec(open('equi_lattice.py').read().split("# short integer vectors")[0].replace("M=int(sys.argv[1])","M=6"))
for key,c in classes.items():
    pts=[ones]+[W[j] for j in c]
    Dm=sp.Matrix([list((p-ones)//2) for p in pts[1:]])
    print('code rows:'); print(np.array(pts))
    for k,l in itertools.combinations(range(D),2):
        for s in (1,-1):
            a=[0]*D; a[k]=1; a[l]=s
            # a in rowspace iff rank unchanged
            if Dm.rank()==Dm.col_join(sp.Matrix([a])).rank():
                x=sp.Matrix([a]).T
                sol=Dm.T.gauss_jordan_solve(x)[0]
                print('  weight-2 vector',a,'= combination of differences with coeffs',list(sol.T.subs({s_:0 for s_ in sol.free_symbols})))
                break
        else: continue
        break

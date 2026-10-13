s = open('qsearch.py').read()
old = s[s.index("    for f, m in sp.factor_list(n)[1]:"):s.index("print('done, found', found)")]
new = '''    for f, m in sp.factor_list(n)[1]:
        if sp.degree(f, var) < 1: continue
        for r in sp.roots(sp.Poly(f, var), filter='Q'):
            p = [complex(sub[x]) if x != var else complex(r) for x in SY]
            try:
                c, res, tv = curve(p); sm = summary(c)
                ok = sm['nplanes'] == 13 and not sm['bad'] and sm['ncontacts'] == 19 and sm['mindist'] > 1e-6
            except Exception:
                ok = False
            if ok:
                print('GENUINE RATIONAL POINT', {str(k): str(v) for k, v in sub.items()}, str(var), '=', r, flush=True)
                found += 1
'''
s = s.replace(old, new)
s = s.replace("import sympy as sp, pickle, random, sys, itertools\n", "import sympy as sp, pickle, random, sys, itertools, warnings\nwarnings.filterwarnings('ignore')\nfrom kapcurve import curve\nfrom kapcheck import summary\n")
open('qsearch.py', 'w').write(s)

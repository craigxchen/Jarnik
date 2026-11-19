# check that F = prod_{0<=a<b<=5}(y_a-y_b) (y_0=0) satisfies all j=3 conditions exactly (over Z)
import sys; sys.argv=['x','3','1000003']
exec(open('h0M.py').read().split("for wp in (True, False):")[0])
import flint
R = flint.fmpz_mpoly_ctx.get(('y1','y2','y3','y4','y5'),'lex') if hasattr(flint,'fmpz_mpoly_ctx') else None
ys = R.gens(); Y=[ys[0]-ys[0]]+list(ys)
F = ys[0]**0
for a in range(6):
    for b in range(a+1,6): F *= (Y[a]-Y[b])
coef = {tuple(e):int(c) for e,c in zip(F.monoms(),F.coeffs())}
vec = [coef.get(m,0) for m in M]
rows = build(True)
bad = sum(1 for r in rows if sum(v*vec[k] for k,v in r.items())!=0)
print("terms of F:",len(coef),"rows violated by F (exact over Z):",bad)

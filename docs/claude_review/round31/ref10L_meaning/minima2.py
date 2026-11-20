import math
exec(open('minima.py').read().split('random.seed(1)')[0])
a,b=83,59; u=a*(a+b); v=b*(a+2*b)
y=[0,u,v,u+v,(a+b)*(a+2*b)]
pr=[59,67,71,83]
for p in pr: print(p, sorted({t%p for t in y}))
rows,q=basis_Nprime(y,pr)
red,gsn=lll(rows)
spread=max(y)-min(y); R=math.sqrt(len(y)-1)*spread
print("q",q,"spread",spread,"C",spread/math.sqrt(q),"GS",[round(g) for g in gsn],"R",round(R))
print("reduced basis",red)

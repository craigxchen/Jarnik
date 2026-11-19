from itertools import combinations
from sympy import Matrix, Rational, symbols, factor

H=tuple(map(Rational,(0,1,3,6,10)))
marks=range(5); edges=list(combinations(marks,2))
def edge(i,j):return tuple(sorted((i,j)))
def make(h,s):
 D=Matrix(5,5,lambda i,j: 1/(h[j]-h[i]) if i!=j else -sum(1/(h[k]-h[i]) for k in marks if k!=i))
 ss=Matrix([sum(c*x**j for j,c in enumerate(s)) for x in h]); z=Matrix.diag(*h)*D*ss-2*ss
 d=D[:,1:].gauss_jordan_solve(z)[0];d=Matrix([0]+list(d))
 return {e:(d[j]-d[i]-h[i]*ss[j]+h[j]*ss[i])/(2*(h[j]-h[i])) for e in edges for i,j in [e]}
def matrix(h,x):
 rr=[]
 for i in range(1,5):
  rr.append([int(i in e)-int(0 in e) for e in edges])
  rr.append([(int(i in e)-int(0 in e))*x[e] for e in edges])
 for tr in combinations(marks,3):
  l,m=set(marks)-set(tr)
  v={e:0 for e in edges}
  for j in range(3):
   a,b,c=tr[j],tr[(j+1)%3],tr[(j+2)%3]
   e,f,g=edge(a,b),edge(c,l),edge(c,m)
   k=h[b]-h[a]
   v[e]+=k*(x[f]+x[g]);v[f]-=k*(x[e]+x[g]);v[g]-=k*(x[e]+x[f])
  rr.append(list(v.values()))
 return Matrix(rr)
if __name__=='__main__':
 for s in [(0,1,2,3),(0,0,1,0),(0,0,0,1),(0,1,0,0)]:
  x=make(H,s); L=matrix(H,x); print(s,L.rank(),L.nullspace()[:1])
def re_quadratic(h,x,Y):
 rr=[]
 for tr in combinations(marks,3):
  l,m=set(marks)-set(tr);v=[0,0,0,0]
  for a,b,c in [tr,tr[1:]+tr[:1],tr[2:]+tr[:2]]:
   e,f,g=edge(a,b),edge(c,l),edge(c,m); k=h[b]-h[a]
   e0,e1=Y[edges.index(e),:]; f0,f1=Y[edges.index(f),:];g0,g1=Y[edges.index(g),:]
   v[0]+=k*(x[e]*x[f]+x[e]*x[g]+x[f]*x[g])
   v[1]+=k*(e0*f0+e0*g0-f0*g0)
   v[2]+=k*(e0*f1+e1*f0+e0*g1+e1*g0-f0*g1-f1*g0)
   v[3]+=k*(e1*f1+e1*g1-f1*g1)
  rr.append(v)
 return Matrix(rr)
if __name__=='__main__':
 for s in [(0,1,2,3),(0,0,1,0),(0,0,0,1),(0,1,0,0)]:
  x=make(H,s); L=matrix(H,x);Q=re_quadratic(H,x,Matrix.hstack(*L.nullspace()))
  print('Q',s,Q.rank(),Q[:,1:].rank(),Q.nullspace())

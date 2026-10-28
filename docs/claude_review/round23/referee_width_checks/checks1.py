import sys; sys.path.insert(0,sys.argv[1])
from regp import T, regp, full_at
from rec_check_funcs import F
# 1. slices (Lemma 2.1) by enumeration
ok=True
for n in range(2,6):
  for a in range(5):
    for b in range(5):
      S=set(T(n,a,b))
      for c in range(-b-2,a+3):
        sl=sorted(tuple(s[1:]) for s in S if s[0]==c)
        if 0<=c<=a: exp=sorted(T(n-1,a-c,b))
        elif -b<=c<0: exp=sorted(T(n-1,a,b+c))
        else: exp=[]
        if sl!=exp: ok=False; print("slice mismatch",n,a,b,c)
print("Lemma 2.1 slices (n<=5,a,b<=4):",ok)
# S_m = mR_5 ∩ M sizes
print("|T_5(m,m)| m=1..5:",[len(T(5,m,m)) for m in range(1,6)])
# 2. reg T_n(a,b) vs F_n (mod p upper bound) and deficiency at F_n-1 (evidence)
for n in (3,4,5):
  lim={3:6,4:5,5:3}[n]
  rows=[]
  for a in range(lim+1):
    for b in range(a,lim+1):
      S=T(n,a,b); f=F(n,a,b)
      up=full_at(S,f); below=full_at(S,f-1) if f>0 else False
      rows.append((a,b,f,up,below))
  print("n=%d: all full at F_n:"%n, all(r[3] for r in rows), " full already at F_n-1 at:",[(r[0],r[1]) for r in rows if r[4]])
print("T_5(2,3): full at 7?",full_at(T(5,2,3),7),"full at 6?",full_at(T(5,2,3),6))
# 3. reg(S_m) for m<=3 directly
for m in (1,2,3):
  S=T(5,m,m); print("m=%d |S|=%d  full at floor(10m/3)=%d: %s ; full at one less: %s"%(m,len(S),10*m//3,full_at(S,10*m//3),full_at(S,10*m//3-1)))

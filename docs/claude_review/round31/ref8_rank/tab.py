import pickle, sys
from math import sqrt, log
from kcfun import Kinit
recs=[]
for f in sys.argv[1:]: recs += pickle.load(open(f,'rb'))
print("windows", len(recs))
for Nmin in [1, 10**4]:
  print(f"-- N >= {Nmin}; thresholds on ARC constant C (theorem hypothesis)")
  print(" T   K_C(T)  maxM  max(M-1-omega)  ex  max(M-1-rank)  ex")
  for T in [0.5,1,sqrt(2),1.7,2,2.5,3,4,5,5.9,6,8,8.3,9]:
    S=[r for r in recs if r[4]<=T and r[0]>=Nmin]
    if not S: continue
    a=max(S,key=lambda r:(r[1]-1-r[2], -r[4])); b=max(S,key=lambda r:(r[1]-1-r[3], -r[4]))
    print(f"{T:.3f} {Kinit(log(T))} {max(r[1] for r in S)} {a[1]-1-a[2]} {a[:4]+(round(a[4],4),)} {b[1]-1-b[3]} {b[:4]+(round(b[4],4),)}")
# smallest constants at which defects occur
for name,idx in [("omega",2),("rank",3)]:
  for d in [1,2,3,4,5]:
    S=[r for r in recs if r[1]-1-r[idx]>=d and r[0]>=10**4]
    if S:
      r=min(S,key=lambda r:r[4]); r2=min(S,key=lambda r:r[5])
      print(f"defect_{name}>={d} (N>=1e4): min arc C {r[4]:.5f} at {r[:4]}; min chord C' {r2[5]:.5f} at {r2[:4]}")

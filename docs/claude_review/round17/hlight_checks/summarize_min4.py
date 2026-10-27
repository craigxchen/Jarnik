# Summarise exact 4-term minima: check against Lemma 2.7 / Prop 2.5 bounds, decide when the global minimum is exact.
import re, math
from fractions import Fraction as Fr
rows=[]
for fn in ['min4_small_out.txt','min4_regime_out.txt']:
    for line in open(fn):
        m=re.match(r'q=(\d+) configs=\d+ min sigma\(n=4,F=0\)=(-?\d+)',line)
        if m: rows.append((int(m.group(1)),int(m.group(2))))
bad=0; exact=0; minratio=(9,None); lowq=[]
for q,s in rows:
    lb4=max(Fr(q,3), Fr(q,2)-11-Fr(3,2)*(math.isqrt(q)+1))
    if s<lb4: bad+=1
    if s<Fr(q-3,2): exact+=1   # then global min over all non-pair lambda = s (Thm 2.8 classes)
    else: lowq.append(q)
    if q>=19 and s/(2*q)<minratio[0]: minratio=(s/(2*q),q)
print(f"primes={len(rows)} (q from {rows[0][0]} to {rows[-1][0]}); violations of Lemma 2.7/Prop2.5 n=4 bound: {bad}")
print(f"global min exactly determined (4-term min < (q-3)/2) for {exact} primes; otherwise (global min in [(q-3)/2,(q-1)/2]): {lowq}")
print(f"smallest 4-term sigma/r over q>=19: {minratio[0]:.4f} at q={minratio[1]}")
for q,s in rows:
    if q in (3251,3259,3271,5479) or q<=47 or q in (103,199,499,1019,1091):
        print(f"  q={q} min4 sigma={s} sigma/r={s/(2*q):.4f} (q-3)/2={(q-3)//2} LB(q)={float(min(max(Fr(q,3),Fr(q,2)-11-Fr(3,2)*(math.isqrt(q)+1)),Fr(q-3,2))):.1f}")

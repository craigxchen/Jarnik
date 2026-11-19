exec(open('filter_ref.py').read().split('out=[]')[0])
from collections import Counter
c=Counter(); mx=Counter()
for S in reps:
    rk,M,ns=kernel(S)
    if rk!=5: continue
    w=find_witness(sp.Matrix.hstack(*ns))
    if w: c[w[0]]+=1; mx[max(abs(x) for x in w[1])]+=1
print(c,mx)

# Paley profile a_ij = 1 + chi(j-i), e_j=2 on columns j != 0, rows i in S subset of A={i: chi(-i)=+1}.
# Check: affine rank = |S|-1 (VIII.1 with K=0), omega = q-1, and pair margin d_ij - W/2 in units of log P (equal primes).
import random
from analyse import rank_mod, affrank
for q in [43,127,383]:
    chi=lambda x:0 if x%q==0 else (1 if pow(x%q,(q-1)//2,q)==1 else -1)
    A=[i for i in range(q) if chi(-i)==1]
    for rho in [1.0,0.9]:
        S=[i for i in A if random.random()<rho]
        rows=[[1+chi(j-i) for j in range(1,q)] for i in S]
        ar=affrank(rows)
        dist=set(sum(abs(x-y) for x,y in zip(rows[a],rows[b])) for a in range(len(S)) for b in range(a))
        print(f"q={q} m={len(S)} omega={q-1} affrank={ar} (m-1={len(S)-1}) pair l1 distances={sorted(dist)} W/logP={2*(q-1)} -> margin (d-W/2)/logP={min(dist)-(q-1)}")

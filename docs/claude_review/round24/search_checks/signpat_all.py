# half-plane sign patterns compatible with the argument principle for ALL ordered pairs (a,b):
# eps_A = +1 if rho_A in upper half plane; eps_{A^c} = -eps_A; need #{A: a in A, b notin A, eps_A=+1} = 4 for all a != b.
import itertools
FULL=frozenset(range(1,6))
TS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]  # A not containing 5
def sign(A,pat):  # pat dict over TS
    return pat[A] if 5 not in A else -pat[FULL-A]
good=[]; npair=0
for bits in itertools.product([1,-1],repeat=15):
    pat=dict(zip(TS,bits))
    ok=True
    for a in range(1,6):
        for b in range(1,6):
            if a==b: continue
            cnt=0
            for k in range(1,5):
                for A in itertools.combinations(range(1,6),k):
                    A=frozenset(A)
                    if a in A and b not in A and sign(A,pat)>0: cnt+=1
            if cnt!=4: ok=False; break
        if not ok: break
    if ok: good.append(bits)
print('patterns satisfying all 20 ordered-pair constraints:',len(good))
for g in good[:20]: print(''.join('+' if x>0 else '-' for x in g))

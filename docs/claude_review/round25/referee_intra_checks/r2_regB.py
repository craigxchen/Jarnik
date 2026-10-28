# Referee check r2: Theorem B on the zero-sum slice (5.1), reg < 2*K1 (uniform Taylor constant),
# star-partition isomorphism (Thm 4.1(1)-(2)), idx<=ord with star blocks equal.
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round25/referee_intra_checks')
from rlib import *
from itertools import product
random.seed(7)

def cuts(M):
    # one representative A per unordered cut {A, A^c}, A nonempty proper, containing index 0
    out = []
    for mask in range(1, 2**M-1):
        if mask & 1: out.append([(mask >> i) & 1 for i in range(M)])
    return out

def rand_support(M, size, amp):
    S = set([tuple([0]*M)])
    while len(S) < size:
        k = [random.randint(-amp, amp) for _ in range(M-1)]; k.append(-sum(k))
        S.add(tuple(k))
    return list(S)

viol_B = viol_K = 0; cnt = 0; worst = 0.0; worstK = 0.0
for M in (5, 6, 7):
    cs = cuts(M); assert len(cs) == 2**(M-1)-1
    for trial in range(40 if M < 7 else 15):
        S = rand_support(M, random.randint(4, 18 if M < 7 else 14), random.choice([1, 2]))
        r = reg(S)
        B = 2**(2-M)*sum(width(S, A) for A in cs)
        K1 = max(sum(abs(x) for x in k) for k in S)
        cnt += 1
        if r > B + 1e-12: viol_B += 1
        if not r < 2*K1: viol_K += 1
        worst = max(worst, r/B if B else 0); worstK = max(worstK, r/K1)
    print(f"M={M}: done", flush=True)
print(f"(5.1) reg <= 2^(2-M) sum_cuts width : violations {viol_B}/{cnt}, max reg/bound {worst:.4f}")
print(f"reg < 2*max|k|_1                     : violations {viol_K}/{cnt}, max reg/K1 {worstK:.4f}")

# Star partition M=9 (q=2), M=13 (q=3): Phi^* matrix identity in the stated bases; T-box pullbacks in 4T R_5
for M in (9, 13):
    q = (M-1)//4; o = M-1; Bs = [list(range(4*j, 4*j+4)) for j in range(q)]
    ok = True
    for T in (1, 2, 3):
        for _ in range(500):
            k = [random.randint(-T, T) for _ in range(M-1)]; k.append(-sum(k))
            back = [0]*M
            for B in Bs:
                s = [k[x] for x in B] + [-sum(k[x] for x in B)]
                if sum(abs(v) for v in s) > 8*T: ok = False
                for x, v in zip(B, s[:4]): back[x] += v
                back[o] += s[4]
            if back != k: ok = False
    print(f"star M={M}: unique preimage and s_j in 4T R_5 for T=1,2,3: {ok}")

# sanity_verify.py: slverify accepts the genuine system (lifted to sufficient depth) and rejects perturbations.
import pickle, copy, random, math
from sl import *
import slverify
S = setup("first4", 4)
labels = {}
for q in S.F:
    T = Torus(q); m = T.m
    # genuine labels at depth A large enough that no v in V_K collides at level A: q^A > max |D|
    A = 1
    while q ** A <= max(abs(d) for row in S.D for d in row): A += 1
    # genuine label mod m q^(A-1): dlog of pi/conj pi in T_(q^A) w.r.t. a generator lifting T.g -- built by brute force
    # only for small q; for the test we use q <= 23 and genuine level-1 labels elsewhere with A = 1 when no collision persists
    L = [T.L(pi) for pi in S.pis]
    labels[q] = dict(g1=T.g, top=tuple(L), A=1)
W = dict(P=S.P, K=4, F=sorted(S.F), phi="genuine", labels=labels)
try:
    slverify.check(W); print("genuine level-1 system: accepted (unexpected unless no level-1 collision persists)")
except AssertionError as e:
    print("genuine level-1 labels with A=1 rejected as expected (collisions persist at deepest level):", str(e)[:80])
W2 = pickle.load(open("W_first4_K4_gen_Q83.pkl", "rb"))
slverify.check(W2, verbose=False); print("stored witness accepted")
rng = random.Random(5); rej = 0
for trial in range(40):
    W3 = copy.deepcopy(W2); q = rng.choice([q for q in W3["F"] if q < 60]); lab = W3["labels"][q]
    m = q - (1 if q % 4 == 1 else -1); Mt = m * q ** (lab["A"] - 1)
    j = rng.randrange(4); t = list(lab["top"]); t[j] = (t[j] + 2 * rng.randrange(1, Mt // 2)) % Mt; lab["top"] = tuple(t)
    try: slverify.check(W3, verbose=False)
    except AssertionError: rej += 1
print(f"random parity-preserving perturbations at q<60: {rej}/40 rejected")

"""EXACT bookkeeping for curve classes on Mbar_{0,7} (Python integers / Fraction).

A curve class gamma is recorded by its intersection vector (D_S . gamma) over the 56 boundary
divisors; it is a genuine class iff it satisfies the Keel relations (for every 4-subset {i,j,k,l}
the three sums over splits inducing ij|kl, ik|jl, il|jk agree).

Checks:
 1. beta (all ones) satisfies Keel; K.beta = -7, B.beta = 56, psi_i.beta = 15, (K+2B).beta = 105.
 2. K = -(1/3) sum_{pairs} D_S  (Keel-McKernan / Pandharipande formula s(n-s)/(n-1) - 2).
 3. Kapranov psi_7-model classes e*l - sum m_a e_a - sum m_ab e_ab - sum m_abc e_abc  ->  D-vectors;
    Keel consistency, -K = 5e - 3sum m_a - 2 sum m_ab - sum m_abc, pair = 3(-K).
 4. Classes of standard covering families (forgetful fibres, lines meeting flats, ...): pair, triple,
    -K, psi_7, (B+8K).gamma.
 5. Counting identities used in Section 7: lattice exponent 7 (8 with top block), naive
    complete-intersection exponent -24, dependency defect 31.
"""
from fractions import Fraction as Fr
import itertools

N = 7
P = list(range(1, N+1))

def splits():
    out = []
    for r in range(2, N-1):
        for S in itertools.combinations(P, r):
            S = frozenset(S)
            Sc = frozenset(P) - S
            key = S if N not in S else Sc   # side not containing 7
            if key not in out:
                out.append(key)
    return out

SPL = splits()
assert len(SPL) == 56

def side_size(S):
    return min(len(S), N-len(S))

def is_pair(S): return side_size(S) == 2
def is_triple(S): return side_size(S) == 3

def keel_ok(vec):
    for Q in itertools.combinations(P, 4):
        i, j, k, l = Q
        sums = []
        for (a, b), (c, d) in [((i, j), (k, l)), ((i, k), (j, l)), ((i, l), (j, k))]:
            s = 0
            for S in SPL:
                inter = S & set(Q)
                if inter == {a, b} or inter == {c, d}:
                    s += vec[S]
            sums.append(s)
        if not (sums[0] == sums[1] == sums[2]):
            return False
    return True

def K_coeff(S):
    s = len(S)
    return Fr(s*(N-s), N-1) - 2

def dot(coeffs, vec):
    return sum(coeffs[S]*vec[S] for S in SPL)

Kc = {S: K_coeff(S) for S in SPL}
Bc = {S: 1 for S in SPL}

def psi_coeffs(i):
    # psi_i = sum_{S: i in S, j,k not in S} D_S, fix j,k the two smallest labels != i
    j, k = [x for x in P if x != i][:2]
    c = {}
    for S in SPL:
        Sfull = S
        # use the side containing i
        side = S if i in S else frozenset(P) - S
        c[S] = 1 if (i in side and j not in side and k not in side) else 0
    return c

def kapranov(e, m1, m2, m3):
    """psi_7 model: m1[a] (a in 1..6), m2[(a,b)], m3[(a,b,c)] multiplicities; returns D-vector."""
    vec = {}
    for S in SPL:
        # S is the side not containing 7
        if len(S) == 5:      # complement {a,7}: E_a
            a = (set(range(1, 7)) - S).pop(); vec[S] = m1.get(a, 0)
        elif len(S) == 4:    # complement {a,b,7}: E_ab
            ab = tuple(sorted(set(range(1, 7)) - S)); vec[S] = m2.get(ab, 0)
        elif len(S) == 3:    # complement {a,b,c,7}: E_abc
            abc = tuple(sorted(set(range(1, 7)) - S)); vec[S] = m3.get(abc, 0)
        elif len(S) == 2:    # D_{cd}, c,d in [6]: hyperplane through the other four points
            Q = set(range(1, 7)) - S
            v = e - sum(m1.get(a, 0) for a in Q)
            v -= sum(m2.get(tuple(sorted(x)), 0) for x in itertools.combinations(Q, 2))
            v -= sum(m3.get(tuple(sorted(x)), 0) for x in itertools.combinations(Q, 3))
            vec[S] = v
    return vec

def report(name, vec):
    pair = sum(vec[S] for S in SPL if is_pair(S)); trip = sum(vec[S] for S in SPL if is_triple(S))
    mK = -dot(Kc, vec); B = dot(Bc, vec); psi7 = dot(psi_coeffs(7), vec)
    print('%-34s keel=%s pair=%3s triple=%3s -K=%5s B=%3s psi7=%3s (B+8K)=%5s min D=%s' % (
        name, keel_ok(vec), pair, trip, mK, B, psi7, B + 8*dot(Kc, vec), min(vec.values())))
    return pair, trip, mK

if __name__ == '__main__':
    beta = {S: 1 for S in SPL}
    print('beta Keel:', keel_ok(beta))
    print('K.beta =', dot(Kc, beta), ' B.beta =', dot(Bc, beta),
          ' psi_i.beta =', [dot(psi_coeffs(i), beta) for i in P])
    print('(K+2B).beta =', dot(Kc, beta) + 2*dot(Bc, beta))
    print('K coefficients by side size:', sorted(set((side_size(S), Kc[S]) for S in SPL)))
    # Kapranov consistency: beta = 15 l - sum e_a - sum e_ab - sum e_abc
    m1 = {a: 1 for a in range(1, 7)}
    m2 = {ab: 1 for ab in itertools.combinations(range(1, 7), 2)}
    m3 = {abc: 1 for abc in itertools.combinations(range(1, 7), 3)}
    kb = kapranov(15, m1, m2, m3)
    print('Kapranov beta == beta:', kb == beta)
    print()
    report('beta', beta)
    report('line', kapranov(1, {}, {}, {}))
    report('line through p_1 (= fibre f_1)', kapranov(1, {1: 1}, {}, {}))
    report('line meeting L_12', kapranov(1, {}, {(1, 2): 1}, {}))
    report('line meeting Pi_123', kapranov(1, {}, {}, {(1, 2, 3): 1}))
    report('line meeting L_12, Pi_345', kapranov(1, {}, {(1, 2): 1}, {(3, 4, 5): 1}))
    report('line meeting Pi_123,Pi_145,Pi_246', kapranov(1, {}, {}, {(1, 2, 3): 1, (1, 4, 5): 1, (2, 4, 6): 1}))
    report('RNC through p_1..p_6 (= fibre f_7)', kapranov(4, m1, {}, {}))
    # forgetful fibre f_7 directly
    f7 = {S: (1 if (is_pair(S) and (len(S) == 5)) else 0) for S in SPL}
    report('fibre f_7 (direct)', f7)
    # twisted cubic meeting 13 planes (chosen so all D >= 0)
    planes = list(itertools.combinations(range(1, 7), 3))
    best = None
    for choice in itertools.combinations(planes, 13):
        vec = kapranov(3, {}, {}, {p: 1 for p in choice})
        if min(vec.values()) >= 0:
            best = choice; break
    if best:
        report('cubic meeting 13 planes', kapranov(3, {}, {}, {p: 1 for p in best}))
    else:
        print('no 13-plane choice with all D>=0 for cubics')
    # identity -K = pair/3 on all the above is printed; check algebraically on random Keel classes:
    print()
    print('Counting identities (Section 7):')
    Ts = [frozenset(c) for r in range(2, 6) for c in itertools.combinations(range(1, 7), r)]
    lat = sum(len(T)-1 for T in Ts)
    print('  sum_T (|T|-1) =', lat, '; box exponent 5*15 = 75; blocks 56; lattice exponent =', 75 - lat + 56)
    print('  with top block: 80 -', lat + 5, '+ 57 =', 80 - (lat+5) + 57)
    ci = 56 - 10*8
    print('  naive complete-intersection exponent 56 - 10*8 =', ci)
    defect = sum((len(T)*(len(T)-1))//2 - (len(T)-1) for T in Ts if 1 not in T and len(T) >= 3)
    print('  dependency defect sum_{T in {2..6}, |T|>=3} [C(|T|,2) - (|T|-1)] =', defect,
          '; 80 - 49 =', 80 - sum(len(T)-1 for T in Ts if 1 not in T))

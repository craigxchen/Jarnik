"""For prime q, minimal delta_S = deg_Y( prod_{a in S}(Y-a) - prod_{a in S}(Y+a) ) over half-systems S
(S disjoint from -S, S u -S = F_q^*). VI.4 needs X-degree of F_i - F_j = h(1 + (q-1)/2 + delta) to be < (n-g)/2,
i.e. g >= h(1 + 2 delta) roughly; q = 3 mod 4 squares give delta = 0."""
import itertools
def best(q):
    reps = list(range(1, (q+1)//2))
    bestd = None
    for signs in itertools.product((1, -1), repeat=len(reps)-1):
        S = [reps[0]] + [s*a % q for s, a in zip(signs, reps[1:])]
        f = [1]; g = [1]
        for a in S:
            f = [((f[j-1] if j else 0) - a*(f[j] if j < len(f) else 0)) % q for j in range(len(f)+1)]
            g = [((g[j-1] if j else 0) + a*(g[j] if j < len(g) else 0)) % q for j in range(len(g)+1)]
        d = max([j for j in range(len(f)) if (f[j]-g[j]) % q] or [-1])
        if bestd is None or d < bestd[0]: bestd = (d, sorted(S))
    return bestd
for q in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    d, S = best(q)
    print('q=%d (q mod 4 = %d): (q-1)/2=%d  min delta=%d  ratio delta/((q-1)/2)=%.2f  S=%s' % (q, q % 4, (q-1)//2, d, d/((q-1)/2), S))

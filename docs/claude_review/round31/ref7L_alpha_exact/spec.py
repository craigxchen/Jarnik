"""Integer specialisation of the five marks, exact Gaussian integer arithmetic (python ints)."""
import sys, json, math
from math import isqrt

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def rdiv(x, m):  # nearest integer to x/m, m>0
    return (2*x + m)//(2*m)
def gdivmod(a, b):
    nb = gnorm(b); num = gmul(a, gconj(b))
    q = (rdiv(num[0], nb), rdiv(num[1], nb))
    r = gsub(a, gmul(q, b))
    assert gnorm(r) < nb
    return q, r
def gexact(a, b):
    q, r = gdivmod(a, b); assert r == (0, 0), (a, b); return q
def gdivides(b, a):
    return gdivmod(a, b)[1] == (0, 0)
def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gdivmod(a, b)[1]
    return a
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
def normalize(g):  # choose associate with re>0, im>=0
    for u in UNITS:
        h = gmul(g, u)
        if h[0] > 0 and h[1] >= 0: return h
    return g
def geval(coeffs, T):
    acc = (0, 0)
    for c in reversed(coeffs):
        acc = (acc[0]*T + c[0], acc[1]*T + c[1])
    return acc
def unit_primary(z):
    """unit epsilon with z == epsilon mod (2+2i) (primary-prime convention); z odd."""
    for u in UNITS:
        d = gsub(z, u)
        if gdivides((2, 2), d): return u
    raise ValueError('even')

def analyse(data, T):
    zs = [geval(c, T) for c in data['coeffs']]
    g = zs[0]
    for z in zs[1:]: g = ggcd(g, z)
    g = normalize(g)
    N = data['bezout_N']
    g_divides_N = gdivides(g, (N, 0))
    q = [gexact(z, g) for z in zs]
    R2 = gnorm(q[0])
    same_norm = all(gnorm(x) == R2 for x in q)
    gg = q[0]
    for x in q[1:]: gg = ggcd(gg, x)
    prim = gnorm(gg) == 1
    distinct = len(set(q)) == 5
    assoc = False; cassoc = False
    for a in range(5):
        for b in range(a+1, 5):
            for u in UNITS:
                if gmul(u, q[b]) == q[a]: assoc = True
                if gmul(u, gconj(q[b])) == q[a]: cassoc = True
    d2 = [gnorm(gsub(q[a], q[b])) for a in range(5) for b in range(a+1, 5)]
    logR = 0.5*math.log(R2)
    e = 0.5*math.log(max(d2))/logR
    emin = 0.5*math.log(min(d2))/logR if min(d2)>0 else float('-inf')
    gamma = 2/(1-e)
    units = [unit_primary(x) for x in q] if R2 % 2 == 1 else None
    # convention-free pair test: zeta_a/zeta_b = +-P/conj(P), P conj-primitive  <=> Re(zeta_a conj zeta_b) odd
    pm = [gmul(q[a], gconj(q[b]))[0] % 2 == 1 for a in range(5) for b in range(a+1, 5)]
    # t_ab := |c_ab|/2, c_ab=(za-zb)/gcd(za,zb): is |c|^2/4 a perfect square integer?
    tint = []
    for a in range(5):
        for b in range(a+1, 5):
            c = gexact(gsub(q[a], q[b]), ggcd(q[a], q[b]))
            nc = gnorm(c)
            tint.append(nc % 4 == 0 and isqrt(nc//4)**2 == nc//4)
    # constant: maxdist / (R''^{1/2} T^{1/2})
    return dict(T=T, gcd=g, gcd_norm=gnorm(g), g_divides_N=g_divides_N, R2_digits=len(str(R2)),
                same_norm=same_norm, primitive=prim, distinct=distinct, associates=assoc,
                conj_associates=cassoc, e=e, emin=emin, gamma=gamma,
                units=units, pm_pairs_all=all(pm), pm_count=sum(pm), t_integer_count=sum(tint),
                const=math.exp(0.5*math.log(max(d2))-0.5*logR-0.5*math.log(abs(T))))

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    Ts = [int(x) for x in sys.argv[2].split(',')]
    n = data['n']; D = data['D']
    print('n=%d seeds=%s D=%d target e=%.6f target gamma=%.6f N=%d' % (n, data['seeds'], D, 0.5+1/(2*D), 4+1/(2*n), data['bezout_N']))
    allg = set()
    for T in Ts:
        r = analyse(data, T)
        allg.add(r['gcd_norm'])
        print('T=%-12d Ng=%-10d g|N=%s normEq=%s prim=%s dist=%s assoc=%s cassoc=%s e=%.6f emin=%.6f gamma=%.6f c=%.3f units=%s pm=%d/10 tint=%d/10' % (
            T, r['gcd_norm'], r['g_divides_N'], r['same_norm'], r['primitive'], r['distinct'], r['associates'], r['conj_associates'],
            r['e'], r['emin'], r['gamma'], r['const'], ''.join({(1,0):'1',(0,1):'i',(-1,0):'-',(0,-1):'j'}[u] for u in r['units']) if r['units'] else None, r['pm_count'], r['t_integer_count']))
    print('gcd norms seen:', sorted(allg))

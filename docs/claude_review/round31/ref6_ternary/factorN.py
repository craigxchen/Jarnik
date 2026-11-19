"""Structure of N for the four-mark example (factorisation of N in F_3[t] by python-flint, used only as a factoring tool)."""
import pickle, collections
from flint import nmod_poly
from f9lib import F9P
d = pickle.load(open('four_mark_L11.pkl', 'rb'))
N = F9P(*d['N']); assert N.is_f3()
Np = nmod_poly([int(x) for x in N.re], 3)
c, fac = Np.factor()
print('deg N', Np.degree(), 'lead', c, '#distinct primes', len(fac))
hist = collections.Counter((f.degree(), e) for f, e in fac)
print('(prime degree, multiplicity): count')
for key in sorted(hist): print('  ', key, hist[key])
odd = sum(f.degree()*e for f, e in fac if f.degree() % 2 == 1)
print('total degree in odd-degree (inert) primes:', odd)
sqf = sum(f.degree() for f, e in fac); print('degree of radical', sqf)
# gcd extra part
G = F9P(*d['G'])
print('G (full gcd) degree', G.deg())

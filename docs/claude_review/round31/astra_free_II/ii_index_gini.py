"""Exact extremal Gini bound for the small odd primes in the balanced-index lemma."""
from fractions import Fraction
from itertools import combinations_with_replacement, product

def gini(u, v):
    return sum(abs(i-j)*a*b for i,a in u.items() for j,b in v.items())

def max_gini(c):
    # P_r is a mixture of the fair law and the fair law conditional on a nonconstant pattern.
    fair, nc = {}, {}
    for bits in product((0,1), repeat=5):
        x = sum(b*k for b,k in zip(bits,c))
        fair[x] = fair.get(x,0) + Fraction(1,32)
        if sum(bits) not in (0,5):
            nc[x] = nc.get(x,0) + Fraction(1,30)
    a,b,e = gini(fair,fair),gini(fair,nc),gini(nc,nc)
    lin, quad = 2*(b-a), a-2*b+e
    vals = [a,e]
    if quad and 0 < -lin/(2*quad) < 1:
        vals.append(a-lin*lin/(4*quad))
    return max(vals)

def centered(x,q):
    return (x+q//2)%q-q//2

def orbit(c,q):
    return {tuple(sorted(centered(k*x,q) for x in c)) for k in range(1,q)}

def check(q):
    seen=set()
    worst, examples, n=0,[],0
    centered_residues=range(-(q//2),q//2+1)
    for c in combinations_with_replacement(centered_residues,5):
        if not any(c) or sum(c)%q or c in seen:
            continue
        orb=orbit(c,q)
        seen.update(orb)
        n+=1
        floor=min(max_gini(kc) for kc in orb)
        if floor>worst:
            worst,examples=floor,[c]
        elif floor==worst:
            examples.append(c)
    assert worst<Fraction(q,2)
    print(q, n, worst, examples)

if __name__=='__main__':
    for q in (3,5,7):
        check(q)

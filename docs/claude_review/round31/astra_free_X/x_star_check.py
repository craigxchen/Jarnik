"""Exact modular and polynomial-scale negative controls for the all-height ray construction."""
from itertools import product
from math import gcd, isqrt, lcm

def primes_to(limit):
    return [n for n in range(2,limit+1) if all(n%j for j in range(2,isqrt(n)+1))]

def allowed(value, modulus, r, parity, hub):
    if r == 2:
        if value % 2 != parity:
            return False
        if modulus < 16:
            return True
        criterion = 16
    else:
        criterion = r*r if r in (3,5) else r
        if modulus < criterion:
            return True
    return all((value-shift)%criterion for shift in (0,hub,-hub))

def verify_sumsets():
    checks=0
    for r in primes_to(43):
        modulus=2
        if r==2:
            powers=(2,4,8,16,32)
        else:
            powers=(r,r*r) if r<7 else (r,)
        for modulus in powers:
            for hub in range(modulus):
                if gcd(hub,modulus) not in (1,2):
                    continue
                parity_classes=(0,1) if r==2 else (0,)
                lists=[[x for x in range(modulus) if allowed(x,modulus,r,parity,hub)] for parity in parity_classes]
                for i,j in product(parity_classes,repeat=2):
                    realized={(x+y)%modulus for x in lists[i] for y in lists[j]}
                    expected=set(range((i+j)%2,modulus,2)) if r==2 else set(range(modulus))
                    assert realized==expected
                    checks+=len(realized)
    print('attainable target sums, over all local component choices:',checks)

def ray_levels():
    # Enumerate orders only. The exhaustive alpha-loop checks all exponent
    # classes accepted by gcd | 120, with an explicit worst word witnessing
    # its first fourth-root-of-unity collision.
    checks=0
    largest_ratio=(0,1)
    max_level=3000
    for ell in primes_to(max_level):
        if ell==2:
            continue
        power=ell
        while power<=max_level:
            order=(ell-(1 if ell%4==1 else -1))*(power//ell)
            for x in range(order):
                g=gcd(order,x)
                if 120%g:
                    continue
                first=order//gcd(order,4*x)
                assert order != 0 and first >= 1 and (4*first*x)%order==0
                assert 4*first*120 % order == 0
                assert power<=720*first
                if power*largest_ratio[1]>first*largest_ratio[0]:
                    largest_ratio=(power,first)
                checks+=1
            power*=ell
    print('admissible exponent classes at all odd prime powers <=3000:',checks)
    print('largest tested collision level / first ray multiplier:',largest_ratio)

if __name__=='__main__':
    verify_sumsets()
    ray_levels()

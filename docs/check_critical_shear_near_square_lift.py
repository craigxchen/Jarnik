"""Exact integer denominators and circle norms for the critical shear lift."""
from itertools import product
from math import gcd, lcm

from check_gaussian_reflection_replacement import gmul, gconj, gnorm, gcd_many
from check_mobius_conductor_transfer import fixtures, least_norm
from check_mobius_joint_anchor_residues import factors, gaussian_prime, gv


def main():
    sources = fixtures()
    sources.append([(1,0),(1,1),(1,2),(3,1),(3,2),(5,1),(5,2),(2,5)])
    count = nontrivial = two_tail = 0
    for A,D,B in product(range(1,5),range(-4,5),range(-5,6)):
        if not D or (B == 0 and A == abs(D)):
            continue
        U,V = (A+D,-B),(A-D,B)
        Gamma = (B*B-A*A+D*D,2*A*B)
        K,T = gnorm(Gamma),A*A+B*B+D*D
        assert K == gnorm(U)*gnorm(V) == T*T-4*A*A*D*D > 0
        for rows in sources:
            N = least_norm(rows)
            ps = set()
            for H in rows:
                ps.update(factors(gnorm(H)))
            Nclip = 1
            expected_h = 1
            both = False
            for p in ps:
                if p % 4 != 1:
                    continue
                pi = gaussian_prime(p)
                barpi = gconj(pi)
                d = gv(V,pi)-gv(U,pi)
                e = gv(V,barpi)-gv(U,barpi)
                left,right = sorted((d,-e))
                ts = [gv(H,pi)-gv(H,barpi) for H in rows]
                clips = [max(left,min(right,t)) for t in ts]
                Nclip *= p**(max(clips)-min(clips))
                g,gbar = gv(Gamma,pi),gv(Gamma,barpi)
                assert -g <= left <= right <= gbar
                lt,rt = max(0,-g-min(ts)),max(0,max(ts)-gbar)
                expected_h *= p**max(lt,rt)
                both |= lt > 0 and rt > 0
            assert K % Nclip == 0 and N % Nclip == 0
            erased, cofactor = N//Nclip, K//Nclip
            Kplus, Kminus = B*B+(A+D)**2, B*B+(A-D)**2
            nplus = gcd(Nclip,Kplus)
            nminus = Nclip//nplus
            assert Kminus % nminus == 0
            uplus, uminus = Kplus//nplus, Kminus//nminus
            assert uplus*uminus == cofactor
            assert uplus*nplus-uminus*nminus == 4*A*D
            assert (4*A*D) % gcd(nplus,nminus) == 0
            h = 1
            numerators = [gmul(Gamma,gmul(H,H)) for H in rows]
            for H,z in zip(rows,numerators):
                h = lcm(h,gnorm(H)//gcd(gnorm(H),*z))
            assert h == expected_h and h % 2 == 1
            assert (N//Nclip) % h == 0
            lifted = []
            for H,z in zip(rows,numerators):
                f = gnorm(H)
                assert h*z[0] % f == 0 and h*z[1] % f == 0
                lifted.append((h*z[0]//f,h*z[1]//f))
            X,Y = h*T,2*h*A*D
            assert h*h*K == X*X-Y*Y
            assert all(gnorm(z) == h*h*K for z in lifted)
            assert lifted[0] == (X-2*h*A*A,2*h*A*B)
            assert h*h*K % N == 0
            common_norm = gnorm(gcd_many(lifted))
            assert common_norm*N == h*h*K
            assert common_norm*erased == h*h*cofactor
            ordinary_content = gcd(*(x for z in lifted for x in z))
            assert gcd(ordinary_content,h) == 1
            assert all(x % ordinary_content == 0 for x in Gamma)
            if gcd(A,B,D) == 1:
                integer_center_content = gcd(ordinary_content,h*T)
                assert 2*A % integer_center_content == 0
            count += 1
            nontrivial += h > 1
            two_tail += both
    print(f'PASS: {count} triangular matrix/source configurations, {nontrivial} nontrivial integer denominators, {two_tail} with missing factors in both orientations.')
    print('PASS: coefficient interval containment, exact minimal integer denominator, and h divides N/Nclip.')
    print('PASS: all lifted integer points lie on one near-square circle, including mixed source units and exact common content.')
    print('PASS: complementary small-cofactor source divisors, minimal-denominator coprimality, and integer-center content.')


if __name__ == '__main__':
    main()

# Exact radius cost of a critical half-angle power map

The rational circle map induced by u → u^d, for odd d, contracts angles
to order d at its actual anchor. Its least integral target radius is a
Gaussian least common multiple, rather than a degree-d height of the source
circle. In particular:

* This canonical map never decreases the primitive radius for odd d.
* For d=3, an explicit primitive four-point endpoint family has
  R_new = R_old^(5+o(1)). Its angular width contracts cubically, but
  its normalized endpoint constant tends to infinity.

The family disproves a general R_new=O(R_old^3) cost bound and a general
bounded-endpoint-constant guarantee for this canonical cubic map. It does
not settle behavior restricted to arbitrarily many points, nor rule out
other nonlinear maps or prove a radius-independent point count.

The [existing fixed-degree inflation theorem](fixed_degree_rational_map_inflation.md)
already treats general fixed-degree rational maps on sufficiently many
extracted rows. The result here supplies an exact least-radius dictionary
and explicit four-point failures, below that theorem's row-count regime.

## Exact least-radius dictionary

Let z_0,...,z_(M-1) be a primitive Gaussian tuple of common modulus R.
Include the actual anchor, and write its rational relative directions as

    z_i/z_0 = (a_i+i b_i)/(a_i-i b_i),
    gcd_Z(a_i,b_i)=1,        (a_0,b_0)=(1,0).

The antipode, if present, uses (0,1). For odd positive d, define

    e_i = 1 if a_i,b_i are both odd, and 0 otherwise,
    H_i^(d) = (a_i^d+i b_i^d)/(1+i)^e_i,
    w_i^(d) = (a_i^d+i b_i^d)/(a_i^d-i b_i^d).

Each H_i^(d) is a Gaussian integer coprime to its conjugate. Indeed,
any common odd Gaussian prime would divide both coprime ordinary
integers a_i^d,b_i^d; and the numerator has exactly one factor 1+i
when both are odd, since its norm is two modulo four. Consequently

    w_i^(d) = i^e_i H_i^(d)/bar(H_i^(d)).

The exact least integral realization is

    L_d = lcm_(Z[i],i) bar(H_i^(d)),
    Z_i^(d)=L_d w_i^(d),
    R_d=|L_d|.                                             (1)

The lcm is defined up to a Gaussian unit. Every integral realization
has Gaussian anchor alpha because w_0=1, and coprimality forces
each bar(H_i^(d)) to divide alpha. Thus (1) is minimal, and its
output tuple is Gaussian primitive. It includes every source denominator
and all common-factor savings. Equivalently, its squared radius is

    R_d^2=product_(Gaussian primes rho)
           Norm(rho)^max_i v_rho(bar H_i^(d)).              (2)

No estimate replacing a maximum over rows by a single degree is implicit.
For d=1, source primitivity gives |L_1|=R.

## No radius descent for any odd power

Put sigma=(-1)^((d-1)/2). Polynomial division gives

    a^d+i b^d is divisible in Z[i] by a+sigma i b.

After removing the common single ramified factor, H_i^(d) is divisible
by H_i^(1) when sigma=1, and by a Gaussian associate of
bar(H_i^(1)) when sigma=-1. The same orientation applies to every
row. Taking Gaussian least common multiples proves the exact bound

                                  R_d >= R.               (3)

For comparison, the elementary upper bound from the product of the row
denominators is only

    R_d <= 2^[(d-1)(M-1)/2] R^[d(M-1)].                   (4)

Here Norm H_i^(1) divides R^2 by the source pair dictionary, and
Norm H_i^(d) <= 2^(d-1) Norm(H_i^(1))^d. Formula (4) is deliberately
crude; it displays the row count which a degree-only bound would omit.

## A primitive four-point endpoint source

Let n be a positive multiple of ten, with n>=20. Take the four
actual Gaussian integers

    z_(sigma,tau)=(n+sigma i)(n+1+tau i)/(1+i),
    sigma,tau in {+1,-1},
    z_0=z_(-1,-1).

They have common radius

    R=[(n^2+1)((n+1)^2+1)/2]^(1/2) ~ n^2/sqrt(2).        (5)

They are Gaussian primitive: the gcd of n+i,n-i is a unit, and the
gcd of n+1+i,n+1-i is an associate of 1+i. Taking all independent
products takes the sum of the two minimum valuations at every prime;
division by 1+i removes the complete common factor.

Their four anchor-relative half-angle parameters are

    0,       1/n,       1/(n+1),       B/A,
    A=n^2+n-1,          B=2n+1.                          (6)

The last fraction is primitive: gcd(A,B) divides five from
4A=B^2-5, while n=0 mod 5 gives B=1 mod 5.
The angular width is

    Delta=2 arctan(1/n)+2 arctan(1/(n+1)) ~ 4/n.

In particular Delta sqrt(R)<4 for this family. One elementary bound
uses Delta<4/n, sqrt(R)<(n+2)/2^(1/4), and
1+2/n<=11/10<2^(1/4) for n>=20.

## Exact cubic denominators and their bounded gcd costs

For the three nonanchor rows in (6), the reduced cubic Gaussian factors
are exactly

    C_1=n^3+i,
    C_2=((n+1)^3+i)/(1+i),
    C_3=(A^3+i B^3)/(1+i).

Thus R_3=|lcm(C_1,C_2,C_3)|; conjugating the entire lcm in (1)
does not change its norm. An alternative exact expression retaining all
gcds is

    R_3 = |C_1 C_2 C_3| |gcd(C_1,C_2,C_3)|
          /[|gcd(C_1,C_2)| |gcd(C_1,C_3)| |gcd(C_2,C_3)|].   (7)

Factor the unreduced cubic numerators as

    D=(n-i)(n+1-i),
    Q_1=n^2+i n-1,
    Q_2=(n+1)^2+i(n+1)-1,
    Q_3=A^2+i AB-B^2,

    C_1=(n-i)Q_1,
    C_2=(n+1-i)Q_2/(1+i),
    C_3=D Q_3/(1+i).                                    (8)

The six pairwise polynomial resultants over Z[i], in the order
(D,Q_1),(D,Q_2),(D,Q_3),(Q_1,Q_2),(Q_1,Q_3),(Q_2,Q_3), are

    6+9i,       6-9i,       225,       -2,       9i,       -9i.

They are all nonzero. The polynomial Bezout identities imply that each
pairwise gcd after any integral substitution n divides its fixed
resultant. Their absolute values have product K=4,264,650.

Let P=D Q_1 Q_2 Q_3. The lcm of these four factors has absolute value
at least |P|/K: prime by prime, the product divided by the lcm divides
the product of the six pairwise gcds. Every factor divides (1+i) times
the cubic lcm, by (8), while every C_i divides P. Therefore all
remaining cancellation costs are bounded explicitly:

    |P|/(sqrt(2) K) <= R_3 <= |P|.                       (9)

The polynomials have respective degrees 2,2,2,4 and leading coefficient
one. Consequently

    R_3 asymp n^10 = R^(5+o(1)).                         (10)

Finally the transformed angular width is exactly

    Delta_3=2 arctan((B/A)^3) ~ 16/n^3.

Equations (9)--(10) give Delta_3 sqrt(R_3) asymp n^2, which tends
to infinity. The critical angular contraction is real, but clearing the
three simultaneous Gaussian denominators destroys the endpoint scale
even on this literal four-point source family.

## A cubic map preserving more old factors still loses the endpoint scale

Consider the improved critical map

    f(u)=2u^3/(1+3u^2),
    K(a,b)=(a+i b)^2(a-2i b)
          =a^3+3ab^2+2i b^3.                            (11)

It has angular vanishing order three and retains a square of the old
Gaussian factor. Restrict the same actual source family to n=20k.
For each of its three nonanchor primitive pairs, b is odd and a is
either odd or divisible by four. The real and imaginary coordinates
in (11) have ordinary gcd exactly two; after division by two they
are coprime and have opposite parity. Thus the three exact reduced
Gaussian numerators are

    T_1=(n+i)^2(n-2i)/2,
    T_2=(n+1+i)^2(n+1-2i)/2,
    T_3=(A+iB)^2(A-2iB)/2.                              (12)

Each is coprime to its conjugate. Hence the exact least target radius
is |lcm(T_1,T_2,T_3)|, with the same complete gcd formula as (7).

Write

    E=(n+i)(n+1+i),       G_1=n-2i,       G_2=n+1-2i,
    G_3=A-2iB,           P_*=E^2 G_1 G_2 G_3.

The six polynomial resultants for E^2,G_1,G_2,G_3, in pair order, are

    72-54i,       72+54i,       2025,       1,       3,       3.

Their absolute values have product K_*=147,622,500. Every one of these
four factors divides twice the lcm in (12), while every T_i divides
P_*. The same exact valuation argument as before gives

    |P_*|/(2K_*) <= R_* <= |P_*|.                       (13)

The four factors have degrees 4,1,1,2, respectively, so

    R_* asymp n^8 = R^(4+o(1)).

The function f is strictly increasing for positive u. Consequently the
output angular width is exactly 2 arctan(f(B/A)), asymptotic to 32/n^3.
Its normalized endpoint constant is therefore asymptotic up to fixed
positive factors to n, and again tends to infinity. Preserving the old
factor twice improves this particular cost from exponent five to four;
it still exceeds the exponent three tolerated by cubic flattening.

[The exact checker](check_nonlinear_half_angle_power_radius_cost.py)
verifies Gaussian primitivity, the half-angle identities, the least-lcm
tuple, the complete three-factor gcd formula, the six polynomial
resultants for both cubic maps, and the explicit cancellation bounds. It does not infer
an unbounded-point endpoint family from these four-point examples.

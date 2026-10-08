# Joint axis lift: triangle gcd and Pell-coefficient audit

Consider three primitive Gaussian pair factors A_1=A_ab, A_2=A_bc, and
A_3=A_ac with

    A_1 A_2 = q A_3,       q an odd rational integer
    d_i = Norm(A_i) | N.

Put g_i=N/d_i and Z_i=g_i A_i^2. The norm identity gives
d_1 d_2=q^2 d_3. Consequently the lifted points satisfy the exact relation

    Z_1 Z_2 = N Z_3.                                  (1)

All three have Gaussian norm N^2. This is stronger than treating the three
axis points independently, but its primewise content has to be tracked before
extracting a common factor.

## Primewise common-gcd formula

Fix a Gaussian prime pi, write n=v_pi(N), and set

    e_i = v_pi(Z_i),       u_i = n-e_i.

From (1),

    u_3 = u_1+u_2.                                    (2)

The conjugate valuation is v_barpi(Z_i)=2n-e_i=n+u_i, since
Norm(Z_i)=N^2. Therefore the common Gaussian gcd
D=gcd(Z_1,Z_2,Z_3) has exponents

    v_pi(D)      = n - max(u_1,u_2,u_3)
    v_barpi(D)   = n + min(u_1,u_2,u_3).               (3)

whenever these displayed exponents are nonnegative, as they must be for an
actual common divisor. After dividing by D, the residual valuation range at
the pair pi,barpi is exactly

    max(u_i)-min(u_i).                                (4)

Thus common-gcd extraction removes the common baseline but leaves the
primewise spread. Equation (2) alone does not bound that spread.

In the primitive fixed-unit source setting, N is odd and has only split
rational prime factors, so inert primes and the prime two do not occur in
this calculation.

## Relation with the Gaussian factors

If alpha_i=v_pi(A_i) and baralpha_i=v_barpi(A_i), then

    u_i = baralpha_i-alpha_i.

Primitivity gives alpha_i*baralpha_i=0 at a split prime, while
A_1 A_2=q A_3 gives the additive relation (2), with the rational valuation
of q cancelling from the centered quantities. Hence the joint lift does
provide a clean valuation tree: the two input centered valuations add to the
output centered valuation. It does not, by itself, force the centered
valuations or their spread to be small.

The residual points after dividing by D have the rational Gaussian relation

    (Z_1/D)(Z_2/D) = (N/D)(Z_3/D),

where N/D need not be an algebraic integer: (3) allows the conjugate
exponent of D to exceed v_barpi(N). Even in cases where D divides N, the
coefficient N/D can retain prime powers allowed by (3). Calling this
coefficient squarefree, integral in general, or bounded by the small
imaginary coordinates would require an additional argument; none follows
from the triangle identity or from the individual bounds h_i=2g_i b_i^2.

## Consequence for a joint counting route

The exact identity (1) is a useful new reduction for triples and cycles.
However, common-gcd normalization does not yield a uniform Pell coefficient
bound under the current hypotheses. Any successful joint route must control
the primewise spreads in (4), likely using simultaneous restrictions at
several rational primes or a new relation involving the odd parts of the
b_i. The large individual gcds g_i do not cancel those spreads.

Accordingly this triangle calculation supplies an exact structural lemma but
no new uniform endpoint count on its own.

## Rational content of an entire pair family

There is a separate, exact statement for a source subset I. Let every pair
lift Z_ij (i<j in I) have integer-coordinate gcd

    c_ij = N/d_ij,

as holds for odd N and primitive pair factors. Define the rational content of
the whole family by

    D_I = gcd over i<j in I of all integer coordinates of Z_ij.

Since the gcd of the two coordinates of each individual Z_ij is c_ij,
taking the gcd over the union of all coordinates gives

    D_I = gcd_{i<j} c_ij
        = gcd_{i<j} (N/d_ij)
        = N / lcm_{i<j}(d_ij).                         (5)

The last equality is the elementary divisor identity valid because every
d_ij divides N. Thus D_I is an ordinary integer and should not be confused
with the Gaussian gcd D of the three complex numbers in the triangle
calculation above.

Take G to be the Gaussian gcd of the source points indexed by I, and write
z_i=G w_i. At a split prime with source allocation exponents e_i in [0,n],
the exponent in lcm(d_ij) is max(e_i)-min(e_i); the exponent in Norm(G) is
n-max(e_i)+min(e_i). Hence

    N / Norm(G) = lcm_{i<j}(d_ij).

Thus (5) gives D_I=Norm(G), and directly

    Z_ij / D_I = w_i * conjugate(w_j).                 (6)

So ordinary joint normalization commutes exactly with source extraction. It
does not create a smaller radius or a new Pell coefficient: after dividing
by D_I one has simply recovered the inherited pair products.

Finally, the Gaussian gcd of the positive-oriented lifts can be strictly
larger than their rational content and can be nonreal. If `D` is a Gaussian
common divisor, then dividing by `D` changes the radius from `N` to

    N / sqrt(Norm(D)),

and the squared radius is `N^2/Norm(D)`. This need not be an integer-radius
circle even when the squared radius is an integer. Therefore a Chan/Turk-
style integer-radius simultaneous-Pell application cannot be inferred from
Gaussian common divisibility alone; it requires a separate proof that
`Norm(D)` is a square (and, for a rational joint coefficient, that `D` is
rational up to a unit).

## Exact four-point orientation experiment

The source blocks

    (-1,2), (13,8), (8,5), (50,31)

with rows `1111`, `0100`, `0010`, `0001` give
`N=358853785`.  For the three edges on source rows `(0,1,2)`, orienting
all lifted products to have positive imaginary part gives ordinary content
`1` and Gaussian gcd `(-144,-1)`, of norm `20737`.  The normalized squared
radius is exactly

    N^2/20737 = 6209964749425,

which is nonsquare, so the normalized radius is not an integer.

There is a useful warning in the other direction.  For the triangle on
source rows `(1,2,3)`, take the mixed orientation
`(Z_12, conjugate(Z_13), Z_23)`.  Its Gaussian gcd is the rational associate
`(-5,0)`, so division by `5` leaves the integer radius
`R=N/5=71770757`.  The three normalized points have coordinates

    (71763835,-996768), (71769155,479532), (71768893,517260),

and satisfy the exact individual Pell equations
`Y^2=h(2R-h)` with `h=6922,1602,1864`.  However, none of the three pair
products equals `+R` or `-R` times the remaining normalized point.  Thus
this orientation preserves an integer radius but changes the form of the
joint relation. It still satisfies the exact cyclic identity
`W_1 W_2 W_3=R^3`, where W denotes the normalized points. The checker
exhausts all eight conjugation choices for this triangle: every orientation
retaining a two-factor relation with coefficient `+N` has nonsquare gcd
norm, while the square-norm choices have no such two-factor relation
without conjugation. This finite example does not exclude other useful
orientations or applications of the retained cyclic identity.

The persistent [checker](check_joint_axis_pell_coefficients.py) verifies
ordinary normalization on all eleven source subsets of size at least two,
including its exact agreement with source Gaussian-gcd extraction.

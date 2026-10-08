# Common-twist trigonometric determinants are scalar-character products

Three balanced source characters give an exact signed determinant
factorization. After the forced rational factor is removed, it is the
product of three ordinary character imaginary coordinates. Thus this
determinant, bounded through those coordinates and their phase widths,
recovers exactly the product of their scalar bounds. This is a limitation
of this determinant argument, not a general exclusion of additive methods
or a resolution of the uniform lattice-arc problem.

The final section extends the identity to every odd number of characters
and the corresponding full trigonometric Vandermonde determinant.

## Literal conventions and theorem

Fix the Gaussian primes and allocation exponents of actual equal-norm
source rows as in [the pair factorization note](six_row_common_twist_determinant_factorization.md).
For any integer row character `L` with coordinate sum zero, set

```
c_L(p)=sum_i L_i a_i(p),
A_L=product_p pi_p^max(c_L(p),0) bar(pi_p)^max(-c_L(p),0).
```

No separate units or sign normalizations are inserted into `A_L`.
Take three distinct balanced sign characters `lambda_1,lambda_2,lambda_3`
on the same even number of rows. Put `c_r(p)=c_(lambda_r)(p)`. These three
integers have a common parity `kappa_p` at each prime. Define

```
t=product_p p^kappa_p,
h_r=product_p p^[(|c_r(p)|-kappa_p)/2],
A_r=x_r+i y_r,             Norm(A_r)=t h_r^2,
L_rs=(lambda_r-lambda_s)/2, B_rs=A_(L_rs),       r<s.
```

Then, in the displayed row order,

```
D := det [(h_r,x_r,y_r)]_(r=1,2,3)
   = -4 F product_(r<s) Im B_rs,                       (1)

F = product_p p^f_p,
f_p = [sum_r |c_r(p)|
       -(max_r c_r(p)-min_r c_r(p))-kappa_p]/2.         (2)
```

Every exponent `f_p` is a nonnegative integer. These identities allow
arbitrary nested allocations at a prime. They also hold when a character
specializes to a unit, in which case the determinant may vanish.

## Proof, including the sign

Choose real arguments `theta_p` for the fixed `pi_p`, and use the
unreduced real lifts `phi_r=sum_p c_r(p)theta_p`. Then, literally,

```
A_r/(sqrt(t)h_r)=exp(i phi_r),
B_rs/|B_rs|=exp(i(phi_r-phi_s)/2).
```

Thus no independent square-root sign is chosen for a pair. The oriented
three-point circle determinant is

```
det [(1,cos phi_r,sin phi_r)]
 = -4 product_(r<s) sin((phi_r-phi_s)/2).             (3)
```

Factoring `h_r` from each row and `sqrt(t)` from the last two columns
of `D`, (3) gives (1) with

```
F = t h_1 h_2 h_3 / product_(r<s)|B_rs|.              (4)
```

At one prime, write `q_rs=(c_r-c_s)/2`. The three-variable identity
`sum_(r<s)|q_rs|=max c_r-min c_r` turns the local exponent of (4)
into (2). For three real numbers,

```
sum |c_r|-(max c_r-min c_r)
 = |median(c_1,c_2,c_3)|
   +2 dist(0,[min c_r,max c_r]).
```

This is nonnegative, has parity `kappa_p`, and is at least one when
`kappa_p=1`. Therefore (2) is a nonnegative integer. This proves both
the signed identity and its exact rational factor without discarding
content or inserting Gaussian units.

## Nonvanishing and the resulting height bound

Suppose the actual distinct source points have one common Gaussian unit
and lie on an arc of length at most `C sqrt(R)`, with `C<=2`.
[Affine allocation rigidity](linear_allocation_affine_rigidity.md) implies
that each nonzero `L_rs` has a nonzero allocation exponent somewhere.
Thus `B_rs` is a conjugate-primitive Gaussian nonunit of odd norm. It
cannot be real or purely imaginary: a nonunit on either axis is divisible
by both primes of some conjugate pair. In particular `Im B_rs` is a
nonzero integer and `D!=0`.

Let `Delta` be the angular width and `ell_rs=||L_rs||_1`. Common source
units cancel in the balanced row characters, and the source angles give

```
dist(arg B_rs,pi Z) <= ell_rs Delta/4,
1 <= |Im B_rs| <= |B_rs| ell_rs Delta/4.              (5)
```

Combining (1) with its phase upper bound, or simply multiplying (5),
gives exactly

```
product_(r<s)|B_rs| >= 4^3/[Delta^3 product_(r<s)ell_rs],

(1/2) sum_(r<s) log Norm(B_rs)
 >= 3 log(4/Delta)-sum_(r<s)log ell_rs.               (6)
```

The forced factor cancels completely through (4), and
`D/(4F)=-product Im B_rs` retains any additional evaluated divisibility
inside the three scalar factors themselves. Thus the determinant supplies
exactly their product bound, with no improved radius exponent from the
common twist. For `Delta<=C/sqrt(R)`, (6) has leading term
`(3/2)log R`, the sum of the three ordinary scalar terms.

## Every odd trigonometric Vandermonde determinant

Take `m=2d+1` distinct balanced sign characters with the same literal
conventions, where `d>=1`. The integer matrix whose row `i` is

```
(h_i^d,
 h_i^(d-1) Re A_i, h_i^(d-1) Im A_i,
 ...,
 h_i^(d-r) Re A_i^r, h_i^(d-r) Im A_i^r,
 ...,
 Re A_i^d, Im A_i^d)
```

has determinant `D_d` satisfying

```
|D_d| = 2^(2d^2) F_d product_(i<j)|Im B_ij|,          (7)
B_ij=A_((lambda_i-lambda_j)/2),

F_d=product_p p^E_p,
E_p=(d/2)sum_i|c_i(p)|
    -(1/4)sum_(i<j)|c_i(p)-c_j(p)|
    -kappa_p d^2/2.                                  (8)
```

Every `E_p` is a nonnegative integer. To see this, the pairwise absolute
value identity reduces (8) to

```
2E_p = sum_(i<j : c_i c_j>0) min(|c_i|,|c_j|)
        -kappa_p d^2.                               (9)
```

If `kappa_p=0`, all terms in the sum are even and nonnegative. If
`kappa_p=1`, all exponents are nonzero odd integers. With `u` of them
positive, there are

```
binom(u,2)+binom(2d+1-u,2)
 = d^2+(u-d)(u-d-1)
```

same-sign pairs. This number is at least `d^2`, and differs from `d^2`
by an even integer. Each minimum is positive and odd. Consequently (9)
is an even nonnegative integer in both parity cases.

For the determinant proof, factor `h_i^d` from row `i` and the powers
of `sqrt(t)` from the trigonometric column pairs. Their combined factor
is `t^[d(d+1)/2] product_i h_i^d`. The ordinary complex Vandermonde on
`exp(i phi_i)`, converted to the real columns
`1,cos phi,sin phi,...,cos dphi,sin dphi`, has modulus

```
2^(-d) product_(i<j)|exp(i phi_i)-exp(i phi_j)|
 = 2^(2d^2) product_(i<j)|sin((phi_i-phi_j)/2)|.
```

This proves (7), with the exact positive factor

```
F_d = t^[d(d+1)/2] product_i h_i^d / product_(i<j)|B_ij|,
```

whose prime exponents are (8). Under the same actual short-arc hypotheses,
all pair characters are nonunits, and multiplying their inequalities (5)
gives precisely the determinant's height bound:

```
product_(i<j)|B_ij|
 >= (4/Delta)^binom(m,2) / product_(i<j)||L_ij||_1.
```

Thus the whole displayed odd trigonometric determinant family, including
its value-dependent integer coefficients involving the `h_i`, factors
into ordinary scalar character coordinates. This does not assert a
no-go theorem for other additive identities, other determinants, or extra
arithmetic information about the evaluated scalar factors.

## Exact algebra audit

[The exact checker](check_common_twist_triple_determinant_factorization.py)
verifies the signed identity and norm/content identity on all bounded
single-prime equal-parity exponent triples and on several-prime fixtures.
It also checks every triple of the ten standard balanced six-row
characters on deterministic nested allocation fixtures, the general
absolute determinant identity at `d=1,2`, and the local exponent parity
argument for those sizes. These algebra fixtures do not assert that the
source points lie on a short arc.

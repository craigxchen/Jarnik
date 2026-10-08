# A finite sector obstruction for the nonorthogonal twelve-factor pattern

This note turns the four-coordinate certificate from
[the critical moment-code note](critical_pell_moment_codes.md) into a
finite Gaussian-integer inequality. The factors may vary freely with the
radius. No fixed Pell template or asymptotic expansion is used.
The sector and prefactor conditions below are substantive; they have not
been established for arbitrary endpoint configurations.

## 1. Exact hypotheses

Let S be the eight-by-twelve matrix in equation (15) of the linked note,
with indices starting at zero. Let H_j be nonzero Gaussian integers and
let E_i be nonzero Gaussian integers with a common modulus. Form

\[
z_i=E_i\prod_{j:S_{ij}=1}H_j
              \prod_{j:S_{ij}=-1}\overline{H_j}.
\]

These points have a common radius R. Write

\[
c=(-1,-1,-1,0,2,2,0,-1),\qquad
cS=6(e_1+e_2-e_9-e_{10}).
\tag{1}
\]

Assume:

* The z_i lie in an arc of angular width Delta < pi/4.
* The four factors H_1,H_2,H_9,H_10 have argument lifts in a common
  interval of width omega < pi/12.
* The Gaussian rational number product_i E_i^(c_i) is positive real.
  This holds automatically if all E_i are one common Gaussian unit.
* The Gaussian integer
  B = H_1 H_2 conjugate(H_9) conjugate(H_10) is not real.

One sufficient arithmetic condition for the last hypothesis is that the
four selected factors have pairwise coprime rational norms, each is
coprime to its own conjugate, and at least one is not a unit. Indeed, a
Gaussian prime dividing a selected nonunit occurs in only one orientation
in B, whereas a real Gaussian integer has equal valuations at conjugate
primes.

## 2. The finite lower bound

Put P = |B| = |H_1 H_2 H_9 H_10|. Then

\[
\boxed{\Delta\geq\frac{3}{2P}.}
\tag{2}
\]

To prove this, choose argument lifts phi_i for the points in one interval
of width Delta and lifts theta_j for the four selected factors in their
sector. Set

\[
\Theta=\theta_1+\theta_2-\theta_9-\theta_{10}.
\]

Since the positive and negative coefficient masses of c both equal four,
sum_i c_i = 0 gives

\[
\left|\sum_i c_i\phi_i\right|\leq4\Delta<\pi.
\tag{3}
\]

The equal positive and negative masses in Theta give
|Theta| <= 2 omega, so |6 Theta| < pi. The multiplicative identity (1)
and the positive-real prefactor condition imply

\[
\sum_i c_i\phi_i\equiv6\Theta\pmod{2\pi}.
\]

Both sides belong to (-pi,pi), so they are equal. Hence
|Theta| <= 2 Delta/3. Finally B is a nonreal Gaussian integer, and thus

\[
1\leq|\operatorname{Im}B|
 =P|\sin\Theta|\leq P|\Theta|
 \leq\frac23P\Delta.
\]

This proves (2). The branch comparison above is needed; reducing the
angle identity modulo 2 pi without it would not justify (2).

## 3. Consequences at the endpoint scale

If the arc length is at most C sqrt(R), its angular width is at most
C/sqrt(R). Whenever the hypotheses hold, (2) gives

\[
\boxed{\sqrt R\leq\frac{2C}{3}P.}
\tag{4}
\]

Writing W = log(R^2) and W_sel = sum_{j in {1,2,9,10}} log N(H_j),
the same inequality is

\[
\boxed{W_{\rm sel}\geq\frac W2+2\log\frac{3}{2C}.}
\tag{5}
\]

If the circle configuration is primitive, W is its conductor weight.
Without primitivity it is the total log norm; the geometric inequality
(4) still holds as stated.

In particular, a uniform bound P <= A R^(1/3), with A fixed, yields

\[
R\leq(2CA/3)^6.
\tag{6}
\]

Thus this pattern with comparable-size factors in the specified sector
is excluded above an explicit radius, even if the factors are selected
independently for each circle. More generally, P <= A R^eta with
eta < 1/2 gives R <= (2CA/3)^(1/(1/2-eta)).

## 4. A stable form of the reciprocal collision

There is also an error-tolerant algebraic version of the four-term
reciprocal lemma. For nonzero real a,b,c,d, put

\[
u=a+b-c-d,\qquad
v=a^{-1}+b^{-1}-c^{-1}-d^{-1}.
\]

The exact identity

\[
\boxed{(a+b)(a-c)(a-d)=a^2u+abcd\,v}
\tag{7}
\]

follows by multiplying out. If the four absolute values are pairwise
separated by at least delta > 0 and all are at most H, then

\[
\boxed{\delta^3\leq H^2|u|+H^4|v|.}
\tag{8}
\]

Every factor on the left of (7) has absolute value at least delta by the
reverse triangle inequality. This proves (8), including arbitrary signs
and near-cancellation of a+b. In particular, if both residuals are at
most epsilon, then epsilon >= delta^3/(H^2+H^4).

This estimate requires bounds on both residuals. A small value in one
embedding of a quadratic field does not by itself bound its conjugate;
no such implication is assumed here.

[ReciprocalCollision.lean](../GaussianChain/ReciprocalCollision.lean)
checks the cubic identity and the quantitative norm estimate as
`reciprocal_cubic_identity` and `stable_reciprocal_bound`. The latter
exposes the three lower bounds on |a+b|, |a-c|, |a-d| and the four upper
bounds on the individual absolute values. The Gaussian sector argument
and the circle-to-residual reduction are not formalized by that module.

## 5. A weaker bound without any sector restriction

Keep the prefactor condition and B nonreal, but drop both the sector
condition and the restriction Delta < pi/4. Then

\[
\boxed{\Delta\geq\frac{3}{8P^2}.}
\tag{9}
\]

Indeed, (1) and (3) still show that the distance from arg(B) to a
multiple of pi/3 is at most 2 Delta/3. Write B = x+iy with integers x,y.
If a closest multiple is zero or pi, nonreality gives |y| >= 1, so the
angular distance is at least 1/P. At any other multiple, the function
sin(theta)^2 - 3 cos(theta)^2 vanishes. Its derivative has absolute value
at most four, while its value at theta = arg(B) has absolute value

\[
\frac{|y^2-3x^2|}{P^2}\geq\frac1{P^2}.
\]

The integer numerator is nonzero because sqrt(3) is irrational and B
is nonzero. The mean value theorem therefore bounds this angular
distance below by 1/(4P^2). This also bounds the axis case, since P >= 1.
Combining the distance bounds proves (9). Only the inequality
|sum_i c_i phi_i| <= 4 Delta from (3) is used here, not its final
comparison with pi.

At endpoint scale the resulting finite constraints are

\[
\sqrt R\leq\frac{8C}{3}P^2,\qquad
\boxed{W_{\rm sel}\geq\frac W4+\log\frac{3}{8C}.}
\tag{10}
\]

Thus this exact pattern cannot place four selected factors of total
log norm substantially below one quarter of the total into an endpoint
cluster, irrespective of their directions. The stronger half-weight
constraint (5) requires the sector assumption.

## Remaining gap

The result supplies a finite weight constraint for one nonorthogonal
pattern. It does not assert that a large arbitrary cluster contains this
pattern or that its selected factors have less than one quarter of the
total weight. To use the stronger half-weight conclusion, a common sector
is still required. The result also does not remove an arbitrary
private prefactor phase. The unrestricted uniform exponent-1/2 theorem
and the general intrinsic eight-point inequality remain unproved.

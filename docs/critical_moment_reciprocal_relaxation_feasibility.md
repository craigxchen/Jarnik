# The first-moment and reciprocal cross-Gram relaxation is feasible

There exists a finite eight-row sign matrix `S`, with `n=1790=4s+2`
and `s=447`, and positive pairwise distinct real magnitudes `a_h`, such
that adding a constant column and a balanced column makes the rows
orthogonal, and

```
S a=0,
beta_i=sum_h S_ih/a_h,       gamma=sum_h a_h^-2,
C_ik=sum_h S_ih S_kh/a_h^2
     =gamma-(beta_i-beta_k)^2/2           (i in A, k in B).
```

All eight `beta_i` are distinct. The within-group reciprocal distance
expressions are positive, agreeing with the necessary signs for this
matrix's balanced differing supports. The full within-group residual
identities are not included in this assertion.

This is feasibility of a **weaker system**. It does not establish the
required critical/successor pair root orders, higher odd moments,
within-group residual-polynomial formulas, or the scalar cross-product
identities. Nor does it identify any prescribed values with the actual
middle polynomial coefficients. In fact, the final section proves that
every cross support violates the required paired sign order and has
nonzero fourth reciprocal elementary symmetric function. Thus this is
neither a full moment-template solution nor a circle construction.

## A positive first-moment base

Take all `256` columns `v in {+1,-1}^8` once, and delete the two columns
`1` and `b=(1,1,1,1,-1,-1,-1,-1)^t`. Write the resulting `8 x 254`
matrix as `S_0`. The two row groups are `A={0,1,2,3}` and
`B={4,5,6,7}`. Sign-cube orthogonality gives

```
G=S_0 S_0^t=256 I-11^t-bb^t,
S_0 1=-1-b,
a*=1+S_0^t(1+b)/248,        S_0 a*=0.
```

Every entry of `a*` lies between `1-2/31` and `1+2/31`. Perturb it
inside `ker S_0` by less than `1/100` in supremum norm to obtain a
rational vector `a^0` with all entries distinct. Such a perturbation
exists: no coordinate-equality functional `e_j-e_k` vanishes on
`ker S_0`. Indeed, its squared projection onto the row space is

```
(v_j-v_k)^t G^-1(v_j-v_k) <=32/248<2=||e_j-e_k||^2.
```

Here the smallest eigenvalue of `G` is `248`. Thus each forbidden
coordinate equality removes only a proper affine hyperplane; finitely
many such hyperplanes cannot exhaust a rational open ball in the kernel.
In particular `9/10<a_h^0<11/10`.

Define the base reciprocal sums `beta^0`, diagonal value `gamma_0`,
and full reciprocal Gram matrix `C^0`. The elementary bounds are

```
gamma_0<314,       |beta_i^0|<283,       |C^0_ij|<314.  (1)
```

## First-moment-neutral triads

For a sign column `v` and positive `p,q`, append the three columns
`v,v,-v` with respective magnitudes `p,q,p+q`. Their first-moment
contribution is zero. Their reciprocal linear and Gram contributions are

```
t v,       t^2 vv^t,
t=1/p+1/q-1/(p+q)>0.
```

The identity in the quadratic contribution is
`t^2=1/p^2+1/q^2+1/(p+q)^2`. Reversing all three column signs realizes
the negative coefficient `-t` without changing the Gram contribution.
Their unweighted count Gram contribution is always `3vv^t`.

Every nonzero real coefficient `t` can be realized, with a free shape
parameter `rho in (1,2)`, by taking

```
p=(1+rho+rho^2)/(rho(1+rho)|t|),       q=rho p
```

and orienting the three columns by `sign(t)`.

## Exact completion of the cross Gram entries

Fix `L=10000` and `gamma=2L^2`. Prescribe the distinct reciprocal sums

```
beta_i=i/L                         (i=0,1,2,3),
beta_(4+j)=2L+j/L                  (j=0,1,2,3).
```

For every cross pair prescribe
`C*_ik=gamma-(beta_i-beta_k)^2/2`. These sixteen rational numbers
have absolute value less than `7`. Put `Delta_ik=C*_ik-C^0_ik`.
For each of the `256` sign columns `v`, define

```
W_v=[gamma-gamma_0+sum_(i in A,k in B) Delta_ik v_i v_k]/256,
T_v=[(beta-beta^0) dot v]/256.                         (2)
```

The elementary sign-cube Fourier sums imply exactly

```
sum_v W_v=gamma-gamma_0,
sum_v T_v v=beta-beta^0,
sum_v W_v v_i v_k=Delta_ik             (i in A,k in B),
sum_v W_v v_i v_j=0                   (distinct i,j in one group).
```

The bounds (1) give the uniform estimates

```
W_v>(2L^2-5450)/256,
|T_v|<(8L+12/L+2264)/256.
```

At `L=10000` the former lower bound exceeds the square of the latter
upper bound. Therefore `W_v>T_v^2` for every `v`. The two real numbers

```
t_(v,+)=[T_v+sqrt(2W_v-T_v^2)]/2,
t_(v,-)=[T_v-sqrt(2W_v-T_v^2)]/2
```

are nonzero, of opposite sign, and have sum `T_v` and sum of squares
`W_v`. Append one triad for each. Equations (2) now prove all the
claimed first-moment, reciprocal-sum, diagonal, and cross-Gram equations
exactly.

There are `512` triads, hence `1536` appended columns. Their count Gram
is `6 sum_v vv^t=1536 I`. The completed matrix consequently has

```
S S^t=1792 I-11^t-bb^t,       n=254+1536=1790.
```

Adding columns `1,b` therefore gives the required eight orthogonal rows.
Each triad shape parameter can be chosen successively to avoid all
previous magnitudes: each of `p(rho),q(rho),p(rho)+q(rho)` is a
nonconstant rational function on `(1,2)`, and so takes any forbidden
value at only finitely many parameters. Inside a triad, `p<q<p+q`.
This proves pairwise distinctness without changing any verified equation.

## The within-group distance signs also hold

The appended within-group off-diagonal Gram entries vanish. Thus for
distinct rows in one group the resulting expression is

```
4L^2-2C^0_ij-(beta_i-beta_j)^2.
```

Within either group, the squared reciprocal-sum difference is at most
`9/L^2`, so this expression is positive by (1).

Indeed, `t_(v,+)>0` and `t_(v,-)<0`, so the two oriented triads at each
`v` have cancelling column sums. Thus the completed row sums remain
`S1=-1-b`: they are constant on each entire group. Every within-group
differing support consequently has equally many positive and negative
signed roots. In a full solution all such pairs would have successor
type A and positive distance expression. Thus the computed signs match
those necessary signs for this actual sign matrix. This does not assert
the full successor residual formulas or pair root orders.

## Every cross support fails the next required conditions

The failure of full-moment admissibility can be proved at the stated
finite `L`, independently of the free triad shapes. Put

```
Q=gamma-gamma_0,       D=sum_cross |Delta_ik|,
M_0=max_h (a_h^0)^-2.
```

For a fixed cross support, let `M_4` be its base contribution to the
fourth reciprocal moment, and write `p_l=sum z_h^-l` for its full
signed reciprocal power sums. Equation (2) gives

```
(Q-D)/256 <= W_v <= (Q+D)/256.
```

Exactly `128` sign vectors `v` differ on this cross pair. All three
columns of either associated triad then belong to its differing
support. Each triad has reciprocal-square mass `t^2`; each individual
reciprocal square is at most `t^2`, and its fourth-moment mass is at
most `t^4`. Consequently

```
p_2 >= (Q-D)/2,
q_max^2 <= max(M_0,(Q+D)/256),
p_4 <= M_4+(Q+D)^2/512,                              (3)
```

where `q_max` is the largest reciprocal magnitude on this support.
The first-moment base bounds give, uniformly over the sixteen supports,

```
Q>2L^2-314,       D<5136,       M_0<2,       M_4<388.
```

At `L=10000` these imply `Q>=max(2D,16M_0)` and `Q^2>=512M_4`.
In particular, (3) yields

```
p_2>=Q/4>=4q_max^2.
```

But a critical paired sign order, together with the already satisfied
identity `p_1^2=p_2`, requires `p_2<4q_max^2`, as proved in
[the reciprocal concentration note](critical_moment_reciprocal_concentration.md).
Thus every cross pair fails that order.

There is also a direct algebraic failure. From (3),

```
rho=p_4/p_2^2 <=16M_4/Q^2+9/128 <=13/128<1/9.
```

With `e_2=0`, Newton's formula reads
`12e_4=4p_1p_3-p_2^2-3p_4`. Cauchy--Schwarz and `p_1^2=p_2`
give `|p_1p_3|<=p_2 sqrt(p_4)`. It follows that

```
12e_4 <=p_2^2(4sqrt(rho)-1-3rho)<0,
```

since `4sqrt(rho)-1-3rho=-(1-sqrt(rho))(1-3sqrt(rho))`.
Every cross support therefore has `e_4<0`, whereas its full critical
residual would require `e_4=0`.

The [exact checker](check_critical_moment_reciprocal_relaxation_feasibility.py)
verifies the finite sign-cube sums, an explicit distinct rational base,
the uniform bounds, triad identities, count Gram, and finite thresholds
in this section. It certifies the construction's algebra and existence
bounds; it does not numerically materialize the algebraic square roots.

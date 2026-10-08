# Four-row nearest-square Gram rigidity

Four simultaneously small nearest-square defects make the positive Gram
matrix rigid.  In the full four-row Boolean profile, the norm, bracket, and
metric equations then supply an **integral near-axis realization** without
separately testing the three Gaussian divisibilities in the one-coordinate
reduction.  The new realization can have integer content.  Removing that
content and trimming only the affected prime-power depths recovers a
primitive balanced profile with subpower losses.  This is a conditional
sufficiency statement: the four small defects have not been proved from the
metric equations or from the first norm alone.

## 1. A finite matrix criterion

Let `P_i=(X_i,Y_i) in Z^2`, `1<=i<=4`, have rank two, and put

```text
N_i=|P_i|^2,       S_ij=P_i dot P_j,
Delta_ij=det(P_i,P_j),
a_i=floor(sqrt(N_i)),       r_i=N_i-a_i^2,
B_ij=S_ij-a_i a_j,          E=max(1,max_(i,j)|B_ij|).
```

All these numbers are integers.  Let `G=gcd_Z[i](P_1,...,P_4)`,
`R_i=P_i/G`, and

```text
D_i=gcd_Z[i](R_j : j!=i),       L=min_i |D_i|.            (1)
```

Gaussian gcds are specified only up to units; their moduli are
unambiguous.  Since `gcd_Z[i](R_1,...,R_4)=1`, one has
`gcd_Z[i](D_i,R_i)=1`.  Therefore every nonzero Gaussian-integer relation
`sum_i c_i P_i=0` has

```text
max_i |c_i| >= L.                                        (2)
```

**Finite rigidity theorem.** If

```text
6 E^3 < L,                                              (3)
```

then there are integers `b_i` satisfying

```text
r_i=b_i^2,       S_ij=a_i a_j+b_i b_j,
Delta_ij=a_i b_j-a_j b_i.                               (4)
```

The choice of all `b_i` is simultaneous; replacing every `b_i` by its
negative reverses all brackets, so (4) uses the orientation of the given
rows.  In particular `a_i+i b_i` is an integral realization of exactly
the same oriented Gram data.  The exact divisibilities

```text
a_1-i b_1 divides S_1j+i Delta_1j,       j=2,3,4,       (5)
```

follow from that realization.

To prove the theorem, write `H` for the `2 by 4` integer matrix whose
columns are `P_i`.  Then `S=H^t H` has rank two, and
`B=S-aa^t` has rank at most three.  If `rank(B)=3`, choose three
independent rows of `B` and take their signed `3 by 3` cofactors.
This gives a nonzero integer `c in ker B` with
`max|c_i|<=6E^3`.  Here `a` cannot belong to `im S`, since otherwise
`im B subset im S` and `rank B<=2`.  But
`Sc=a(a dot c)` lies in `im S`; hence `a dot c=0`, `Sc=0`, and
`Hc=0`.  This contradicts (2)-(3).

If `rank(B)=2`, then `a in im S`: if `a` had a component outside that
two-dimensional space, the symmetric rank-one update `S-aa^t` would
have rank three.  Consequently `im B=im S` and `ker B=ker S`.
Two independent rows of `B` yield a nonzero three-supported cofactor
vector in `ker B`, with entries at most `2E^2<=6E^3`.
It again satisfies `Hc=0`, contradicting (2)-(3).

Thus `rank(B)<=1`.  Rank zero would give `S=aa^t`, contrary to
`rank S=2`.  Rank one implies `a in im S`; write `a=H^t u` for
`u in Q^2` (two independent integer columns of `H` determine `u`).
Then

```text
B=H^t (I-u u^t) H.
```

Since `H` maps `R^4` onto `R^2`, `rank B=rank(I-u u^t)=1` forces
`u dot u=1`.  Let `v=(-u_2,u_1)`.  The matrix with rows `u^t,v^t`
is a rational rotation of determinant one, and
`B=H^t v v^t H`.  Its second coordinates `b_i=v dot P_i` are
rational, while `b_i^2=B_ii=r_i` are integers.  A rational number
whose square is an integer is itself an integer.  This proves (4),
and `(S_1j+i Delta_1j)/(a_1-i b_1)=a_j+i b_j` proves (5).

## 2. The balanced Boolean threshold

Assume the original `P_i` are conjugate-primitive Gaussian rows with
the complete fifteen-block factorization

```text
P_i=K_i product_(T contains i) H_T,
n_T=Norm(H_T),
Delta_ij=t_ij product_(T contains i,j)n_T !=0,
log n_T=w+o(w),
log Norm(K_i), log max(1,|t_ij|)=o(w),
```

where the odd split-core norms have pairwise disjoint rational-prime
supports.  These original rows need not be near the real axis.  Let
`S_ij>0` for all pairs, as in the endpoint chamber, and suppose

```text
max(1,max_i r_i)=exp(o(w)).                             (6)
```

Then `N_i=exp(8w+o(w))`, `a_i=exp(4w+o(w))`, and
`|Delta_ij|=exp(4w+o(w))`.  For `i!=j`, the exact identity

```text
(S_ij-a_i a_j)(S_ij+a_i a_j)
 =a_i^2 r_j+a_j^2 r_i+r_i r_j-Delta_ij^2               (7)
```

has denominator `S_ij+a_i a_j=exp(8w+o(w))`.  Its numerator is at most
`exp(8w+o(w))`; hence every `|B_ij|<=exp(o(w))`,
including `B_ii=r_i`.  Thus `E=exp(o(w))`.

The co-singleton gcds in (1) have `L>=exp(w/2-o(w))`.
Indeed `H_([4]\{i})` divides the other three rows.  Its part absorbed
by their common gcd `G` also divides `P_i`; by disjoint core supports
that part divides `K_i`.  Hence

```text
|D_i| >= |H_([4]\{i})|/|K_i|=exp(w/2-o(w)).              (8)
```

For large `w`, (3) follows.  The finite theorem forces every defect
in (6) to be a square and gives the integral near-axis rows (4).
In particular, once all four nearest-square defects are small,
the first-row candidate and the three Gaussian congruences of (5)
hold automatically.  Actual extracted near-axis profiles already
satisfy (6).  To obtain a four-row obstruction for a uniform endpoint
bound, one would need to **exclude** (6) for large admissible metric/core
data, or show that it forces a correction or determinant residual to
have fixed positive exponent.  An isolated small first defect does not
imply (6) without the divisibility or a near-axis integral realization.

## 3. Primitive recovery after a small content loss

Put `Q_i=a_i+i b_i` from (4).  Every `max(1,|b_i|)=exp(o(w))`.
In fact `b_i!=0` for all sufficiently large `w`.  If `b_i=0`, then
`N_i=a_i^2` and, for any `j!=i`,
`Delta_ij=a_i b_j`.  The singleton norm `n_{ {i} }` divides `N_i`
and is coprime to the pair core
`G_ij=product_(T contains i,j)n_T`.  Squaring the determinant identity
and using `Delta_ij=G_ij t_ij` gives
`n_{ {i} } | t_ij^2`, contrary to
`log n_{ {i} }=w+o(w)` and `log|t_ij|=o(w)`.

Let `c_i=gcd(a_i,|b_i|)` and `C=product_i c_i`.
Then `c_i<=|b_i|=exp(o(w))`.  The rows `Q_i/c_i` have primitive
integer coordinates and still have subpower imaginary coordinates.
To retain the entire Boolean core, define at each **prime-power depth**

```text
n'_T=n_T/gcd(n_T,C^2).                                  (9)
```

The `n'_T` remain pairwise coprime odd split norms, with
`log n'_T=w+o(w)`.  For each incident row,
`n'_T | Norm(Q_i/c_i)`; for each pair `i,j in T`,
`n'_T | det(Q_i/c_i,Q_j/c_j)`.  At a prime `p^e||n_T`, this follows
from the explicit retained exponent
`max(0,e-2v_p(C))`, which is at most both
`e-2v_p(c_i)` and `e-v_p(c_i)-v_p(c_j)` whenever positive.

The primitive orientation argument of the
[one-coordinate divisor reduction](four_row_one_coordinate_divisor_reduction.md)
now supplies Gaussian factors `H'_T` of norms `n'_T`, coherently dividing
every incident new row.  The remaining correction norms and pair
determinant residuals have logarithms `o(w)`: each lost block factor
has norm at most `C^2`, and there are only fifteen blocks.  This
recovers an actual primitive balanced profile, though its exact norms
and blocks can differ from the input by subpower factors.

## 4. Why four columns enter the finite argument

Three columns admit bounded nonsquare defects while every nonzero ordinary
integer relation has unbounded height.  This concerns the abstract Gram
argument; it does not supply a balanced Boolean profile.

Let `U^2-2V^2=1`, with positive Pell solutions beginning at `(3,2)`.
The integral matrix

```text
J = [ (U+1)/2  (U-1)/2  V ]
    [ (U-1)/2  (U+1)/2  V ]
    [ V        V        U ]
```

preserves `x^2+y^2-a^2`.  Apply it to the columns
`(4,1,2)`, `(4,3,2)`, `(5,4,2)`, and retain their first two
coordinates as `P_1,P_2,P_3`.  Their last coordinates are
`a=(5V+2U,7V+2U,9V+2U)`.  The Gram defect is constant:

```text
S-aa^t = [ 13  15  20 ]
         [ 15  21  28 ]
         [ 20  28  37 ],       det(S-aa^t)=-16.
```

The diagonal defects `13,21,37` are nonsquares and are smaller than
`2a_i+1`, so each `a_i` is exactly `floor(sqrt(N_i))`.
The three oriented minors are
`Delta_12=8U+4V`, `Delta_13=11U+4V`, `Delta_23=U`.
Their gcd is one: `U` is odd and coprime to `V`.  Consequently the
primitive integer relation is
`(U,-11U-4V,8U+4V)`, whose height tends to infinity.

Thus bounded Gram error plus absence of a bounded ordinary linear
relation cannot force integral rotation with only three columns.  No
claim about the Boolean gcd threshold `L` is made for this family.

The [exact checker](check_four_row_nearest_square_gram_rigidity.py)
verifies the matrix rank/cofactor identities, an example satisfying
`6E^3<L` with a nontrivial rational rotation, prime-power trimming,
and twelve members of the three-column Pell family.  The finite checks
do not replace the cofactor or primewise proofs above.

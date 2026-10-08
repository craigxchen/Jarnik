# Exact integrality and coordinate barriers for a pairwise Vieta descent

For the full four-row profile, the alternate positive binary form through
two prescribed norm values is not generally integral. Its integrality has
an exact criterion: the actual primitive pair imaginary residue must be
`+/-1` or `+/-2`, and the actual Gaussian pair-gcd norm must divide the
primitive pair real coordinate. This retains all correction primes.
In a positive endpoint chamber the alternate form also increases the
leading coefficient. An integral change of binary coordinates that
preserves small second coordinates cannot reduce that coefficient.

These are barriers to this specified descent, not a proof excluding the
four-row system. In particular, they do not rule out a coordinated change
of several block norms or a different arithmetic construction.

## 1. The actual integral form and pair data

Use the full profile and its exact norm interpolation from
[the many-row norm system](many_row_norm_barycentric_system.md).
Let `H=H_[4]=A+iB`, `g=Norm(H)`, and choose

```text
r=A B^(-1) mod g,     s=(r^2+1)/g,
a_i=(X_i-rY_i)/g.
```

The full block is odd and conjugate-primitive, so these are integers.
The binary vectors `v_i=(a_i,Y_i)` are primitive: a common integer
divisor would divide both `X_i=g a_i+rY_i` and `Y_i`. The actual
Gaussian quotient and its norm are

```text
Q_i=P_i/H=bar(H) a_i+((r+i)/H)Y_i,
f(a,y)=g a^2+2ray+s y^2,       gs-r^2=1.                 (1)
```

The two Gaussian coefficients in (1) form an integral oriented basis
of determinant one. Thus every `Q_i` is conjugate-primitive, and

```text
delta_ij=det(v_i,v_j)=Im(bar(Q_i)Q_j)=Delta_ij/g,
S_ij=Re(bar(Q_i)Q_j) in Z.
```

Let `D_ij` be the norm of the actual Gaussian gcd of `Q_i,Q_j`.
The Gaussian integer

```text
u_ij+i v_ij=bar(Q_i)Q_j/D_ij                             (2)
```

is conjugate-primitive. Indeed, after removing the common Gaussian
gcd, a split prime can occur in the product in only one orientation:
equal original orientations cancel to their valuation difference, and
opposite original orientations reinforce the same product orientation.
Consequently

```text
gcd(u_ij,v_ij)=1,     D_ij odd,
S_ij=D_ij u_ij,       delta_ij=D_ij v_ij,
v_ij!=0.
```

This uses the actual gcd, not merely the prescribed common core.
More explicitly, put

```text
G_ij=product_(T containing i,j, T!=[4])n_T,
k_i=Norm(K_i).
```

Then

```text
D_ij=G_ij e_ij,       e_ij in Z_(>0),
e_ij divides k_i k_j,   t_ij=e_ij v_ij.                 (3)
```

For (3), a shared core prime has no extra gcd valuation beyond the
two correction valuations; at every other prime, at least one row has
only its correction valuation available. Disjoint core supports make
this a primewise proof. Thus `log e_ij=o(w)`, but `e_ij` must not
be discarded from an exact integrality test. In the balanced four-row
profile,

```text
log D_ij=3w+o(w),       log|v_ij|=o(w).                  (4)
```

## 2. The two forms through a prescribed pair

Fix `i,j`, and write `delta=delta_ij`, `S=S_ij`. The values
`f(v_i),f(v_j)` and determinant `det f=1` determine the mixed
value to be `S` or `-S`, since

```text
S^2=f(v_i)f(v_j)-delta^2.
```

When `S!=0`, these are exactly two distinct rational positive forms.
If `F=[[g,r],[r,s]]`, the other form is

```text
F'=F+(2S/delta^2) C,
C=[[2Y_iY_j, -(a_iY_j+a_jY_i)],
   [-(a_iY_j+a_jY_i), 2a_i a_j]].                       (5)
```

To see this, use the invertible matrix with columns `v_i,v_j`.
Its Gram matrix has diagonal entries `f(v_i),f(v_j)` and mixed
entry `S`. Replace that mixed entry by `-S` and change coordinates
back. Both Gram matrices are positive with determinant `delta^2`,
which proves positivity and determinant one in (5).

Let `c` be the gcd of the integer entries of `C`. Since both binary
vectors are primitive, the product
`(Y_i x-a_i y)(Y_j x-a_j y)` has coefficient content one.
It follows that `c` divides two. Its middle coefficient has the parity
of `delta`; hence

```text
c=1 if delta is odd,       c=2 if delta is even.          (6)
```

All three coefficients of `F'` are integers if and only if

```text
delta^2 divides 2c S.                                   (7)
```

This follows by an integer Bezout combination of the entries of `C`;
checking only the leading coefficient would give a weaker condition
that incorrectly allows factors of `Y_iY_j` to clear the denominator.

Substitute `delta=Dv` and `S=Du` from (2). Since `D` is odd
and `gcd(u,v)=1`, (7) is equivalent to the particularly simple test

```text
F' is integral  <=>  D divides u and |v| in {1,2}.        (8)
```

In fact (7) says `D v^2 | 2c u`. Its odd part forces `D|u`,
while coprimality forces `v^2|2c<=4`. Conversely, for `|v|=1`
one has `c=1`, and for `|v|=2` one has `c=2`; these conditions
and `D|u` suffice.

For four balanced rows, (8) requires an additional actual divisibility
of logarithmic size `3w+o(w)` in the primitive pair real coordinate,
as well as a residue of exactly one or two. The given subpower bound
on `t_ij`, even together with (3), does not by itself provide either
requirement. If `Y_i=Y_j=1`, the original odd primitive rows force
`delta` even, so the only possible primitive residue is `|v|=2`.

## 3. This obstruction actually occurs in a balanced family

In the genuine
[balanced three-row polynomial family](balanced_three_row_gaussian_polynomial_family.md),
remove the common full block `J_g=1+it`. The primitive pair imaginary
residues are `(-2,18,16)`. Thus pairs `13` and `23` fail (8)
immediately, for every progression value.

For the remaining pair `12`, the exact pair-gcd norm and real coordinate
are

```text
D=t^2-12t+61,
30u=t^4-19t^3+151t^2-559t+870.
```

The remainder of `30u` modulo the displayed quadratic is

```text
504-60t.
```

For every `t>=64` this remainder is nonzero and has absolute value
less than `D`. All values of that family's progression exceed 64.
Therefore `D` does not divide `30u`, and hence cannot divide `u`.
Pair `12` fails (8) as well.

Thus an actual unbounded balanced, prime-disjoint, conjugate-primitive
Gaussian family can have **no integral pairwise Vieta alternate at all**.
This three-row example does not decide whether an additional fourth
row forces a special pair satisfying (8). It does show that such a
four-row assertion would need a new simultaneous argument.

## 4. Leading-coefficient direction in an endpoint chamber

After choosing row signs so `X_i>0`, the small-imaginary four-row
profile has `S_ij>0` for all pairs once `w` is sufficiently large.
The leading coefficient of (5) is

```text
g'=g+4S_ij Y_iY_j/delta_ij^2.                            (9)
```

If the imaginary coordinates have the same sign, this increases `g`.
An actual endpoint extraction can choose its anchor at one end of the
selected short arc; after changing individual row signs, all `X_i`
and `Y_i` are then positive. Every pair alternate in that chamber
therefore increases the full-block norm `g`, even when (8) holds.
The increment has logarithm `w+o(w)`, so it creates no strict
block-height gap either.

For opposite-sign imaginary coordinates the difference in (9) is
negative. The statement is not that all changes after arbitrary
reanchoring increase `g`. Nor does the pair alternate automatically
preserve the other two prescribed norm values: three distinct binary
directions already determine a quadratic form uniquely.

## 5. Small-coordinate rigidity rules out a free Gauss-reduction step

Let a simultaneous rational coordinate change have second coordinate

```text
Y'_i=c a_i+dY_i.
```

For every pair the exact relation is

```text
Y_jY'_i-Y_iY'_j=c delta_ij.                              (10)
```

Put `E=max_i |Y_i|`, `E'=max_i |Y'_i|`. If `c` is a
nonzero integer, (10) gives

```text
E'>=|delta_ij|/(2E)=exp(3w-o(w)).                        (11)
```

Thus an integral unimodular change preserving subpower second
coordinates has `c=0`. Its diagonal entries are then `+/-1`, so
it is only a translation/shear of the first coordinate and sign
changes. Under the corresponding inverse change of the form, its
leading coefficient remains `g`.

More generally, if `q` is any positive common denominator of the
rational coordinate-change matrix and `c!=0`, then `|c|>=1/q`.
Preserving `E,E'=exp(o(w))` requires

```text
q>=|delta_ij|/(2EE')=exp(3w-o(w)).                       (12)
```

A rational reparametrization with subpower denominator cannot bypass
the integral rigidity. Gauss reduction of the positive unimodular
form is valid as a lattice operation, but a step reducing `g` must
pay either the large second-coordinate cost (11) or the denominator
cost (12). This makes explicit why a reduced form with bounded
coefficients does not bound the original marked small-imaginary rows.

The [exact checker](check_four_row_integral_vieta_descent_barrier.py)
checks (5)--(8) on integral positive forms and primitive odd-norm rows,
including integral and nonintegral alternates, and checks the balanced
family's exact remainder obstruction. It also verifies (10) on
integer and rational changes. No extra triple Heron equation is
counted as an independent condition.

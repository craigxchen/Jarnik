# A clean-prime frame can enlarge Vieta content

Irreducibility and avoidance of all boundary nodes do not make the
extra Vieta content a local resultant of the normalized quadratic
coefficients. A fixed section with coefficient sum `1884` has the
following exact family at the split prime `29`:

```text
v_29(g_i)=2m,       v_29(F_i)=m,
v_29(Norm(K_i))=0,  v_29(E_i)=m,
max |t_jk|<=6,     all v_29(b_jk)=0.                  (1)
```

Here all old Gaussian rows are conjugate-primitive, all old directions
are distinct, and every old row is unramified on the section. The
normalized row quadratic has coefficient gcd `120`, a `29`-adic
unit. All fifteen boundary-node coefficients and all five normalized
old gradient multipliers are also `29`-adic units, and the chart
discriminant is squarefree modulo `29`. The extra `29^m` in (1)
comes from the frame and a collision of the **other** root with a
retained row.

This is a local obstruction, not a counterexample to a global bound
`E_i<=C_Q^A exp(Bsigma) product b_jk^D`. In this family

```text
sigma=(m/2)log29+O(1),
log E_i=2m log29+O(1)=4sigma+O(1).                    (2)
```

Only one cut block grows; the other thirty nonempty blocks are units.
Thus the family also lacks the full-fair lower weight bounds. These
limitations are part of the result. It rules out charging every
clean-core extra valuation to a coefficient resultant or a vanishing
boundary slope, but leaves the proposed full-fair height inequality
open.

## 1. A fixed irreducible nonnodal section

In the integral six-graph basis of
[the gradient calculation](five_row_gradient_discriminant_core_count.md),
take coefficients

```text
(35,50,-4,-4,-4,-4).
```

Its original polynomial coefficient sum is `C_Q=1884`. In the chart
with rows `(0,1),(1,0),(1,1),(1,a),(1,b)`, it is

```text
F(a,b)=(81-27b)a^2+(19b^2-84b-39)a+15b^2+35b.        (3)
```

At `(a,b)=(2,3)` one has

```text
F(a,3)=-120a+240,
F(2,3)=0,       F_a(2,3)=-120.                        (4)
```

Thus the second root in the `a` direction is infinity, namely the
first retained row. The five integral gradient multipliers at these
old rows are

```text
(-222,34,9,-120,77),
```

so the old point is unramified in every row, with every displayed
multiplier a unit modulo `29`.

The discriminant of (3) in `a` is

```text
361b^4-1572b^3+4494b^2-4788b+1521.                    (5)
```

It is squarefree over `Q` and modulo `29`, as checked by exact
polynomial gcd with its derivative in both fields. The checker
derives (5) as `B^2-4AC` from the actual coefficient polynomials
`A=81-27b`, `B=19b^2-84b-39`, and `C=15b^2+35b` in (3).
The three coefficients of (3) have polynomial gcd one:
the only root of `81-27b` is `3`, and the middle coefficient is
`-120` there. Gauss's lemma and the nonsquare discriminant therefore
prove geometric irreducibility of the chart curve.

In the fixed Petersen edge ordering, the fifteen boundary-node
coefficients are

```text
(8,-4,-4,-35,39,-4,27,-81,54,-19,-15,34,8,-50,-42).
```

Every entry is nonzero and is a unit modulo `29`. Hence the section
has no boundary component and avoids every boundary node, including
after reduction modulo `29`. The usual open quotient and boundary
factor argument extends the chart irreducibility to the original
invariant polynomial. In particular, neither reducibility nor a
boundary-node degeneration explains (1).

## 2. Literal conjugate-primitive Gaussian rows

Put `pi=5+2i`, `n=29^m`, and `H=pi^m`. Choose an even integer `r`
with

```text
r == i mod H,
0<=r<58n,
29 does not divide Norm((-r-2nt+i)/H), t=0,1,2,3.     (6)
```

These choices always exist. If `H=A+iB`, then `B` is a unit modulo
`n` and `r=-A/B mod n` gives the first congruence. Among fifty-eight
successive representatives `r_0+kn`, twenty-nine are even. For each of
the four values of `t`, exactly one class among those twenty-nine gives
an additional factor `pi`; their union excludes at most four classes.
The conjugate factor cannot occur because `Y=1` and `r=i mod pi`.

Take the five rows

```text
P_0=-1,
P_1=-r+i,
P_2=-r-2n+i,
P_3=-r-4n+i,
P_4=-r-6n+i.                                         (7)
```

The even real parts and imaginary part one make the four nonconstant
rows conjugate-primitive, including at `1+i`. They are distinct, and
the constant row has a different direction. Let

```text
U={1,2,3,4},
H_U=H,        H_V=1 for every other nonempty cut V,
K_0=-1,       K_j=P_j/H for j=1,2,3,4.                (8)
```

The only nonunit core has odd split support, and every `Norm(K_j)`
is a `29`-adic unit by (6). All supports and corrections here are
literal Gaussian integers.

The old normalized rows used in Section 1 are transformed by

```text
M=[ -r  -1 ] [ 1  0  ],       det M=2n,
  [  1   0 ] [ 0  2n ]
```

and the first row is divided by `2n` to make it primitive. Relative
`GL_2` invariance and row homogeneity preserve the numerical zero of
the same fixed `Q`. They also preserve nonvanishing of its five
projective row derivatives.

For two finite labels corresponding to `s,t in {0,1,2,3}`,

```text
|Delta_st|=2n|s-t|,
Norm(gcd_G(P_s,P_t))=n,
|t_st|=b_st=2|s-t|<=6.                               (9)
```

The pair gcd in (9) has no extra prime: its extra norm would divide
`2|s-t|`, whereas the rows are odd and conjugate-primitive, so any
extra Gaussian gcd prime would have an odd split norm at least five.
Pairs with `P_0` have bracket and gcd norm one. Thus
`product b_jk=768`, independently of `m`, and every pair factor is a
`29`-adic unit.

## 3. Exact content and the lost frame factor

For the moving row `i=3`, direct transformation of (4) gives

```text
q_3(X,Y)=480n^2 Y(X+(r+4n)Y).
```

Consequently its ordinary coefficient content and old root are

```text
g_3=480n^2,
P_3=(-r-4n,1),       lambda_3=-480n^2.                (10)
```

The other primitive root is `(1,0)`, projectively the retained
`P_0`. This explains the large content without requiring ramification
at the old point. The clean old pair residues in (9) are unaffected.

At the only nonunit cut, row three is inside and the other-row count
is three. Thus the forced divisor of the
[Vieta content lemma](five_row_vieta_swap_content_and_residue_audit.md)
is `F_3=n`. Its excess is exactly

```text
E_3=g_3 Norm(K_3)/F_3=480n Norm(K_3),
v_29(E_3)=m.                                         (11)
```

On the other hand, normalizing the retained four rows returns `b=3`
in (3). The three quadratic coefficients then are
`0,-120,240`, whose gcd is a `29`-adic unit. A Bezout resultant for
those normalized coefficients therefore does not account for the
factor in (11). Pulling its identity back through the nonunimodular
frame costs precisely the kind of source-prime factor at issue.

Finally, (6)--(7) give the explicit bounds

```text
36n <= max_j Norm(K_j) <= 4097n,
16n <= Norm(K_3) <= 3845n,
7680n^2 <= E_3 <= 1845600n^2.                         (12)
```

For example, the largest absolute real part is `r+6n<64n`, and
the moving row has `r+4n<62n`; divide their squared moduli by `n`.
Thus `max|K_j|=sqrt(n) exp(O(1))` and `Norm(K_3)=n exp(O(1))`,
with absolute constants. This proves (2).
The example does not disprove a height estimate that charges the
full correction size. Establishing such an estimate in the full-fair
setting still requires control of approaches to the Vieta inverse
images of boundary points, in addition to the normalized coefficient
resultant. Boundary nonnodality alone does not supply that control.

## Verification

[check_five_row_vieta_clean_prime_content.py](check_five_row_vieta_clean_prime_content.py)
checks the fixed polynomial and its coefficient sum, derives the
discriminant from the chart coefficients, and verifies squarefreeness
over `Q` and modulo `29`. It checks that all fifteen boundary
coefficients, all five normalized old gradient multipliers, and the
normalized coefficient gcd `120` are `29`-adic units. For six
prime-power depths it constructs the
rows (7), verifies Gaussian primitivity and all ten actual pair
residues, and checks (10)--(12), including the exact ordinary content
and the growing correction size. The transformed old gradients remain
nonzero; the unit assertions refer to the normalized base point.
It asserts no full-fair family, correction-independent height bound,
or endpoint descent.

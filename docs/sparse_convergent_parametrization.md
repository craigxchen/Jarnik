# Exact continued-fraction parametrization of a sparse Gaussian window

## Status

This is an independent audit and refinement of the sparse completion
principle in
[sparse_phase_lattice_reconstruction.md](sparse_phase_lattice_reconstruction.md).
A fixed Gaussian
divisor `P` and bounds `U,T` with `2UT<|P|` determine at most one
conjugate-primitive multiple `h` up to sign in the window
`|h/P|<=U`, `0<|Im h|<=T`.

The refinement identifies that possible multiple exactly: it comes from
the last continued-fraction convergent below denominator `T` of a
specific square root of minus one divided by `N(P)`. Its denominator
is the small imaginary residue itself, rather than the substantially
larger remaining Gaussian multiplier. A large next-denominator gap
is necessary.

This supplies an explicit sparse parametrization. It does not bound the
number of possible fixed divisors, and does not prove the uniform
endpoint count.

## 1. The Gaussian ideal and its modular root

Let `P=a+ib` be a nonunit Gaussian integer coprime to its conjugate,
and put `Q=N(P)>1`. In the present odd split-prime setting, its
integer coordinates are coprime and its norm is odd. In particular

```text
gcd(a,Q)=gcd(b,Q)=1.
```

Choose the integer `r` in `[0,Q)` satisfying

```text
r=-b a^(-1) mod Q.
```

Then

```text
r^2=-1 mod Q,
P Z[i] = {x+it : x,t in Z and x=r t mod Q}.       (1)
```

For the second assertion, multiplication by `P` has integer matrix
`[[a,-b],[b,a]]`, of determinant `Q`. If
`x+it=P(u+iv)`, direct substitution gives the congruence in (1).
Both integer lattices have index `Q`: the first by this determinant,
the second because `(x,t) -> x-r t mod Q` is surjective. The
containment is therefore equality.

The modular root is unchanged when `P` is multiplied by a Gaussian
unit, since the Gaussian ideal is unchanged. There is no choice of
coordinate rotation or an associated factor of `sqrt(2)` in (1).

## 2. Uniqueness throughout the entire window

Assume `U,T>=1` and

```text
2UT<|P|.                                          (2)
```

If `h=x+it` and `h'=x'+it'` are multiples of `P` with

```text
|h/P|, |h'/P| <= U,
0<|t|, |t'| <= T,
```

their integer determinant is divisible by `Q`. Its absolute value
is bounded by

```text
|xt'-x't| <= 2|P|UT < Q.
```

Thus the determinant vanishes, and all lattice points in this
window lie on one real line. More generally the same argument
applies to the entire rectangle

```text
|x|<=|P|U,        |t|<=T.                         (3)
```

If both `h,h'` are coprime to their conjugates, their integer
coordinate pairs are primitive. Two collinear primitive integer
pairs differ only by sign. Hence there is at most one such
Gaussian integer up to sign across all allowed multipliers and
all nonzero imaginary residues in the window.

This is stronger than fixing one target integer and proving uniqueness
of its multiplier. It also permits the full remaining singleton block
to vary inside that multiplier.

## 3. The only possible continued-fraction index

Take the sign of a possible primitive `h=x+it` so that `t>0`.
By (1),

```text
k=(r t-x)/Q
```

is an integer. Also `gcd(k,t)=1`, since a common divisor would
divide both `x` and `t`. Equation (2) gives

```text
|r/Q-k/t|=|x|/(Qt)
 <= U/(|P|t)
 < 1/(2Tt)
 <= 1/(2t^2).                                   (4)
```

The usual elementary Legendre criterion therefore makes `k/t`
a regular continued-fraction convergent of `r/Q`.

Use the canonical finite regular expansion, with last partial
quotient at least two. Let its convergents be
`p_j/q_j`, beginning with `p_0/q_0=0/1`. Define

```text
j_* = the largest index j for which q_j<=T.       (5)
```

This index is well defined because `T>=1`. The rational endpoint
has denominator `Q`, which exceeds `T`: (2) and `U>=1` give
`T<|P|/2<Q`. Thus `j_*` is a proper convergent and has a next
convergent. If the initial convergents have repeated denominator
one, the largest-index convention in (5) resolves that case.

The possible primitive `h` necessarily corresponds to **exactly**
the index `j_*`.

To prove this, suppose it corresponds to an earlier index `j`.
The next convergent has `q_(j+1)<=T`. The positive Euclidean
remainders

```text
R_j=|r q_j-Q p_j|
```

strictly decrease along the proper convergents. Hence the next
lattice vector

```text
h_next=(r q_(j+1)-Q p_(j+1))+i q_(j+1)
```

also lies in rectangle (3). Its determinant with `h` is exactly
`+Q` or `-Q`, by the consecutive-convergent determinant identity.
This contradicts the strict rectangle bound in Section 2.
The next vector need not satisfy the original circular modulus
bound `|h_next/P|<=U`; its membership in the rectangle is all
the argument uses.

Consequently the decision procedure is exact. Compute the single
convergent `p_(j_*)/q_(j_*)`, put

```text
t_*=q_(j_*),
x_*=r t_*-Q p_(j_*),
h_*=x_*+i t_*.
```

There is an admissible conjugate-primitive multiple in the window
if and only if `h_*` is coprime to its conjugate and

```text
x_*^2+t_*^2 <= Q U^2.                            (6)
```

If it exists, it is `+h_*` or `-h_*`. Integrality of `h_*/P`
is already guaranteed by (1). No search over target integers,
small coefficients, or earlier convergents is needed.

The primitivity check itself requires no Gaussian factorization.
Because `gcd(p_(j_*),t_*)=1`,

```text
gcd(x_*,t_*)=gcd(Q,t_*).
```

Thus `h_*` is coprime to its conjugate exactly when
`gcd(Q,t_*)=1` and `x_*+t_*` is odd. The norm of its Gaussian
multiplier is the integer `(x_*^2+t_*^2)/Q`.

## 4. Necessary denominator and partial-quotient gaps

For the candidate, write `t=t_*`, `R_j=|x_*|`, and let
`t_next=q_(j_*+1)`, `R_next=R_(j_*+1)`. The continued-fraction
remainder identity is

```text
Q=|x_*| t_next+R_next t,
0<=R_next<|x_*|.
```

Therefore

```text
Q/|x_*|-t < t_next <= Q/|x_*|,
t_next > |P|/U-T.                               (7)
```

The exact next rational endpoint is allowed; then `R_next=0`
and the upper bound is attained. The nonzero `x_*` needed here
is automatic: a conjugate-primitive purely imaginary Gaussian
integer is a unit and cannot be a multiple of the nonunit `P`.

Let `a_next` be the next partial quotient. Since
`t_next=a_next t+q_(j_*-1)`, with
`0<=q_(j_*-1)<=t`, we obtain

```text
a_next > |P|/(U T)-2.                            (8)
```

There is also the exact necessary coprimality

```text
gcd(t,Q)=1.                                     (9)
```

Indeed, if a rational prime dividing `Q` divided `t`, its
oriented Gaussian divisor in `P` would divide both `h` and
`it`, hence would divide the rational integer `x`. That rational
prime would then divide both coordinates of the primitive `h`.

Equations (7)--(9) give an explicit continued-fraction denominator
restriction. The index itself is `O(log(2+T))` by the usual
Fibonacci growth of convergent denominators. This index bound
can still grow with the radius through `T`; it is not a uniform
count of possible shared products.

## 5. Audit of the shared-block and nonprivate-block consequences

In the fixed-anchor system, put

```text
P_i=product_(S contains i, |S|>=2) H_S,
U_i=K_i H_{ {i} },
h_i=P_i U_i.
```

Fixing all nonsingleton blocks therefore fixes `P_i`. Whenever
common bounds on `|U_i|` and `|Im h_i|` satisfy (2), the entire
primitive `h_i`, not only a correction-to-target ratio, is fixed
up to sign. The continued-fraction formula above makes that
conditional parametrization explicit.

The symmetry investigation further fixes only the nonprivate cuts,
meaning cuts whose two sides both have at least two vertices.
For `k` selected vertices, write `B=2^(k-1)`. Each pair is
separated by `B/2` cuts, exactly two of which are private cuts.
If `D_i` denotes the block oriented toward its unique vertex,
then, with `h_ij/bar(h_ij)=z_j/z_i`,

```text
h_ij=K_ij D_j bar(D_i) Q_ij,
```

where `Q_ij` is fixed by the nonprivate blocks and contains
`B/2-2` such separating factors. The orientation convention is
important: reversing the definition of `h_ij` conjugates the
displayed private product.

If the block log norms are approximately `w`, the private
multiplier has log modulus approximately `w`, while

```text
log|Q_ij| approximately (B/4-1)w.
```

With negligible coefficient and residue budgets, (2) therefore
holds for all pairs once `B>8`, that is, `k>=5`. The relative
block errors may be retained explicitly; they tend to zero in
the extracted system and do not alter this strict inequality.
Thus fixed nonprivate blocks determine every primitive pair
numerator up to sign across the entire admissible private-block
and correction window.

Consequently they determine all labeled relative phases. If two
primitive Gaussian tuples have those same relative phases, then
`z'_i=lambda z_i` for every `i`, with `lambda in Q(i)`. The
common Gaussian gcd of the first tuple being a unit forces
`lambda` to be a Gaussian integer; primitivity of the second
forces it to be a unit. The two tuples therefore differ only by
a global Gaussian unit. For selected tuples having a common
divisor, first pass to their intrinsic primitive normalization;
the pair phase numerators are unchanged.

More generally, the same argument applies whenever the fixed
Gaussian factors separating each required pair have log modulus
larger than the unfixed factor budget plus the correction and
residue budgets. The conditional uniqueness is a statement about
sufficient shared divisor mass, not a property restricted to the
private-cut labels.

## 6. Scope and finite verification

Across rows, shared Gaussian divisors make the corresponding
modular roots agree on their common norm moduli. Combining the
resulting candidate lattice vectors recovers the normalized
determinant identities already recorded in
[all_anchor_primitive_transition.md](all_anchor_primitive_transition.md).
No new cross-row height inequality follows here from the
continued-fraction index or the next-denominator gap alone.

An independent read of the companion reconstruction note checked its
shared-window quantifiers, connected-graph normalization, private-cut
orientations, and targeted second/fourth-moment extraction. In particular,
central-layer removal costs `2s` in the signed frozen-minus-unfrozen
weight and the primitive pair correction costs another `2s`, giving
the stated `3 zeta W/4+2s+log 2` sufficient slack after the residue
bound is included. Its 49-block anchor and 63-block all-pair incidence
counts also agree with the displayed constructions.

Exact integer computations checked 2,072 windows, using primitive
Gaussian divisors with both coordinates between minus twelve and
twelve, and integer bounds `1<=U,T<=5` satisfying (2). There
were 600 nonempty primitive windows. Direct enumeration agreed
with the single-convergent decision rule in every case; the
remainder identity and denominator coprimality also held. These
checks supplement the proof and are not a Lean formalization.

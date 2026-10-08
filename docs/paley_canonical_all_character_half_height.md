# Every saturated character is above half-height in canonical Paley profiles

Let `q=3 mod 4` be prime, `M=q+1`, and `H` the normalized Paley
Hadamard matrix with one infinity row and column. Take `b>=5` physical
copies of each nonconstant column. At each finite row `i`, flip one
distinct physical entry of label `i`; at the infinity row, flip a
second physical copy of a fixed finite label `a_*`. The following
exact theorem holds for the full saturated integer character lattice:

```text
min_(c in Z^M, sum c=0, c!=0) B(c) = bM/2-2,          (1)
B(c)=sum_(physical j) |c^T s_j|/2.
```

Since there are `r=b(M-1)` physical columns, the minimum exceeds
`r/2` by `b/2-2>0`. With actual distinct nearby split rational primes,
every nonzero primitive Gaussian character therefore has log norm
strictly above half the total prime log weight. This is an **unbounded
all-character obstruction to using the individual primitive half-height
inequalities alone** at
the archived `W=O(M log M)` radius scale. The actual source arguments
are uncontrolled; no endpoint realization, exclusion by simultaneous
phases, or growth improvement for actual endpoint arcs is asserted.
The separate [prime-field polynomial theorem](paley_primefield_polynomial_character_gap.md)
obtains an all-character half-height obstruction for **every**
capacity-two flip assignment at sufficiently large prime Paley orders.
The theorem here needs the diagonal assignment but covers every Paley
order and gives the exact minimum.

The subsequent [compatible-phase relaxation](paley_compatible_phase_coordinate_gap.md)
uses the support bounds below to satisfy every scalar coordinate gap
with exact compatible endpoint phases at the same `O(M log M)` scale.
Its factors are complex numbers with prime squared moduli, not Gaussian
integers. Joint Gaussian integrality remains an additional requirement.

## Skew-Hadamard reduction

Use the signs

```text
H_(infinity,j)=H_(i,infinity)=1,
H_(i,i)=-1,
H_(i,j)=-chi(i-j)  (i!=j),
```

where `chi` is the quadratic character of `F_q`. Since `chi(-u)=-chi(u)`,
the row-switched matrix `K=D H`, with `D=diag(1,-1,...,-1)`, satisfies

```text
K=I+J,       J^T=-J,       J_ij in {+1,-1} for i!=j.
```

Only this skew-Hadamard identity and orthogonality are used in the
character proof. The same exact minimum holds for any normalized
skew-Hadamard matrix with this diagonal row-to-label assignment and
one extra assignment at the distinguished row.

For an integer zero-sum row character `c`, put `d=Dc` and
`T=H^T c=K^T d=d-Jd`. The constant-column coefficient is `T_infinity=0`.
Hadamard orthogonality gives `sum_a T_a^2=M sum_i c_i^2`.

First take `c` ternary, so `h=|supp(c)|` is even and
`sum_i c_i^2=h`. Put `t_a=|T_a|` and
`L=sum_a t_a/M=1+delta`; Hadamard inversion gives `delta>=0`.
Every nonzero `t_a` is even and lies between `2` and `h`. Parseval
then gives the exact defect identity

```text
sum_a t_a(h-t_a)=Mh delta.                          (2)
```

Let `k` count the supported row flips that lower a physical character
coefficient by one. Every other supported flip raises its absolute
coefficient by one, so

```text
B(c)=bM L/2+h-2k,
B(c)-r/2=b/2+h-2k+(bM/2)delta.                   (3)
```

At a finite row `i`, the assigned physical column has old sign
`H_(i,i)=-1`, so the flip is beneficial precisely when
`d_i T_i>0`. On the support of `d`,

```text
d_i T_i=1-d_i(Jd)_i.
```

The sum `d_i(Jd)_i` has `h-1` terms, each `+1` or `-1`. A beneficial
finite flip with `|T_i|=h` requires
`d_i(Jd)_i=-(h-1)`: every other supported vertex points toward `i`
in the skew sign tournament. At most one vertex can have this property,
because the signs on the edge between any two proposed vertices are
opposite. The infinity flip supplies at most one further beneficial
entry. Thus at least `(k-2)_+` beneficial finite flips have
`2<=|T_i|<=h-2`. Their labels `i` are distinct, and each contributes

```text
t_i(h-t_i) >= 2(h-2)
```

to (2). If `h>=4`, this proves

```text
delta >= 2(h-2)(k-2)_+/(Mh).                       (4)
```

For `k<=2`, (3) gives `B-r/2>=b/2+h-4>0`. For `k>=2`, substitute
(4) into (3):

```text
B(c)-r/2
 >= b/2+h-4+(k-2)[b(1-2/h)-2].                 (5)
```

The bracket is at least `1/2` when `b>=5,h>=4`, so
`B(c)>=bM/2` for every ternary character with `h>=4`.
If `h=2`, at most one finite supported vertex can be beneficial,
and the infinity flip adds at most one. Thus `k<=2`, `delta=0`, and
(3) gives `B(c)>=bM/2-2`.

If `A=max_i|c_i|>=2`, Hadamard inversion gives `L(c)>=A`, while the
`M` row flips change the unflipped coefficient sum by at most
`sum_i|c_i|<=MA`. Therefore

```text
B(c)>=M A(b/2-1)>=2M(b/2-1)>=bM/2-2.            (6)
```

Together these cases prove the lower bound in (1). It is attained by
`c=e_infinity-e_(a_*)`: the two supported flips both lower their
coefficients, so `B=bM/2-2`.

## Saturation and actual prime heights

Each finite label has at least `b-2>=3` unflipped physical copies.
Comparing a flipped copy with an unflipped copy forces twice each
rational row coefficient to be integral. Conversely every zero-sum
integer `c` gives the integral physical vector `c^T S/2`. Thus (1)
covers every nonzero character in the saturated rational-row lattice,
not merely a chosen short-support family.

Choose `r` distinct rational primes `p_j=1 mod 4` in a common
logarithmic interval `ell<=log p_j<ell+1/r`, as in the existing
[nearby-prime construction](strict_obtuse_prime_box_countermodels.md).
For each character, its conjugate-primitive Gaussian composite has
exact height `V=sum_j |c^T s_j| log p_j/2`, while
`W_0=sum_j log p_j`. From (1),

```text
V-W_0/2 > ell(b/2-2)-1/2 >0                     (7)
```

once `ell>1/(b-4)`. For every fixed `b`, such actual clusters exist
at `ell=4log M+O(1)`, so `W_0=O(M log M)`. The literal Gaussian
rows are primitive and distinct, but their angles need not lie on
an endpoint arc. Equation (7) says that no individual saturated
character can contradict a fixed-constant endpoint assumption via
the usual half-height grid gap. It does not constrain simultaneous
relations among all Gaussian phases.

The [exact checker](check_paley_canonical_all_character_half_height.py)
verifies the skew matrix, defect and flip identities, and the sharp
minimum on every zero-sum ternary vector at order twelve; it samples
additional characters at orders twenty and twenty-four. The proof
above covers all orders and arbitrary integer amplitudes.

# The tensor coefficient floor is sharp at a generic coalition

Fix `n,d>=2`, `m=nd`, and a source monomial whose coordinate-color
blocks are `B_1,...,B_n`, each of size `d`. For a coalition `T`, put

```text
f=#{a : B_a subset T},       e=#{a : B_a intersect T is empty}.
```

The tensor argument gives a local coefficient valuation at least
`min(f,e)` after the row-correction loss. This note constructs a
section of the actual CB-annihilator kernel attaining that order at a
generic one-parameter coalition. Thus the universal local bound on one
source coefficient cannot be strengthened. Special arithmetic residues
or simultaneous conditions at different primes may still force more.

## An explicit level-one kernel polynomial

Take `2n` labels with binary directions `(P_i,Y_i)` and vector variables
`v_i=(v_i^1,...,v_i^n)`. Define the `2n`-by-`2n` determinant

```text
L_(P,Y)(v)=det [ P_i v_i^a | Y_i v_i^a ]_(i=1..2n, a=1..n).
                                                               (1)
```

It is multilinear in the vector rows and has `GL_n` relative weight
`det^2`. For a source monomial assigning exactly two rows `{i_a,j_a}`
to each color `a`, its coefficient is

```text
+/- product_a (P_(i_a)Y_(j_a)-P_(j_a)Y_(i_a)).          (2)
```

Indeed, the coefficient determinant splits into the `n` two-by-two
blocks belonging to the colors, after row and column permutations.
For distinct directions these coefficients are nonzero.

In the affine chart `P_i!=0`, write `z_i=Y_i/P_i`. Up to the scalar
`product_i P_i`, (1) is `det[v_i | z_i v_i]`. For any distinct colors
`a,b`, let `D_(z,ab)=sum_i z_i v_i^a partial_(v_i^b)`. Since the
polynomial has total color-`b` degree two,

```text
D_(z,ab)^2 L_z=2! L_z(v_i^b=z_i v_i^a)=0.             (3)
```

The substitution makes column `b` of the first block equal column
`a` of the second block. This proves actual kernel membership directly;
no inference from a fusion dimension is needed.

## The local normalization and the chosen coefficient

Work over the formal DVR `K[[t]]`, where `K` has characteristic zero
and contains independent generic row parameters. For `i in T`, set
`P_i=t u_i`, `Y_i=1`; give the outside rows distinct generic constant
directions with `P_j!=0`. For the selected `2n` labels, let `s` be the
number in `T`. Each bracket has order one if both its labels belong
to `T`, and order zero otherwise. Formula (2) therefore shows that
the common coefficient order of `L_(P,Y)` is exactly

```text
g=max(0,s-n).                                         (4)
```

Every perfect matching has at least `s-n` internal `T` pairs, and a
matching attaining this minimum exists. Any matching is a color-pair
assignment for a source coefficient. Consequently `t^(-g)L_(P,Y)`
has all coefficients in `K[[t]]` and at least one unit coefficient.
This is a local primitive frame; it is not an assertion that the
globally primitive polynomial tuple has no collision zero.

Now choose two labels from every `B_a`. If `B_a` meets both `T` and
its complement, choose one label on each side. A full color contributes
two labels in `T`, and an empty color contributes none. Thus

```text
s=2f+(n-f-e)=n+f-e.
```

For the chosen source coefficient, (2) has precisely `f` internal
`T` pairs. Its order after the common normalization (4) is exactly

```text
f-max(0,f-e)=min(f,e).                                (5)
```

## Extension from two labels per color to arbitrary degree

There remain `d-2` labels in each color block. Partition them into
`d-2` rainbow groups, each containing one label of every color, and
let `H(v)` be the product of their `n`-by-`n` vector determinants.
The coefficient specified by the original color assignment is `+/-1`
in `H`. All its rows are disjoint from the `2n` rows in (1). Set

```text
Q(v)=t^(-g)L_(P,Y)(v) H(v).                           (6)
```

This is an integral formal-DVR source polynomial of relative weight
`det^d`. It belongs to the full degree-`d` CB annihilator. In fact,
`H` has color-`b` degree `d-2`, so the Leibniz expansion of
`D_(z,ab)^d(LH)` has every term zero: terms with at least two
derivatives on `L` vanish by (3), and the others have more than
`d-2` derivatives on `H`. The computation takes place in `K((t))`,
where every `z_i` is defined, and proves the required generic-fiber
kernel identity. Its source coefficients are already in `K[[t]]`.

Because the row sets are disjoint, the selected source coefficient
in (6) is exactly its coefficient in `t^(-g)L`, up to sign. Its
valuation is (5). This completes the sharpness construction for
every occupancy pattern and every `n,d>=2`.

The construction is local at one coalition. It supplies neither a
single global integer section with all these valuations simultaneously
nor a circle-point configuration with the extracted core profile.
In particular, it does not rule out a stronger bound using simultaneous
prime conditions, additional special residue relations, or the small
real coordinates. It gives no improvement to the general growth bound.

The [exact checker](check_local_tensor_coefficient_floor_sharpness.py)
expands the determinant, verifies the repeated-column annihilator
coefficientwise, verifies the product construction, and checks the
local orders across a finite occupancy grid. The proof above is
independent of that finite range.

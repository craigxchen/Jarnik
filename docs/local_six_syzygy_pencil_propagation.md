# Propagation of twelve-row local syzygies along one direction

The localized six-row polar relations generate the terminal coherent
kernel for **every** configuration

```text
(b_0,b_1,...,b_11)=(t,1,2,...,11),   t in P^1 minus {1,...,11},
```

with `t=infinity` interpreted as binary row `P_0=(0,1)`, over any
characteristic-zero extension of `Q`. This upgrades the one
sample and generic statement of
[the twelve-row calculation](generic_terminal_coherent_syzygies.md)
along one explicit one-dimensional family. It does not establish
generation for arbitrary distinct twelve-tuples. The argument uses the
all-distinct rank theorem for the coherent map, the exact six-row
identities, and one finite modular rank certificate.

The later [determinantal generation theorem](distinct_configuration_determinantal_kernel_generation.md)
settles generation for every distinct twelve-tuple by a separate argument.
The pencil proof and its additional constant-generator rank conditions
are retained here with their original scope.

## 1. A degree-one kernel propagation lemma

Let `F` be a characteristic-zero field, `E(t)=E_0+tE_1` an affine
matrix pencil, and suppose `rank E(alpha)=rho` for every `alpha` in an
allowed set. Write

```text
C=ker E_0 intersection ker E_1,    dim C=c,
```

and let `g_j(t)=a_j+t b_j` be polynomial kernel vectors. Suppose their
`F(t)`-span is the full generic kernel, and the constant vectors among
the displayed family span `C`. Then their specializations span
`ker E(alpha)` at every allowed `alpha`.

Indeed, put `l=dim ker E(t)=dim(source)-rho`. Choose `l-c` of the
degree-at-most-one vectors that are independent modulo `C` over
`F(t)`. Their leading coefficients are independent modulo `C`:
otherwise some nonzero constant combination of the chosen vectors,
after subtracting `t` times a vector of `C`, would be a constant
polynomial kernel vector, hence would lie in `C`, a contradiction.
Now extend scalars to contain an allowed `alpha`. If a linear
combination of the chosen vectors at `alpha` lay in `C`, subtract that
constant vector. The resulting polynomial combination would equal
`(t-alpha)v` for a constant vector `v`. Since it belongs to the
polynomial kernel of `E(t)`, the integral-domain identity
`(t-alpha)E(t)v=0` gives `v in C`, contradicting independence of the
leading coefficients. Thus the selected vectors and a basis of `C`
remain independent at `alpha`; their number is `l`, so they fill the
specialized kernel.

The same argument covers algebraic `alpha`, because linear
independence over `F` persists after scalar extension. It needs
constant rank at the specialization; the local-generator rank alone
does not supply that hypothesis.

## 2. Application to the explicit twelve-row pencil

Use the integral 462-vector source basis and 4620 localized generators
`Q_I(P)G_j(I^c)` of the
[generic syzygy note](generic_terminal_coherent_syzygies.md). Keep
`b_1,...,b_11=1,...,11` and set `b_0=t`. The coherent matrix `E(t)` is
affine in `t`. The [Nagata balancing theorem](coherent_nagata_balancing.md)
gives `rank E(t)=121` for every `t` distinct from `1,...,11`, hence
the kernel dimension is always `341` there.

Each localized generator has degree at most one in `t`: its matching
coefficient on a six-set `I` is independent of `t` if `0 not in I`,
and is linear if `0 in I`. The `binom(11,6)*5=2310` generators with
`0 not in I` are constant vectors in `C`. The exact modular
[checker](check_local_six_syzygy_pencil_propagation.py) certifies rank
`286` for these 2310 vectors over `F_101` at the fixed eleven
directions. Therefore their rank over `Q` is at least `286`.

The stacked pencil has rank `176` over `Q`, by the exact upper bound
and modular lower bound in the generic syzygy note (and independently
by [the all-size determinant-line proof](terminal_kernel_determinant_line_height.md)).
Consequently `dim_Q C=462-176=286`. Since all 2310 displayed
constant relations lie in `C`, they span `C` exactly. The earlier
certificate gives rank `341` for all 4620 local relations at `t=0`,
so they span the generic kernel over `Q(t)`. Section 1 applies and
proves

```text
ker E_(t,1,...,11)
  = span{Q_I(t,1,...,11) G_j(I^c): |I|=6, 1<=j<=5}
```

for every `t not in {1,...,11}` in characteristic zero.

The same proof includes `t=infinity`. Its selected degree-one
generators have leading coefficients independent modulo `C`, as shown
in Section 1. For a six-set containing label `0`, the corresponding
homogeneous matching coefficient is linear in `P_0=(X_0,Y_0)`, so its
value at `P_0=(0,1)` is exactly the leading coefficient in the affine
chart `P_0=(1,t)`. For a six-set avoiding `0`, the relation is one of
the constant vectors in `C`. These 341 specialized vectors are
independent, and all-distinct coherent rank is 121 also at infinity
after changing the common projective chart.

The certificate consists only of exact finite-field ranks. Its
characteristic-zero conclusions use dimension upper bounds and the
polynomial argument above; it is not a numerical extrapolation from
sampled `t` values.

## 3. Scope of the pencil argument

For arbitrary fixed distinct `b_1,...,b_11`, the
[all-size pencil theorem](terminal_kernel_determinant_line_height.md)
still gives `dim C=286`, and the coherent map still has rank `121` at
every allowed `t`. Thus the same propagation proof would work if the
2310 six-set relations avoiding label `0` spanned `C`, and if the full
local family spanned at one parameter on that fiber. Those two
localized-family rank conditions are proved here for the explicit
eleven-tuple and on a Zariski-open set of eleven-tuples, but not for
every distinct eleven-tuple.

The Schubert formulation describes the step this pencil method leaves
unproved. In the
degree-four Pluecker coordinate ring of `Gr(3,12)`, the coherent map
restricts to the linear Schubert space of 3-planes containing
`span(1,b)`. Its ideal is generated by linear incidence equations in
the full coordinate ring: for each four-set `A`, they are the two
equations `sum_(i in A) (-1)^i p_(A minus i)=0` and
`sum_(i in A) (-1)^i b_i p_(A minus i)=0`, with signs determined by the
order of `A`. The row-squarefree weight projection is not
preserved by a change of basis adapted to `span(1,b)`, so linear
generation there does not by itself prove that the squarefree kernel
is generated by six-row polar relations. Likewise, pointwise
generation after inverting pair brackets would give an algebraic
denominator statement, not a coefficient-height floor for arbitrary
integer sums of local generators.

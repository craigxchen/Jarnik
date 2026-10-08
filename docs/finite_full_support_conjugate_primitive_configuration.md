# A finite full-support conjugate-primitive three-row configuration

This note gives a finite oriented Gaussian configuration with three rows
`P_i=X_i+i`, all `X_i` even. Thus each row is primitive in `Z[i]` and is
coprime to its conjugate. The example has all seven nonempty Boolean core
blocks nonunit, pairwise coprime, and split-supported. It is a finite
configuration only: its block logarithms are not shown to be `w+o(w)`, and
it is not shown to arise from the original short-circle extraction.

## Exact rows and blocks

Take

```text
X_1=1,469,178,    X_2=31,686,    X_3=153,548,
Y_1=Y_2=Y_3=1.
```

These can be obtained from the second three-row fixture in
[three_row_discriminant_interpolation_scope.md](three_row_discriminant_interpolation_scope.md),
whose real coordinates were `(14573,31686,33103)`, by adding the row-specific
shears

```text
delta_1 = 109*157*85 = 1,454,605,
delta_2 = 0,
delta_3 = 109*85*13 = 120,445.
```

Each shear is divisible by every common or pair block incident to that row.
Consequently those oriented common and pair divisibilities persist. The old
private factors `73,4513,4549` do not all persist: the new private factors are
computed directly below.

The row norms factor as

```text
N_1 = X_1^2+1 = 5*17*89*109*157*16673
    = 109 * 157 * 85 * 1,483,897,
N_2 = X_2^2+1 = 13*109*157*4513
    = 109 * 157 * 13 * 4,513,
N_3 = X_3^2+1 = 5*13*17*61*109*3209
    = 109 * 85 * 13 * 195,749.
```

The seven Boolean core norms, in subset order, are

```text
T={1,2,3}: 109                 T={1,2}: 157
T={2,3}:   13                  T={1,3}: 85
T={1}:     1,483,897           T={2}:   4,513
T={3}:     195,749.
```

They are pairwise coprime. Their prime factors are all split primes:
`109,157,13,5,17,89,16673,4513,61,3209`, each congruent to `1 mod 4`.
The factors `109,157,13,89,16673,4513,61,3209` are prime; `85=5*17`.
The checker uses trial division through each square root.

Actual oriented Gaussian divisors are

```text
common {1,2,3}:  10+3i       norm 109
pair {1,2}:     -11+6i       norm 157
pair {2,3}:       3-2i       norm 13
pair {1,3}:      -9+2i       norm 85
private {1}:  1109+504i      norm 1,483,897
private {2}:    -47-48i      norm 4,513
private {3}:   -385-218i     norm 195,749.
```

Exact Gaussian division gives the private quotients displayed above after
dividing each row by its incident common and pair divisors. Each divisor is
primitive, and its norm is coprime to every other block norm. Since all row
real coordinates are even and `Y_i=1`, the rows have opposite coordinate
parity and are coprime to their conjugates. In particular, all three row
norms are odd; the ramified prime 2 occurs only in the pair determinants.

The exact pair norm gcds and determinant quotients are

```text
gcd(N_1,N_2)=109*157=17,113,    D_12=X_1-X_2= 1,437,492,
gcd(N_2,N_3)=109*13 = 1,417,    D_23=X_2-X_3=  -121,862,
gcd(N_1,N_3)=109*85 = 9,265,    D_13=X_1-X_3= 1,315,630,
t_12=D_12/gcd(N_1,N_2)=84,    t_23=-86,   t_13=142.
```

Thus no pair determinant vanishes. The even bounded quotients are permitted;
one must not impose odd pair gaps when the row norms are odd.

## Exact extension conditions and the four-row obstruction

For an oriented Gaussian block `H=a+bi` of norm `n`, with `gcd(a,n)=1`,

```text
H | (X+i)    iff    X = -b*a^(-1) (mod n).
```

Conjugating `H` changes the residue to `+b*a^(-1)`. For a proposed fourth
row, factor each old block norm into coprime split-supported pieces
`n_T=n_T^- n_T^+`, where `n_T^-` remains on the old subset `T` and `n_T^+`
is assigned to `T union {4}`. Full Boolean support requires both pieces to
be greater than 1 for every one of the seven old subsets. The oriented
divisibility requirements on row 4 then give one root congruence modulo each
selected prime power in the `n_T^+` pieces. As the old block supports are
pairwise coprime, the Chinese remainder theorem combines them into one
residue class modulo `Q=product_T n_T^+`. A full extension additionally
requires an integer in that class whose norm `X_4^2+1` has exactly the
prescribed split-prime allocation, with no accidental common block; CRT
alone does not ensure this final norm-factorization condition.

This particular configuration has no full 15-block refinement. Its old
common block `109` and pair blocks `157` and `13` are prime. Each would have
to contribute a nontrivial factor both to `T` and `T union {4}`, which is
impossible with its single rational prime. This is a fixed-fixture
obstruction only; it says nothing about a different three-row fixture whose
seven old core norms each have at least two split-prime components.

The exact checks are in
[check_finite_full_support_conjugate_primitive_configuration.py](check_finite_full_support_conjugate_primitive_configuration.py).

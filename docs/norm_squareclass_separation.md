# Six-wise separation of the primitive norm squareclasses

## Status

The centrally truncated reduction gives a new finite restriction on
the **actual primitive anchor norms**. In an extracted family of
near-uniform endpoint rows, no nonempty product of at most six such
norms can be a rational square, once the original cluster is large
enough. Thus their squareclasses have the separation property of
columns of a binary parity-check matrix of minimum distance at least
seven.

This does not prove the uniform endpoint bound. The resulting lower
bound on the correction conductor is logarithmic in the number of
selected rows and is compatible with its present upper bound. The
argument uses a square-root separation estimate as well as the
full cut profile; it is not just the pairwise exclusion of a common
squareclass.

## 1. Setup and precise finite conditions

Use the notation and common-unit selection of
[endpoint_central_truncation.md](endpoint_central_truncation.md).
The selected original points are `z_0,...,z_(k-1)`, with

```text
z_i=epsilon product_p pi_p^a_ip bar(pi_p)^(e_p-a_ip),
W=log R^2,
|arg(z_i/z_0)|<=C exp(-W/4).
```

Choose the arguments in the inherited short interval about zero.
Write `s=log D`, `W_core=W-2s`, and `B=2^(k-1)`. The retained core
blocks have weights

```text
w'_S=log N(H'_S),
|w'_S-W_core/B|<=eta W_core/B,
eta=eta_core.
```

The primitive anchor numerators are

```text
h_i=product_p pi_p^((a_ip-a_0p)_+)
                  bar(pi_p)^((a_0p-a_ip)_+),
z_i/z_0=h_i/bar(h_i),
n_i=N(h_i),                 1<=i<=k-1.          (1)
```

We will prove that a squareclass relation with support size
`1<=q<=6` is impossible under the following explicit sufficient
conditions:

```text
min_S w'_S > 2q s,                              (2)
F_q=(1+eta)W_core m_q+2q s,
(q C/2) exp((F_q-W)/4) < 1,                    (3)
```

where

```text
q       1,2       3,4       5,6
m_q     1/2       3/4       15/16.              (4)
```

For the growing selection `k=floor(c log_2 M)`, any fixed
`0<c<1/2`, central truncation gives

```text
eta -> 0,
s/min_S w'_S -> 0,
W_core/W -> 1,
W -> infinity.
```

Thus (2)--(3) hold simultaneously for every `q<=6`. The same
conclusion holds for a fixed sufficiently large number of selected
rows with these limiting properties.

## 2. A square norm forces separation from a real direction

Let `g` be a Gaussian integer coprime to its conjugate and suppose
`N(g)` is a positive rational square. Since every split-prime
exponent of `g` then has even parity, and no inert or ramified prime
can occur in this primitive situation,

```text
g=epsilon w^2,    epsilon in {1,-1,i,-i},
H=|g|=|w|^2.
```

If `Im g != 0`, then

```text
|Im g|/|g| >= H^(-1/2).                         (5)
```

Indeed, write `w=a+i b`. For `epsilon=+/-1`, a nonzero imaginary
part has magnitude `2|ab|>=sqrt(a^2+b^2)`. For `epsilon=+/-i`,
it has magnitude `|a^2-b^2|`. When this is nonzero,
`||a|-|b||` is a positive integer, so
`|a^2-b^2|>=|a|+|b|>=sqrt(a^2+b^2)`. This also covers the cases
where one coordinate is zero.

The unit check matters: a square norm alone does not choose the
unit multiplying `w^2`. The bound (5) is valid for all four units.
With the literal prime representatives used below the numerator
is in fact a square without an additional unit, but this stronger
observation is unnecessary.

## 3. A hypothetical short squareclass relation

Suppose a nonempty subset `J subset {1,...,k-1}`, of size `q<=6`,
has `product_(i in J) n_i` a rational square. Choose coefficients
`c_i in {+1,-1}` on `J` with the signs as balanced as possible:
`ceil(q/2)` positive and `floor(q/2)` negative coefficients. Set

```text
beta_p=sum_(i in J) c_i(a_ip-a_0p),
g=product_p pi_p^((beta_p)_+)
             bar(pi_p)^((-beta_p)_+).           (6)
```

This is a Gaussian integer coprime to its conjugate. Since the
original selected rows have one common unit, exactly

```text
g/bar(g)=product_(i in J) (z_i/z_0)^c_i.         (7)
```

Modulo two, `c_i=1` and `|a_ip-a_0p|=a_ip-a_0p`. The assumed
squareclass relation therefore makes every `beta_p` even.
Consequently `N(g)` is a square.

Write the central decomposition of the original exponents as

```text
a_ip=b'_ip+d_ip,     0<=d_ip<=2r_p,
s=sum_p r_p log p.
```

For a prime in core block `S`, the core contribution to `beta_p`
is, up to its common orientation relative to `pi_p`,

```text
e'_p L_S,         L_S=sum_(i in S intersect J) c_i.
```

Its correction has magnitude at most `2q r_p`. Hence

```text
log N(g)=sum_p |beta_p| log p
 <= sum_S w'_S |L_S|+2q s.                     (8)
```

Uniformly averaging `|L_S|` over the `B` subsets `S` amounts to
averaging the absolute value of a sum of `q` independent Bernoulli
bits with the chosen balanced signs. Its expectation is `m_q` in
(4). Thus (8) proves

```text
log N(g) <= (1+eta)W_core m_q+2q s=F_q,
log|g| <= F_q/2.                               (9)
```

For odd `q`, the coefficients do not sum to zero, so the anchor
cannot be omitted when working with the original oriented layer
bits. In that notation the linear combination is
`sum_i c_i(b_i-b_0)`. Its expectation agrees with (4): conditioning
on `b_0` either gives the displayed Bernoulli sum or its negative
after complementing all other bits. Equivalently, the odd support
sizes `1,3,5` have the same absolute expectation as the next even
sizes `2,4,6`.

## 4. Full cut support excludes an exact real relation

We must also rule out `g/bar(g)=1`; the angular bound alone would
not do so. Fix `j in J` and choose a core cut `S` with
`S intersect J={j}`. Such a block exists because every cut is
retained. Its `|L_S|` is one. Restricting the reverse triangle
inequality to its primes gives

```text
sum_(p in S) |beta_p| log p
 >= w'_S-2q sum_(p in S) r_p log p
 >= w'_S-2q s > 0                              (10)
```

by (2). Therefore `g` is nonunit. A real Gaussian integer
coprime to its conjugate must be `+/-1`; hence `Im g != 0`.
This is where independent positive cut blocks and the small
correction budget exclude exact cancellation. No unproved
multiplicative-independence assertion about arbitrary angles is
being used.

Let `theta_i` be the original small argument of `z_i/z_0` and put
`Theta=sum_(i in J) c_i theta_i`. Then (7) gives

```text
|Im g|/|g|=|sin(Theta/2)|
 <= (q C/2) exp(-W/4).                         (11)
```

On the other hand, (5) and (9) give

```text
|Im g|/|g| >= exp(-F_q/4).                     (12)
```

Combining (11)--(12) contradicts (3). This proves the claimed
six-wise separation. At the largest support size, the limiting
gap is explicit:

```text
F_6/W -> 15/16,
(W-F_6)/4=(1/64+o(1))W.
```

## 5. Binary-code consequence and the square-core case

Let `d` be the dimension over `F_2` of the subgroup of
`Q*/Q*^2` generated by `n_1,...,n_(k-1)`. Products over subsets
of at most three distinct indices have different squareclasses:
equality between two such products would give a nonempty relation
supported on their symmetric difference, of size at most six.
Therefore

```text
2^d >= sum_(j=0..min(3,k-1)) binom(k-1,j).       (13)
```

This is the usual sphere-packing count for a parity-check matrix
with no dependence of weight at most six. No coding-theoretic
existence theorem is needed.

The central corrections satisfy a useful exact divisor bound.
In `h_i=K_i A'_i`, each rational prime exponent of `N(K_i)` is
between zero and `2r_p`; hence

```text
N(K_i) divides D^2.                             (14)
```

If every core exponent `e'_p` is even, then all `H'_S` and `A'_i`
are Gaussian squares. In particular

```text
n_i = N(K_i)  in Q*/Q*^2,
d <= omega(D).
```

The same evenness holds whenever the original squared radius has
only even prime exponents, since `e'_p=e_p-2r_p` preserves parity.
Equations (13)--(14) give

```text
2^omega(D) >= sum_(j=0..min(3,k-1)) binom(k-1,j),
log D >= (log 5/log 2)
         log(sum_(j=0..min(3,k-1)) binom(k-1,j)). (15)
```

The second inequality uses that every prime dividing `D` is an
odd split prime, and hence at least five. For large `k`, it says

```text
omega(D) >= 3 log_2(k-1)-log_2(6)+o(1),
log D >= (3 log 5/log 2)log k-O(1).             (16)
```

All pair corrections have the analogue of (14), by the same
primewise argument with any two anchors. This common support does
not bound the number of possible squareclasses independently of
`D`.

## 6. Remaining gap

The new restriction does not yet force `log D` to be a positive
fraction of a core block weight. The lower bound (16) is compatible
with `log D<=kW/sqrt(M)` even though the latter is negligible
relative to every block. Binary spaces can also contain large
sets with no dependencies of weight at most six once their
dimension grows. Thus a rank contradiction would need an
additional upper bound on the norm-squareclass rank, and no such
bound is proved here.

The threshold six is intrinsic to this particular estimate: for
balanced support sizes seven and eight, the absolute Bernoulli
expectation is `35/32>1`. Then (9) and the elementary square-root
gap no longer contradict the original angular scale. Extending
the conclusion to those support sizes requires another argument,
not a change in the constant of the present proof.

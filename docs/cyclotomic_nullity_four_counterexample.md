# Counterexamples to the proposed weighted nullity bound of three

The proposed universal bound `dim_Q ker M<=3` is false, even for
squarefree, divisor-closed supports. This note gives exact counterexamples
and an infinite family. The support calculation alone does not construct
an endpoint fiber. A separate [five-row construction](cyclotomic_five_row_subendpoint_family.md)
checks the additional width and parity requirements and realizes actual
points. Neither result disproves a finite uniform bound for the one-class
model or for lattice points on circles.

For an odd positive integer `h` and a finite set `E` of odd positive
integers, use the matrix and charged support cost

```text
M[d,e]=mu(e/d) if d|e, and 0 otherwise,   d odd, d<h,
D(E)=2*[1 in E]+sum_(e in E, e>1) phi(e).
```

The charge at `1` is the minimum positive base width compatible with
even base-layer parity. Thus the examples meet the proposed minimum-width
budget with `W_1=2` and every other `W_e=1`. This is a statement about a
permitted support budget, not the existence of an endpoint row family
whose exact active widths are those numbers.

## 1. Divisor-closed supports have an exact nullity formula

Suppose `E` is divisor closed. A matrix row indexed by `d notin E` is
zero: a nonzero entry would require `d|e` for some `e in E`. The square
submatrix on the low indices `d,e in E`, `d,e<h`, ordered increasingly,
is upper triangular with diagonal entries `mu(1)=1`. Consequently

```text
rank_Q M = #{e in E:e<h},
dim_Q ker M = #{e in E:e>=h}.                         (1)
```

For every `n in E` with `n>=h`, there is also an explicit kernel vector

```text
b^(n)_e=1 if e|n, and 0 otherwise.
```

Indeed its row `d` is zero if `d` does not divide `n`, and otherwise is
`sum_(f|n/d) mu(f)`, which vanishes because `d<h<=n`. In the examples
below the four high orders are mutually nondividing, so their four
kernel vectors restrict to the identity matrix on the four high
coordinates.

## 2. A small exact counterexample and an infinite family

For every odd integer `k>=15`, set

```text
n_1=(k-2)(k+2)(k+4),
n_2=k^2(k+4),
n_3=k(k+2)^2,
n_4=k(k+2)(k+4),
h=n_1,
E={positive divisors of at least one n_j}.
```

These orders satisfy `h=n_1<n_2<n_3<n_4<3h`. Since they are odd,
every proper divisor of any order is at most `n_4/3<h`. Hence the
four displayed orders are exactly the high indices of `E`; (1) gives
nullity exactly four.

The identity `sum_(e|n) phi(e)=n`, followed by inclusion-exclusion,
computes the cost without an enumeration. The six pairwise gcds, in
lexicographic order, are

```text
k+4, k+2, (k+2)(k+4), k, k(k+4), k(k+2).
```

The four triple gcds are `1,k+4,k+2,k`, and the fourfold gcd is `1`.
These formulas follow from pairwise coprimality of `k,k+2,k+4` and
`gcd(k^2,k^2-4)=1` for odd `k`; no primality assumption is used.
Thus

```text
U=sum_(e in E) phi(e)=4k^3+15k^2-4k-24,
D(E)=U+1,
h=k^3+4k^2-4k-16,
D(E)-4h=-k^2+12k+41=77-(k-6)^2<0.                  (2)
```

At `k=15`, this gives

```text
h=4199,
(n_1,n_2,n_3,n_4)=(4199,4275,4335,4845),
D(E)=16792 < 4h=16796,
dim_Q ker M=4.
```

This is not claimed to be the smallest counterexample. The earlier
`k=21` instance has orders `(10925,11025,11109,12075)` and charged
cost `43552<43700=4h`.

## 3. The failure persists for squarefree supports

Let

```text
(p_1,p_2,p_3)=(181,193,199),
Q=p_1 p_2 p_3=6951667,
n_j=(p_j-2)Q/p_j for j=1,2,3,
n_4=Q.
```

The six numbers `179,181,191,193,197,199` are distinct primes.
Consequently all four orders, and all their divisors, are squarefree.
Explicitly,

```text
(n_1,n_2,n_3,n_4)=(6874853,6879629,6881801,6951667),
h=n_1.
```

Let `E` again be the union of their divisor sets. It has twenty
indices. Each proper divisor is below `h`, so its nullity is four.
For `i!=j<=3`, `gcd(n_i,n_j)=Q/(p_i p_j)`, which divides `Q`, while
`gcd(n_j,Q)=Q/p_j`. Inclusion-exclusion therefore simplifies to

```text
sum_(e in E) phi(e)
 =Q+sum_(j=1)^3 (n_j-Q/p_j)
 =4Q-3sum_(j=1)^3 Q/p_j.
```

The charged cost is

```text
D(E)=27478592 < 27499412=4h.                          (3)
```

In particular the desired inequality “a divisor downset with four
indices at least `h` has weight greater than `4h`” is itself false.

## 4. Squarefree compression is valid, but does not prove the conjecture

There is a valid rank reduction behind the attempted downset argument.
For squarefree indices, multiplying each matrix row and column by its
Möbius sign turns its nonzero entries into the incidence entries
`1_(d|e)`. Rows with nonsquarefree indices are zero and can be ignored.

Fix a prime `p` and compress the column support down in its Boolean
coordinate: replace an upper-only index `ps` by `s`; retain both if
both are present. Split the evaluation rows according to whether their
indices contain `p`. The lower and upper columns above the same
remaining subset have the forms `(u_s,0)` and `(u_s,v_s)`.

Both the old and compressed spans project onto the same space spanned
by all `u_s`. The old projection kernel contains every `(0,v_s)` from
a retained pair. The compressed span is exactly the direct sum of the
space spanned by all `(u_s,0)` and those paired vertical vectors.
Thus compression cannot increase rank. It preserves the number of
columns, so cannot decrease nullity.

It also cannot increase the charged cost: away from the base it removes
a factor `p-1>=2` from the totient; replacing `p` by `1` changes cost
`p-1` to `2`, again without an increase. Iterating decreases the
indices until a divisor downset is obtained. Formula (1) then applies.
The squarefree example in Section 3 shows exactly why this valid
reduction does not imply nullity at most three.

## 5. Consequences and verification scope

The exhaustive small-budget checks through odd `h=31` remain valid.
They do not extend to all odd thresholds. A claimed universal count
derived by substituting `rho<=3` into a nullity-to-cardinality bound
must therefore be withdrawn.

The four-dimensional kernel alone does not exhibit many endpoint
rows. Such rows must still be integer allocations inside one common
width box, satisfy even base parity and equal-norm phase alignment,
and meet the contact inequality for the **actual** total width.
These conditions do not follow from a minimum charged support cost.
The examples establish neither unbounded nullity under `D(E)<=4h`
nor an unbounded endpoint fiber.
The subsequent [uniform rank theorem](cyclotomic_uniform_rank.md)
proves nullity at most 1125 for every finite odd support under this
budget. A separate [five-row construction](cyclotomic_five_row_subendpoint_family.md)
uses explicit allocations, beyond the rank-four support alone, to
produce five actual subendpoint points.

There is now a positive replacement for the refuted constant three:
[squarefree compression and an intersection graph](cyclotomic_squarefree_uniform_rank.md)
prove `dim ker M<=1125` for all squarefree active supports within the
same budget, with no bound on their number of prime divisors. The
intermediate downset estimate also permits prime powers; extending the
compression to arbitrary nonsquarefree supports remains open.

The extended [exact checker](check_cyclotomic_weighted_nullity.py)
verifies the cost, four kernel vectors, divisor closure and exact
nullity of the counterexamples, along with a finite range of the
infinite family. It retains the separate exhaustive small-threshold
enumeration. The displayed formulas prove the infinite assertions.

# An optimal all-split `q=2` four-norm contact profile over `Q`

The contact degrees `14,14,14,8` admit a determinant-one frame on
`P^1/Q` with common pole degree **seven**. All four prescribed divisors
are defined over `Q`, reduced, pairwise disjoint, and supported at
closed points whose residue fields contain `i`. Seven is the minimum
possible pole degree for this contact problem.

This is an exact example for the four norm-contact conditions only.
It does not assert the additional source, cut, integer-specialization,
or endpoint-circle conditions of the lattice-point problem, and does
not settle its uniform bound `M(C)`.

## 1. The rational function and frame

On the affine line with coordinate `t`, put

```text
Q=t^7+1,                 S=t^4-t^3=t^3(t-1),
P=S-2Q=-2t^7+t^4-t^3-2,
f=P/Q=-2+S/Q,

U=((f,f-1),(1,1)).                                  (1)
```

The determinant is one. Since the only roots of `S` are `0,1`, and
`Q(0)=1`, `Q(1)=2`, we have `gcd(P,Q)=gcd(S,Q)=1`. The denominator
`Q` is squarefree, has degree seven, and `f` is regular at infinity
with value `-2`. Thus the least common pole divisor of the frame is
the degree-seven divisor `Z(Q)`:

```text
delta=7.                                             (2)
```

For column vectors `z,w`, the four required norms are

```text
F_1=||z||^2        =(P^2+Q^2)/Q^2,
F_2=||w||^2        =((P-Q)^2+Q^2)/Q^2,
F_3=||z+w||^2      =((2P-Q)^2+4Q^2)/Q^2,
F_4=||-3z+2w||^2   =((P+2Q)^2+Q^2)/Q^2.              (3)
```

Write `N_j` for the four displayed numerators. Each has degree
fourteen, with leading coefficients `5,10,29,1`, respectively, and
each is coprime to `Q`. In particular there are no norm zeros at
infinity or at the frame poles.

## 2. The contact divisors

The key factorization is

```text
N_4=S^2+Q^2
   =t^14+t^8+t^6+1
   =(t^8+1)(t^6+1).                                  (4)
```

Define effective divisors over `Q` by

```text
D_1=Z(N_1),   D_2=Z(N_2),   D_3=Z(N_3),
D_4=Z(t^8+1).                                        (5)
```

Their degrees are exactly `14,14,14,8`, and `D_j<=Z(F_j)`.

The three first numerators are squarefree over `Q`. A short exact
certificate is reduction modulo seven: each keeps degree fourteen,
and polynomial Euclidean division in `F_7[t]` gives

```text
gcd(N_j mod 7, N'_j mod 7)=1       for j=1,2,3.        (6)
```

The attached standard-library verifier reproduces (6). The fourth
numerator is squarefree directly from (4): its factors are squarefree
in characteristic zero and coprime, since a common root would satisfy
`t^2=1` and `t^8=-1`. Hence all four selected divisors are reduced.

The eight possible values of `f` at zeros of the four norms are

```text
{i,-i}, {1+i,1-i}, {1/2+i,1/2-i}, {-2+i,-2-i}.        (7)
```

These sets are pairwise disjoint. Since `Q` is invertible at every
norm zero, (7) proves pairwise disjointness of the four complete norm
zero divisors, and therefore of the selected divisors (5).

Each closed point in `D_1` has a residue-field element `f` squaring
to `-1`; the corresponding elements for `D_2,D_3,D_4` are
`f-1`, `f-1/2`, and `f+2`. Thus every contact is split. More explicitly,
`t^8+1` is the sixteenth cyclotomic polynomial, and `t^4` in its residue
field squares to `-1`.

## 3. Optimality and precise consequence

For any rational determinant-one frame with least common pole divisor
`Delta`, each norm function has pole divisor bounded by `2Delta`.
It is nonzero, so its zero divisor has degree at most `2deg Delta`.
The first degree-fourteen contact condition alone therefore forces

```text
14 <= deg Z(F_1) <= 2delta,       delta>=7.             (8)
```

Our frame attains this lower bound. Consequently the scaled
four-norm profile at `q=2` cannot exclude the threshold `delta<=4q=8`,
even with reducedness, disjointness, and descent to `Q` imposed. The
odd-degree parity obstruction for `7,7,7,4` does not extend to its
doubled profile using these conditions alone.

The construction uses a rational function with seven finite poles;
it is not an example within a polynomial-only frame ansatz. Additional
restrictions on pole support, or other conditions from the original
lattice problem, would need a separate argument.

## 4. A specific failure of the individual cut refinement

In this example each of `N_1,N_2,N_3` is irreducible over `Q`. Exact
certificates are furnished by irreducibility of the degree-fourteen
reductions modulo `43,127,59`, respectively. The verifier checks the
degree-fourteen case of Rabin's finite-field criterion: for the
corresponding polynomial `N` and prime `p`,

```text
t^(p^14) = t mod N,
gcd(N, t^(p^7)-t) = 1,
gcd(N, t^(p^2)-t) = 1,                               (9)
```

all in `F_p[t]`. The polynomial degrees are preserved on reduction.
Since `2,7` are the prime divisors of fourteen, these are a complete
irreducibility certificate; irreducibility of the reductions implies
irreducibility over `Q`.

Consequently each of `D_1,D_2,D_3` is a single degree-fourteen closed
point over `Q`. In particular it cannot be decomposed as a sum of seven
effective, `Q`-defined degree-two contact divisors. This frame therefore
fails that individual-cut refinement, despite satisfying its aggregate
degree requirements. This is a precise limitation of this construction,
not a claim that every actual lattice profile must have degree-two
cut divisors on a rational source.

## Reproduction

Run from the repository root:

```text
python3 docs/check_scaled_split_contact_profile_q2.py
```

The verifier uses integer polynomial arithmetic for the identities
and exact finite-field polynomial arithmetic for squarefreeness,
coprimality, and irreducibility certificates; it performs no numerical
root search.

An independent audit also computed the gcds directly over `Q[t]`,
checked all residue-field claims, and checked regularity at infinity.
The root orchestrator independently verified the factorization and
ran the standard-library certificate.

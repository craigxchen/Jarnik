# A denominator-aware split-prime capacity bound

This note supplies a local restriction involving both the common cotangent
denominator and the exact all-edge norm lcm. It does **not** prove the
fourth-power endpoint inequality.

Let `L>0`, let `X_1,...,X_k` be distinct integers, and assume every
`Q_ij=(X_i X_j+L^2)/(X_j-X_i)` is an integer. Include the anchor edge
cotangents `Q_0i=X_i`. For each edge reduce `Q_e/L=a_e/b_e`, put
`epsilon_e=2` when both reduced entries are odd and `1` otherwise, and set

```text
N = lcm_e ((a_e^2+b_e^2)/epsilon_e).
```

For every rational prime `p=1 (mod 4)`, write `e=v_p(L)` and `s=v_p(N)`.
Then

```text
k+1 <= (s+1)(p-1)p^e.                                    (1)
```

Consequently

```text
v_p(N) >= max(0, ceil((k+1)/((p-1)p^v_p(L)))-1).            (2)
```

The product of the prime powers on the right of (2), over split primes,
divides `N`. Only primes `p<=k+1` can contribute a positive exponent.

## Proof with the pair denominators retained

Choose a root `r` of `r^2=-1` in `Z_p`, and work in `Q_p`. Define

```text
q_0=1,                 q_i=(X_i+rL)/(X_i-rL).
```

These are nonzero and distinct. For an oriented edge, the corresponding
ratio is `(a+rb)/(a-rb)`, where `a,b` are its primitive cotangent
coordinates. Since `p` is odd and `gcd(a,b)=1`, the two factors `a+rb`
and `a-rb` cannot both have positive valuation. Both are integral, so

```text
v_p((a^2+b^2)/epsilon)
    = |v_p((a+rb)/(a-rb))|.
```

Writing `h_i=v_p(q_i)`, the exact all-edge lcm therefore gives

```text
s = max_(i,j)|h_i-h_j| = max_i h_i-min_i h_i.               (3)
```

There are at most `s+1` occupied integer valuation levels.

Fix a level `h`, and write its points as `q_i=p^h u_i` with
`u_i in Z_p^*`. For two points on this level, the half-angle cotangent
identity is

```text
Q_ij/L = r(u_i+u_j)/(u_i-u_j),                             (4)
```

up to an irrelevant sign depending on the edge orientation. If
`u_i=u_j (mod p^(e+1))`, then `u_i+u_j` is a unit, because `p` is odd.
Equation (4) would give `v_p(Q_ij/L)<=-e-1`, contradicting the integrality
of `Q_ij`. This also applies to edges involving the anchor.

Thus the units on a level occupy distinct classes modulo `p^(e+1)`.
There are `(p-1)p^e` such classes. Multiplication by the number of levels
in (3) proves (1).

## Sharpness, including positive cotangents

For the unrestricted integer clique

```text
L=6,       X=(-18,-6,-3,0,2,6,12),
```

all edge norms are `1` or `5`, and `N=5`. At `p=5`, equation (1) is
the equality `8=(1+1)(5-1)`.

There is also a positive example satisfying `X_1>=L`:

```text
L=12,
X=(264,288,294,300,304,312,324),
N=2409083501772215645.
```

Here every pair quotient is integral, `v_5(L)=0`, and `v_5(N)=1`, so the
same local equality holds. These slopes arise by adding `25` to each
slope of the preceding example. The displayed lcm and all pair
divisibilities were checked with exact integer/rational arithmetic.

One can retain `L=6` as well. For the original signed integer coordinates,
the lcm of all pair differences is `D=360`. If `T` is divisible by `D`,
translation `X_i -> X_i+T` preserves every pair divisibility: the
numerator changes by `T(X_i+X_j)+T^2` and its denominator is unchanged.
Taking `T=lcm(D,25)=1800` gives

```text
L=6,
X=(1782,1794,1797,1800,1802,1806,1812),
N=854184548609224587096599042100005.
```

Again `v_5(N)=1`. Indeed the original anchor Gaussian factors have
5-adic valuations at most one, and translation by a multiple of `25`
leaves both their valuations and their orientations unchanged. Equation
(3) then preserves the all-edge exponent. This also gives an infinite
family with all coordinates positive, by taking arbitrary positive
multiples of `1800`; its least cotangent tends to infinity while the
local data `v_5(L)=0`, `v_5(N)=1` remain fixed. The finite displayed
example was independently checked by taking the exact all-edge lcm.

## Scope and the remaining gap

The earlier inert-prime and dyadic estimates in
[integer_cotangent_local_bounds.md](integer_cotangent_local_bounds.md)
bound cardinality using `L` alone. At a split prime, unbounded valuation
levels prevent such an estimate; the exact lcm supplies their range and
gives (1). No matching split-prime capacity statement was found in the
cotangent/local/residue notes during this bounded audit.

For a fixed number of points, (2) only forces powers of finitely many
small split primes. It contains no archimedean parameter `X_1/L`.
Therefore it provides no fourth-power lower bound without an additional
argument linking the cotangent size to these valuation ranges, or to
contributions from larger split primes. That additional argument is
not supplied here.

The independent exact checker
[`check_split_prime_cotangent_capacity.py`](check_split_prime_cotangent_capacity.py)
verifies all three displayed norm lcms, the sharp local equalities, and
seven positive translations retaining the same five-adic norm exponent.

# Divisor estimates and the new near-real region

The general reduction in
[kurlberg_wigman_shrinking_arc_audit.md](kurlberg_wigman_shrinking_arc_audit.md)
produces at least `ceil(m(m-2)/4)` distinct Gaussian divisor classes of one
integer N, with

```text
(4/C^2)sqrt(N) <= Norm(A) <= N^(1/2+1/m),
1 <= |Im A| <= (C/2)N^(1/(2m)),       0<C<=2.           (1)
```

The classes come from pairs of an actual circle cluster. A uniform bound
for all conjugate-coprime divisors in (1) would suffice, but is a stronger
target: an arbitrary collection in this region need not have the pairwise
compatibilities of a circle cluster. This audit finds no applicable bound
of that strength. In particular, an algorithm for one residue class does
not supply it.

## 1. The residue-class algorithm cannot collect these divisors

Jonathon Hales, [*Divisors in Residue Classes Revisited*](https://arxiv.org/abs/2410.05030),
Theorem 1.1, extends Lenstra's algorithm to Gaussian integers and other
imaginary quadratic Euclidean rings. For coprime input data, it lists
divisors of an input B in one prescribed class modulo S using
`O(log |B|)` elementary operations when `|S|^3>|B|`. This is an algorithmic
statement, not a uniform count for a strip. The preprint is from 2024;
the publisher lists DOI
[10.1016/j.jnt.2026.03.014](https://doi.org/10.1016/j.jnt.2026.03.014)
in its January 2027 issue.

Here B would be the ordinary integer N, considered as a Gaussian integer.
Its complex absolute value is N, while its Gaussian norm is N squared.
Consequently the algorithm's threshold is

```text
|S| > N^(1/3).                                        (2)
```

There is a direct obstruction to using (2) to group the divisors in (1).
Their complex absolute values obey

```text
|A| <= X := N^(1/4+1/(2m)).                            (3)
```

For m>6 and

```text
N^((m-6)/(12m))>2,                                    (4)
```

equation (2) implies `|S|>2X`. If two such divisors A,A' were congruent
modulo S, their nonzero difference would be a Gaussian multiple of S,
and hence have absolute value at least |S|. But

```text
|A-A'| <= |A|+|A'| <= 2X < |S|.
```

Therefore A=A'. **Every divisor in this region already occupies a separate
residue class for every modulus satisfying the algorithmic threshold.**
Knowing how to list divisors in one such class gives no bound on how many
classes are occupied. The small imaginary part prescribes no common
large Gaussian modulus. Choosing a smaller modulus would leave the
hypothesis of the cited theorem.

For an actual cluster with m>=24 and C<=2, condition (4) follows
automatically. Summing the m-1 consecutive chord lower bounds gives
`C sqrt(R)>=(m-1)sqrt(2)`, hence `N=R^2>=(m-1)^4/4`.
At m=24, the exponent in (4) is 1/16 and
`23^4/4>2^16`. Both that exponent and the lower bound for N increase
with m, proving (4) for every m>=24. Thus there is no residual
large-radius assumption for the hypothetical unbounded clusters.

This argument does not exclude different applications of the algorithm
after an additional arithmetic transformation. Such an application would
have to construct and justify its new input, modulus, and shared residue.
Neither (1) nor the known pair construction currently provides them.

## 2. Every pair has a distinct ordinary norm

The injectivity holds for all `binom(m,2)` pairs, whenever `m>=3`, rather
than only for the near-critical pairs in (1). Retain the individual
angular bound before replacing it by a common height bound. For
`a=|Re A|`, `b=|Im A|`, and `d=Norm(A)`, it gives

```text
1 <= b^2 <= C^2 d/(4 sqrt(N)) <= d/sqrt(N) <= sqrt(d). (5)
```

The last two inequalities use `C<=2` and `d|N`, hence `d<=N`.
Since `d>1`,

```text
sqrt(d)-a = b^2/(sqrt(d)+sqrt(d-b^2)) < 1.             (6)
```

Thus `a=floor(sqrt(d))`, and `d` cannot be a square. The value of `d`
determines `b` up to sign, so it determines the full Gaussian
associate-and-conjugacy class. The B2 argument in Section 4 of the
linked reduction distinguishes the classes for every unordered pair;
it does not require the upper norm cutoff in (1). Consequently **all
`binom(m,2)` pair norms are distinct ordinary divisors of N**, satisfying

```text
d-floor(sqrt(d))^2 = b^2,
1 <= b^2 <= C^2 d/(4 sqrt(N)),
gcd(floor(sqrt(d)),b)=1.                               (7)
```

At least `ceil(m(m-2)/4)` of these distinct norms also lie in (1), and
their imaginary coordinates obey `b<=H=(C/2)N^(1/(2m))`.

The gcd in (7) follows from conjugate coprimality. Conversely, for odd N
whose prime factors split, every ordinary divisor d satisfying (7)
gives the conjugate-coprime Gaussian divisor `floor(sqrt(d))+ib` of N.
Indeed its coordinates are coprime and of opposite parity; at each
split prime only one conjugate Gaussian factor can occur, with exponent
no greater than its exponent in N. This converse concerns the displayed
arithmetic inequalities. The exact angular condition is
`b/sqrt(d)<=sin(Delta/2)`; replacing the sine by `Delta/2` gives a necessary
inequality, not an equivalent one. Neither converse says that an arbitrary
collection of such divisors comes from one circle cluster.

For completeness, the unrestricted region (1), which omits the individual
angular bound (5), also has norm injectivity for `m>=4`:
`H^2<=N^(1/m)<=N^(1/4)<=sqrt(d)`, so the same argument applies. In that
region its classes are exactly the divisors with `d-floor(sqrt(d))^2`
a positive square at most `H^2` and coprime coordinates.

The narrowness in (6) concerns representations of each individual d.
It does not put the near-critical distinct d in one short additive
interval: their allowed range is still `[sqrt(N),N^(1/2+1/m)]`.
Replacing that range by one of length `N^(1/4)` would be an unsupported
change of scale when m may grow arbitrarily slowly.

## 3. The available sector average has a different quantifier

O. V. Savastru, [*Divisor problem in special sets of Gaussian integers*](https://lib-repo.pnu.edu.ua/bitstream/123456789/1251/1/1088-3169-1-PB.pdf),
studies a summatory function of sector-restricted divisor counts over
Gaussian inputs in a progression. The displayed function on page 306 is

```text
T(x,gamma,omega_0,Sector)
  = sum_(omega=omega_0 mod gamma, Norm(omega)<=x) tau_Sector(omega).
```

Its asymptotic controls an average over omega. It does not bound the
divisor count for every individual ordinary integer N in the
N-dependent region (1). A potentially very sparse sequence of exceptional
N is sufficient to violate the endpoint goal, so an average cannot simply
replace that missing quantifier.

## 4. What an upper bound would have to accomplish

Define D(N,m,C) to count associate-and-conjugacy classes of nonreal
conjugate-coprime Gaussian A with `Norm(A)|N` satisfying (1). The proved
reduction gives

```text
ceil(m(m-2)/4) <= D(N,m,C)                              (8)
```

for every primitive fixed-unit m-point cluster. Thus, for example, any
proved estimate `D(N,m,C)=o(m^2)` uniform in N as m tends to infinity
would settle the goal after splitting into unit classes and subdividing
the original fixed-C arc into arcs of constant at most two.

A radius-dependent estimate such as `D<=N^epsilon`, with epsilon fixed,
does not do this for arbitrarily slow m growth. Nor does the statement
that the imaginary parts in (1) are `N^o(1)`: that quantity can still be
arbitrarily large. Any use of o-notation must retain the linked exponents
`1/m` and `1/(2m)` and all constants.

No divisor upper bound or new general point-count growth rate is proved
here. The exact conclusions are the large-modulus separation in Section 1
and the injective ordinary-norm reduction in Section 2. The former rules
out a direct invocation of the residue-class algorithm; the latter
expresses the remaining arithmetic count using ordinary integer divisors.

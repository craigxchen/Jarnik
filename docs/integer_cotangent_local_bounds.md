# Local bounds for integer cotangent cliques

Let \(L\ne0\) be an integer and let \(A\subset\mathbb Z\) be finite with
the pair condition

```text
(ab+L^2)/(a-b) is an integer                         (1)
```

for every distinct \(a,b\in A\).  This is the integer-coordinate form
of the cotangent condition.  Put \(e=v_p(L)\).  The following bounds are
purely local:

```text
|A| <= (p+1)p^e-1       for every p = 3 (mod 4),       (2)
|A| <= 4*2^e-1          for p = 2.                     (3)
```

The first bound is \(p^{e+1}+p^e-1\).  For \(L\) odd, (3) gives
\(|A|\le3\), recovering the sharp parity bound for \(L=1\).

## The divisibility valuation

For \(a\ne b\), (1) implies

```text
a-b | a^2+L^2.                                           (4)
```

Indeed \(a\equiv b\pmod {a-b}\).  If \(p=3\pmod4\), write
\(s=v_p(a)\).  Since \(-1\) is not a square in \(\mathbb F_p\),

```text
v_p(a^2+L^2) = 2s       if s<e,
                2e      if s>=e.                        (5)
```

The second line follows after factoring \(p^{2e}\): if \(s>e\), the
remaining sum is a unit, and if \(s=e\), it is a sum of two unit squares,
which is nonzero modulo \(p\).

## Proof for an inert odd prime

Partition \(A\) by \(s=v_p(a)\).  In a stratum \(s<e\), write
\(a=p^s x\), with \(x\) a \(p\)-adic unit.  For two members \(a,b\) in
this stratum, (4)--(5) give

```text
s+v_p(x-y) = v_p(a-b) <= 2s,
```

so \(x\) are pairwise distinct modulo \(p^{s+1}\).  There are exactly
\((p-1)p^s\) unit residue classes modulo \(p^{s+1}\).  Hence the stratum
has at most \((p-1)p^s\) members.

For the remaining stratum \(v_p(a)\ge e\), write \(a=p^e x\), allowing
\(x=0\).  The same argument gives

```text
e+v_p(x-y) <= 2e,
```

so the \(x\)'s are distinct modulo \(p^{e+1}\).  This stratum has at most
\(p^{e+1}\) members.  Summing the disjoint strata gives

```text
|A| <= sum_(s=0)^(e-1) (p-1)p^s + p^(e+1)
     = (p^e-1)+p^(e+1) = (p+1)p^e-1.
```

No cross-stratum restrictions were used, so this is a robust local bound;
the full integer condition can only lower the maximum.

## The dyadic bound

For \(p=2\), (5) is replaced by

```text
v_2(a^2+L^2) = 2s       if s<e,
                2e+1    if s=e,
                2e      if s>e.                           (6)
```

The first line is again separation of unequal valuations.  For \(s=e\),
the odd parts have squares congruent to \(1\bmod8\), so their sum has
valuation exactly one.  For \(s>e\), the odd part of \(L\) remains after
factoring \(2^{2e}\).

For \(s<e\), the normalized odd parts are distinct modulo \(2^{s+1}\),
giving \(2^s\) possibilities.  In the high stratum write \(a=2^e x\).
If \(x\) is even, (6) gives \(v_2(x-y)\le e\), hence at most \(2^e\)
even classes modulo \(2^{e+1}\).  If \(x\) is odd, it gives
\(v_2(x-y)\le e+1\), hence at most \(2^{e+1}\) odd classes modulo
\(2^{e+2}\).  Therefore

```text
|A| <= sum_(s=0)^(e-1) 2^s + 2^e + 2^(e+1)
     = (2^e-1)+2^e+2^(e+1) = 4*2^e-1.
```

The \(e=0\) case says there is at most one even member and at most two
odd members: for odd \(a,b\), \(a-b\mid a^2+L^2\equiv2\bmod4\), so every
odd pair differs by \(2\bmod4\).
For every odd \(L\), the set \(\{-L,0,L\}\) attains this bound.

## Scope and finite sharpness checks

These are \(L\)-dependent local cardinality bounds. They do not imply a
uniform bound as \(L\) varies, and they do not by themselves give a lower
bound for the all-edge primitive lcm \(N\) studied in
[integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md).

The local bounds are intentionally independent of the other primes and
need not be attained by integer cliques.  A divisor-based Bron--Kerbosch
search was run as follows: for a fixed minimum \(a\), candidates are
\(a+d\) with \(d\mid a^2+L^2\), and candidate pairs are then checked
against (1).  The search found the seven-element clique

```text
{-18,-6,-3,0,2,6,12}.
```

Every pair in this displayed set satisfies (1), and \(v_2(6)=1\), so it
attains the dyadic bound \(4\cdot2^1-1=7\), proving that the global
maximum at \(L=6\) is exactly seven. The inert-prime bound for
\(p=3\) and \(e=1\) is \(11\). No classification at the other values
of \(L\) is asserted by this finite example.

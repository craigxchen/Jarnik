# Good-cut contributions to moving weighted evaluation content

This note computes the exact contribution of one good two-level cut to the
moving-base denominator clearer and weighted quartic evaluation content.  It
then sums the table over a hypothetical full fair six-label core.  Finally it
gives an actual family of positive rational cotangent cliques for which the
primitive Gram form is fixed with determinant one while the evaluation
content is unbounded.

The last family rules out a bound for `C_eval` in terms of `det(Q0)` or
`Q0` alone.  It is not an endpoint counterexample: its least finite
cotangent is one, and its moving normalized nodes have growing height.  Thus
it does not rule out an endpoint inequality that also uses the node heights,
residues, or their global compatibility with the Gram form.

## 1. Setup

Use the notation of
[`moving_base_gram_content_audit.md`](moving_base_gram_content_audit.md).
Label the anchor by zero and the five finite integer cotangent numerators by
`X_1,...,X_5`, with common denominator `L`.  Select

```text
A=X_1,       B=X_2,       delta=B-A,
d=gcd(delta,A,(A^2+L^2)/delta).
```

The normalized finite coordinates are

```text
t_j=(X_j-A)/delta,
```

so `t_1=0`, `t_2=1`, and the other three coordinates are denoted `x,y,z`.
Write

```text
Q0=(a,b,c),
P0(t)=a t^2+2b t+c,
D_j=product_(k!=j) (t_j-t_k),
W_j=P0(t_j)^2/D_j,       W_infinity=-a^2.
```

If `ell` is the least positive denominator clearer for all entries
`t_j^k/D_j`, `0<=k<=4`, then

```text
Nhat f = ell (W_1,...,W_5,W_infinity),
C_eval = gcd of the six entries of Nhat f.             (1)
```

The following identity will be used repeatedly:

```text
P0(t_j)=(X_j^2+L^2)/(delta*d),
W_j=(delta^2/d^2)
    (X_j^2+L^2)^2/product_(k!=j)(X_j-X_k).             (2)
```

## 2. Exact contribution of one good cut

Let `p` be an odd prime not dividing `L`.  Suppose it gives an exact
two-level cut `S subset {1,...,5}` of positive height `e`, in the sense that

```text
v_p(X_j^2+L^2) = e if j is in S, and 0 otherwise,
v_p(X_j-X_k)   = e if j,k are both in S, and 0 otherwise.   (3)
```

These are the scalar cotangent valuations associated with a good split
prime whose Gaussian allocations have one exact two-level cut.  At a prime
outside `2L`, the upper and lower cotangent difference bounds agree, so no
unrecorded equal-level excess occurs.

Put

```text
s=|S|,       h=1 if {1,2} is contained in S, and 0 otherwise.
```

Then `v_p(delta)=h e`.  Also `d` is a `p`-unit because `d|L`.  In
particular

```text
v_p(det(Q0))=v_p((L/d)^2)=0.                           (4)
```

Formula (2) gives

```text
v_p(W_j) = [2h+(3-s) 1_(j in S)] e,
v_p(W_infinity)=2h e.                                 (5)
```

It remains to compute the denominator clearer.  If `h=0`, then
`v_p(D_j)=(s-1)e` on `S` and zero off `S`.  No power of `t_j` produces a
worse denominator, while the `k=0` entry on `S` attains this one.  Hence

```text
v_p(ell)=(s-1)e             when h=0.                 (6)
```

If `h=1`, then

```text
v_p(D_j)=(s-5)e on S,       v_p(D_j)=-4e off S.
```

The nonzero `t_j` are units on `S`, while off `S` they have valuation
`-e`.  Thus every `t_j^k/D_j`, `0<=k<=4`, is `p`-integral. If `S` is
proper, the `k=4` entry off `S` has valuation zero; if `s=5`, a `k=0`
entry has valuation zero. Therefore

```text
v_p(ell)=0                  when h=1.                 (7)
```

Since every entry in (1) is integral,

```text
v_p(C_eval)=v_p(ell)+min_j v_p(W_j),                 (8)
```

where the minimum includes the infinity row.  Equations (5)--(8) give the
complete table below.  Each entry is `v_p(C_eval)/e`.

```text
                 s=1   s=2   s=3   s=4   s=5
{1,2} not in S     0     1     2     2     --
{1,2} in S        --     2     2     1      0
```

For example, the balanced cut `S={3,4,5}` contributes `p^(2e)` to
`C_eval`, although (4) says that `p` is a unit of the Gram determinant.

## 3. The full fair-core sum

Suppose hypothetically that there is one core block of logarithmic norm
`w+o(w)` for every nonempty `S subset {1,...,5}`.  Applying the preceding
table to the exact good part gives the following contributions.

For `ell`, only cuts not containing both selected labels contribute:

```text
size 2:  9*1 =  9,
size 3:  7*2 = 14,
size 4:  2*3 =  6,
total:          29.                                  (9)
```

For `C_eval`, the size-two, size-three, and size-four totals are

```text
size 2:  9*1+1*2 = 11,
size 3:  7*2+3*2 = 20,
size 4:  2*2+3*1 =  7,
total:                  38.                          (10)
```

Thus the exact core portions have logarithmic sizes

```text
log ell_core     =29w+o(w),
log C_eval,core  =38w+o(w).                          (11)
```

This transfer is stable under the usual fixed-dimensional private-factor
errors.  Indeed, at an arbitrary prime put

```text
n_j=v_p(X_j^2+L^2),
q_jk=v_p(X_j-X_k),
q=v_p(delta).
```

Then exactly

```text
v_p(W_j)=2q-2v_p(d)+2n_j-sum_(k!=j) q_jk,
v_p(W_infinity)=2q-2v_p(d),                          (12)
v_p(D_j)=sum_(k!=j) q_jk-4q,
v_p(t_j)=q_1j-q                                      (13)
```

for the nonzero normalized entries.  The valuation of `ell` is the maximum
of zero and the finitely many forms

```text
v_p(D_j)-k v_p(t_j),       0<=k<=4,                  (14)
```

with zero entries omitted.  The valuations of `ell` and `C_eval` are
therefore fixed minima and maxima of integer linear forms in the data
(12)--(14).  Perturbing every allocation by at most `r_p`, and allowing the
cotangent error at `p|2L`, changes them by
`O(r_p+v_p(2L))`, with an absolute constant because there are only six
labels.  The central truncation bounds make the sum of these errors `o(w)`.

Equation (11) concerns the contribution of the core primes.  Extra primes
may increase the global denominator clearer or content.  Without a separate
all-edge denominator estimate, (11) must not be promoted to a global upper
bound or global asymptotic equality.

The same local calculation recovers the known primitive central profile.
After subtracting the common raw valuation in (8), each finite coefficient
gets two units of cut weight from its singleton cut, four from its incident
pair cuts, and one from the size-four cut whose complement is the anchor and
that label.  The anchor gets one unit from each of the five size-four cuts
and two from the size-five cut.  Every coefficient therefore has total
height `7w+o(w)`.  The large values in (10) are the content that is removed;
they do not supply an additional positive power in the primitive vector.

## 4. A fixed-Gram family with unbounded content

The determinant obstruction is unconditional, rather than merely a feature
of the hypothetical fair profile.  Fix `p=13`, whose square root `5` of
`-1` lifts uniquely to every power of `13`.  For `e>=1`, put `P=13^e` and
choose `rho_e` with

```text
rho_e=5 mod 13,       rho_e^2+1=0 mod P,
0<=rho_e<P.
```

There is one residue `c_bad mod 13` for which

```text
v_13((rho_e+c P)^2+1)>e.
```

Choose three consecutive integers `c_0,c_0+1,c_0+2` avoiding that residue,
and set

```text
r_e=rho_e+c_0 P,
(s_1,...,s_5)=(1,2,r_e,r_e+P,r_e+2P).               (15)
```

These are five distinct positive rational cotangents.  Let `L_e` be the
least common denominator of all ten pair cotangents

```text
(s_i s_j+1)/(s_j-s_i).
```

For two of the last three entries, the pair-cotangent denominator has exact
`13`-valuation `e`, while its numerator has valuation at least `e` because
the two entries are congruent to the same root of `-1` modulo `13^e`.
Every difference joining `1` or `2` to the last three is a `13`-unit.  Hence

```text
v_13(L_e)=0.                                         (16)
```

After putting `X_i=L_e s_i`, all pair cotangents have common denominator
`L_e`, so these are actual positive cotangent cliques.  Selecting the first
two finite nodes gives

```text
A=L_e,       B=2L_e,       delta=L_e,       d=L_e,
Q0=[[1,1],[1,2]],       det(Q0)=1.                   (17)
```

The normalized tuple is

```text
T=(0,1,r_e-1,r_e+P-1,r_e+2P-1),
P0(t)=t^2+2t+2=(t+1)^2+1.                            (18)
```

The last three nodes form the exact good cut `S={3,4,5}` of height `e`.
The table therefore gives

```text
v_13(ell)=2e,       v_13(C_eval)=2e.                 (19)
```

Thus `C_eval` is unbounded while both `Q0` and its determinant are fixed.
The growing factor comes from the moving weighted evaluation lattice, not
from the Gram determinant.  The family has `min s_i=1`, so (19) says
nothing against a lower bound restricted to shrinking endpoint
configurations.

The companion checker
[`check_moving_base_fair_content.py`](check_moving_base_fair_content.py)
verifies the local table symbolically, constructs (15) for several powers of
13, checks every pair cotangent and the least denominator, and verifies
(16)--(19) by exact integer and rational arithmetic.

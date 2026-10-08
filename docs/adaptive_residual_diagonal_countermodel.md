# Adaptive residual diagonal countermodel

## Status

This note gives an exact **abstract** countermodel to the idea that the
projective-residual lower bound can be amplified down the prime-reveal tree
using only conductor supports, oriented common cores, and non-reuse between
levels. It is not a configuration of Gaussian points on one circle: the
integer residuals are deliberately decoupled from the split-prime character
code. Its role is to identify the extra arithmetic input that a proof must
use.

---

## 1. A critical character code

Let

```text
m = 2^r,                 w = 4 log(m-1),
H = (m-1)w.
```

Index the points by the binary-reflected Gray ordering

```text
x_0, ..., x_(m-1) in F_2^r.
```

For each nonzero `a in F_2^r`, introduce one squarefree conductor coordinate
of logarithmic weight `w`, oriented by

```text
epsilon_a(x) = (-1)^(a dot x).
```

For distinct `x,y`, exactly `m/2` nonzero characters separate them. Thus the
weighted symmetric-difference distance is

```text
D(x,y) = (m/2)w = H/2+w/2.                         (1.1)
```

Attach to the ordered point `x_j` the honest primitive integer vector

```text
v_j = (1,j).
```

Its forced-reduced determinant is simply

```text
delta_ij = |det(v_i,v_j)| = j-i.                  (1.2)
```

For every rational prime `ell`, this already has the full projective
collision mechanism:

```text
ell | delta_ij  iff  [v_i] = [v_j] in P^1(F_ell). (1.3)
```

Moreover

```text
log delta_ij <= log(m-1) = w/4
             = D(x_i,x_j)/2-H/4.                 (1.4)
```

This is exactly the endpoint residual allowance with `C=2`.

---

## 2. Exact stability under renormalization

At depth `d`, write

```text
q = 2^d,                 n = m/q.
```

Each aligned `n`-point Gray block is an affine subspace. Its annihilator has
`q` characters, so exactly `q-1` of the nonzero conductor characters are
constant on it. Removing those constant factors gives

```text
L_d = log P_d = (q-1)w,
H_d = (m-q)w,
C_d = 2 exp(-L_d/4).                              (2.1)
```

Every distinct pair inside the block still has character distance `mw/2`.
Consequently its descendant endpoint allowance is

```text
(1/2)(mw/2)-H_d/4+log(C_d/2) = w/4.               (2.2)
```

Thus arbitrarily deep renormalization produces no improvement at all in the
upper bound for `log delta_ij`.

---

## 3. A coherent oriented core still does not descend

Inside a depth-`d` Gray block, the character that is the parity of all active
coordinates changes sign on every consecutive Gray edge. It has `q` ambient
extensions, giving an effective conductor class `K_d` of total weight

```text
weight(K_d) = qw > L_d.                            (3.1)
```

Hence every two consecutive-transition blocks have the same large oriented
compensator: according to their orientations, `K_d` divides either their
ordinary Gaussian gcd or the gcd after conjugating one block. Splitting on
this class removes it in both children. The next level has a disjoint core of
weight `2qw`, and

```text
sum_(d=0)^(r-1) weight(K_d)
  = sum_(d=0)^(r-1) 2^d w
  = (m-1)w
  = H.                                             (3.2)
```

So the cores can be perfectly coherent and can be charged without repetition
between levels. The obstruction is their quadratic reuse among pairs *at the
same level*.

---

## 4. The diagonal slack is saturated

For `n` consecutive integer vectors put

```text
S(n) = sum_(0 <= i < j < n) log(j-i).
```

Stirling summation gives

```text
S(n) = (1/2)n^2 log n-(3/4)n^2+O(n log n).         (4.1)
```

Summing the exact residuals over the whole binary tree gives

```text
sum_(d=0)^(r-1) 2^d S(m/2^d)
  = m^2 log m+O(m^2).                              (4.2)
```

On the other hand, summing the endpoint allowance `w/4` over the same pairs
gives

```text
(w/4) sum_(d=0)^(r-1) 2^d binom(m/2^d,2)
  = m^2 log m+O(m(log m)^2).                       (4.3)
```

The two leading terms agree. Projective collisions at every node, exact
support accounting, a growing oriented common core, and even levelwise
non-reuse are therefore all compatible with the critical diagonal slack.

---

## 5. Consequence for the proof search

This model is intentionally not a Gaussian endpoint configuration: the
residual determinants `(1.2)` do not arise from the Gaussian products encoded
by the conductor characters. That is precisely the point. Any successful
argument must use an identity coupling

```text
Gaussian phase / residual determinant
```

to

```text
the specific split-prime factors that create the character code.
```

Support widths, projective collision counts, adaptive renormalization, and
oriented compensator bookkeeping cannot prove the uniform bound in isolation,
even when all four are combined. The next live target is therefore a
phase-sensitive cocycle or continuant statement that prevents one conductor
core from paying for quadratically many residual pairs at a fixed level.

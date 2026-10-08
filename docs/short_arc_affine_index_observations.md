# Bounded affine-index observations for five and six point windows

This is a finite computational supplement to the affine-shape and conic
notes. It does not extract a bounded affine type from an arbitrary endpoint
cluster.

For a selected `m`-point subset of one integer circle, let `D_I` be its
triangle determinants, `g=gcd_I |D_I|`, and `H=max_I |D_I/g|`. The existing
divisibility theorem gives

```text
N^b_m | product_I |D_I/g|,
b_m=binom(floor(m/2),3)+binom(ceil(m/2),3).       (1)
```

The checker [check_short_arc_affine_index.py](check_short_arc_affine_index.py)
enumerates complete circles for several small norms, scans cyclic windows of
five and six points whose chord span is at most `20 sqrt(R)`, and verifies (1)
exactly. It also reports the smallest observed pair `(H,g)`.

The finite output is:

```text
N=5:   m=5 (H,g)=(5,2),  m=6 (6,2)
N=25:  m=5 (8,2),        m=6 (15,2)
N=65:  m=5 (13,2),       m=6 (18,2)
N=85:  m=5 (17,2),       m=6 (21,2)
N=125: m=5 (25,2),       m=6 (39,2)
N=325: m=5 (15,2),       m=6 (30,2)
N=425: m=5 (20,2),       m=6 (30,2)
```

The repeated index `g=2` in these fixtures reflects the common parity of
integer circle chords; it is bounded here but no universal bounded-index claim
is made. The values of `H` do not remain bounded as the norm grows, and the
finite scan gives no evidence for a radius-independent affine shape theorem.
At the actual endpoint scale `C sqrt(R)` with fixed small `C`, the existing
triangle bound is `H <= C^3 sqrt(R)/(8g)`, while (1) has exponent `b_m` large
enough that substitution leaves a positive power of `R`. Thus five or six
points alone do not close the endpoint count problem through this route.

The fixed Pell quadrilateral is qualitatively different: its unimodular frame
and two fixed conics isolate the whole endpoint neighborhood. The present
calculation only checks the general triangle-content identity on dense finite
windows and does not imply that an arbitrary five- or six-point cluster has a
small-index frame or fixed conic pencil.

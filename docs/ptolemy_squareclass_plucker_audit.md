# Ptolemy squareclasses and global Plucker coordinates

## Status

This note classifies the all-roots squarefree compatibility carried by any
cyclic family of rational lattice points on one circle.

The classification is exact and global.  It shows that Ptolemy forces every
chord squareclass to factor into one common class and one label per vertex.
After removing those labels, the positive chord radicals satisfy the rank-two
Plucker relations.

The result does **not** prove the uniform endpoint bound.  The final audit shows
that the squareclass/Plucker system is equivalent to rational half-angle
parametrization unless one adds further divisor or size information.

---

## 1. Rooted primitive coordinates

Let the chosen root point be

```text
z_0 = h(A+iB),
N = A^2+B^2,
gcd(A,B)=1.
```

Thus the actual circle has squared radius `h^2 N`.

For another point write its displacement from the root as `(m_i,n_i)` and set

```text
k_i = -(A m_i+B n_i),
s_i = A n_i-B m_i.
```

The equal-circle equation gives

```text
k_i^2+s_i^2 = 2 h N k_i.                         (1.1)
```

Put

```text
g_i = gcd(k_i,s_i),
k_i = g_i q_i,
s_i = g_i r_i,
gcd(q_i,r_i)=1.
```

Equation (1.1) implies `q_i | g_i`; write `g_i=lambda_i q_i`.  Then

```text
k_i = lambda_i q_i^2,
s_i = lambda_i q_i r_i,
D_i = q_i^2+r_i^2 = 2 h N/lambda_i.                (1.2)
```

The pair `(q_i,r_i)` is the primitive half-angle vector of the point relative
to the root.

---

## 2. Exact chord formula

For two non-root points define

```text
Delta_ij = q_i r_j-q_j r_i.
```

A direct half-angle calculation gives

```text
d_ij := |z_i-z_j|^2
      = 4 h^2 N Delta_ij^2/(D_i D_j)
      = (lambda_i lambda_j/N) Delta_ij^2.           (2.1)
```

For the root one may formally take

```text
(q_0,r_0)=(0,1),
D_0=1,
lambda_0=2 h N,
```

and the same formula gives

```text
d_0i = 2 h lambda_i q_i^2.                          (2.2)
```

Hence, in the squareclass group `Q^*/Q^{*2}`,

```text
[d_ij] = [N D_i D_j].                               (2.3)
```

Thus all edge squareclasses factor through one common label `[N]` and one
vertex label `[D_i]`.

---

## 3. Abstract squareclass classification

Let `V` be any vector space over `F_2`.  Label each edge of a complete graph by
`v_ij in V`.  Assume that for every four distinct vertices

```text
v_ik+v_jl = v_ij+v_kl = v_il+v_jk.                 (3.1)
```

Then there exist one common class `c in V` and vertex labels `x_i in V` such
that

```text
v_ij = c+x_i+x_j.                                   (3.2)
```

Proof.  Choose vertices `1,2,3` and put

```text
c = v_12+v_13+v_23,
x_1=0,
x_i=c+v_1i.
```

The four-point relations first show that every triangle containing vertices
`1,2` has total label `c`; applying (3.1) to `(1,2,i,j)` then gives

```text
v_ij=c+v_1i+v_1j=c+x_i+x_j.
```

Therefore (2.3) is not an accidental consequence of the chosen root: it is the
unique general form of every Ptolemy-compatible squareclass system.

---

## 4. Ptolemy becomes Plucker

For cyclically ordered points, Ptolemy says

```text
sqrt(d_ik d_jl)
  = sqrt(d_ij d_kl)+sqrt(d_il d_jk),                (4.1)
```

for `i<j<k<l`.

Substituting (2.1), all vertex factors cancel, leaving

```text
Delta_ik Delta_jl
  = Delta_ij Delta_kl+Delta_il Delta_jk.            (4.2)
```

This is the rank-two Plucker relation for the primitive vectors

```text
v_i=(q_i,r_i).
```

The polynomial identity (4.2) is formalized in

```text
GaussianChain/PtolemyPlucker.lean.
```

---

## 5. Converse reconstruction

The Plucker relations by themselves contain no hidden higher-dimensional
constraint.  Suppose positive rational numbers `Delta_ij` satisfy (4.2) and
`Delta_12 != 0`.  Set

```text
v_1=(1,0),
v_2=(0,Delta_12),
v_i=(-Delta_2i/Delta_12, Delta_1i).
```

Then (4.2) gives

```text
det(v_i,v_j)=Delta_ij
```

for every pair.  Clearing one common denominator produces integral vectors,
at the cost of one global scale.

Thus the full positive Plucker system is equivalent to a family of rational
slopes.  It is exactly the algebra behind rational half-angle parametrization
of one circle.

---

## 6. What remains arithmetic

The information not captured by abstract Ptolemy/Plucker compatibility is:

1. each primitive vector satisfies

   ```text
   q_i^2+r_i^2=D_i,
   D_i | 2 h N;
   ```

2. the scaling factor is fixed by

   ```text
   lambda_i D_i=2 h N;
   ```

3. every chord square is the integer

   ```text
   d_ij=(lambda_i lambda_j/N) Delta_ij^2;
   ```

4. the points lie in one endpoint sector, so the slopes `q_i/r_i` occupy an
   interval of length `O(N^{-1/4})` in angular coordinates.

Any successful continuation must combine several of these conditions.  The
Ptolemy equations alone cannot bound the number of points because arbitrary
finite rational slope sets satisfy them.

---

## 7. Verdict

The all-roots Ptolemy analysis produces a useful exact normal form:

```text
chord squareclass = common class * vertex class_i * vertex class_j,
normalized chord radical = rank-two Plucker coordinate.
```

It does not create the missing endpoint slack.  Its value is organizational:
it reduces every proposed global chord invariant to the primitive vectors
`(q_i,r_i)` and the divisor labels `D_i | 2hN`.  A genuinely new theorem must
control many primitive sum-of-two-squares divisors `D_i` whose rational slopes
lie in one `N^{-1/4}` sector.

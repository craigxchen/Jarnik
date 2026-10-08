# Isotropic marked directions: exact contacts and toric heights

This note retains the two conjugate isotropic directions of an integral determinant-one frame. It gives their contact dictionary and an exact height formula. The formula shows that standard toric heights and linear tropical balancing of this lift are already determined by the original Gaussian cut weights. No curve extraction or fixed prime support is assumed.

## 1. The two additional marks

Let `U=[[a,c],[b,d]] in SL_2(Z)`, write `z=a+ib`, `w=c+id`, and let `v_j=(r_j,s_j)` be distinct primitive integer source directions, for `1<=j<=n`. Put

```text
P_j=r_j z+s_j w,
A=(-w,z),             B=(-bar(w),bar(z)),
alpha=[A],             beta=[B].
```

Use `[V,W]=det(V,W)`. Direct calculation gives

```text
[A,B]=-2i,
[A,v_j]=-P_j,          [B,v_j]=-bar(P_j),               (1)
P_i bar(P_j)-P_j bar(P_i)=-2i[v_i,v_j].                (2)
```

The sign in (2) follows, for instance, by taking `U=identity`: then `P_j=r_j+i s_j`. Both `A` and `B` are primitive Gaussian rows. Indeed, `d z-b w=1`, so `gcd_G(z,w)=1`. The real rows are also Gaussian-primitive by ordinary integer Bezout. The marks `alpha,beta` are distinct and nonreal, and conjugation interchanges them while fixing each source direction. These data therefore define a point of `M_(0,n+2)(Q(i))`.

The rational coordinate

```text
t_j=[A,v_j]/[B,v_j]=P_j/bar(P_j)                       (3)
```

sends `alpha,beta` to `0,infinity` and satisfies `bar(t_j)=t_j^(-1)`. The remaining ambiguity is simultaneous multiplication of all `t_j` by one nonzero scalar. Thus the corresponding torus is `(G_m)^n/G_m`, with character lattice

```text
L={a in Z^n : sum_j a_j=0}.
```

In particular the ordinary marked cross-ratios `t_i/t_j` are exactly the original relative Gaussian phases. Equation (2) shows that the non-toric source-collision functions retain the original row determinants as well.

## 2. Prime-power boundary dictionary

At an odd Gaussian prime `pi`, let

```text
n_j=v_pi(t_j)=v_pi(P_j)-v_bar(pi)(P_j).
```

Because `U` is integral unimodular and `v_j` is primitive, `P_j` and `bar(P_j)` cannot both vanish at `pi`. Since `[A,B]` is a unit here, `n_j>0` means that `v_j` meets `alpha` to depth `n_j`, while `n_j<0` means that it meets `beta` to depth `-n_j`.

Modulo simultaneous scaling, the valuation vector is the class of `(n_1,...,n_n)` in `R^n/R(1,...,1)`. Its successive distinct levels give nested subsets

```text
S_t={j:n_j>=t}.
```

The gap between consecutive levels is the length of the corresponding segment separating the two distinguished marks in stable reduction. For a nonempty proper subset `S`, the associated boundary partition is

```text
D_S = D_( {alpha} union S | {beta} union S^c ).         (4)
```

This can be checked by rescaling the `t` coordinate across one gap: precisely the marks in `S`, together with `alpha=0`, lie in the vanishing residue disk; the remaining marks, together with `beta=infinity`, lie on the other side. Further coincidences among equal-valuation source marks may create additional branches. They do not change these segments.

In the exact independent-block model

```text
P_j=product_(S containing j) H_S,
```

with pairwise disjoint odd split Gaussian prime supports, also disjoint from all conjugate supports, a prime power `pi^e` in `H_S` gives the vector `e 1_S`. At `bar(pi)` it gives `-e 1_S`, equivalently `e 1_(S^c)` modulo the constant vector. Thus the block produces `D_S` at one conjugate prime and `D_(S^c)` at the other, each to depth `e`.

Singleton cuts are genuine boundaries after the two marks are retained: `D_{j}` separates `alpha` and source mark `j` from the other marks. Empty and full cuts are different. Their valuation vectors are constant, so they disappear under simultaneous scaling; the purported partition would have a singleton side and would not be a stable boundary. A full-support Gaussian factor rotates every `t_j` by the same factor `H/bar(H)`. It changes none of the characters below.

The label count matters. Here there are exactly `n` source marks and two isotropic marks. If an original reference row with Gaussian value `P_0=1` is retained in addition to `n` active rows, the space is instead `M_(0,n+3)`. The cut containing all active rows but excluding that reference is then a proper cut and remains visible; its character values relative to `t_0=1` do change. No full-support weight relative to an omitted external reference may be declared irrelevant to the original circle or radius merely because it disappears from the smaller marked space.

## 3. Exact toric-height theorem

Let `mathcal A` be any fixed finite nonempty subset of `L`. Consider the projective monomial map

```text
F_A(t)=[t^a]_(a in mathcal A),
width_A(v)=max_(a in mathcal A) a dot v
             -min_(a in mathcal A) a dot v.
```

For the exact independent-block model, put `w_S=log Norm(H_S)`. Then the absolute logarithmic projective height satisfies the exact identity

```text
h(F_A(t)) = (1/2) sum_S w_S width_A(1_S).              (5)
```

**Proof.** At the complex place, all `|t_j|=1`, so all monomial coordinates have absolute value one and the local height is zero. At `pi^e || H_S`, the coordinate valuations are `e(a dot 1_S)`; at the conjugate prime they are their negatives. In the absolute height over `Q(i)`, each finite place contributes its logarithmic norm with weight `1/2`. Pairing the two places therefore contributes

```text
(e log Norm(pi)/2) width_A(1_S).
```

Summing proves (5). There are no ramified or inert block primes under the stated hypotheses.

Since every exponent has coordinate sum zero,

```text
width_A(1_empty)=width_A(1_full)=0,
width_A(1_(S^c))=width_A(1_S).                         (6)
```

Thus the formula retains exactly the complementary paired weight and discards exactly the common-scaling weight; it does not count an empty or full cut as an extra boundary.

For the standard projective torus coordinates, take

```text
mathcal A={0,e_1-e_n,...,e_(n-1)-e_n}.
```

Every nonempty proper cut has width one, giving

```text
h([1:t_1/t_n:...:t_(n-1)/t_n])
    = (1/2) sum_(empty != S != full) log Norm(H_S).    (7)
```

For equal weights `w`, this is `(2^(n-1)-1)w`. In contrast, a single relative phase has height `2^(n-2)w`, obtained by using just the exponents `0,e_i-e_j` in (5).

The formula also applies asymptotically to the actual corrected extraction. If

```text
P_j=K_j product_(S containing j) H_S,
```

with nonzero Gaussian integers `K_j`, then, without any disjoint-support
assumption on the `K_j`,

```text
|h(F_A(t))-(1/2)sum_S w_S width_A(1_S)|
    <= C_A sum_j log|K_j|.                            (8)
```

To prove this, add the signed correction valuation vectors and use
`|width_A(v+u)-width_A(v)|<=width_A(u)`. The weighted absolute signed valuations of each `K_j` are at most `log Norm(K_j)`; the ramified and inert primes have zero signed valuation. At infinity the corrected ratios still have modulus one. Hence the coefficients in (8) depend only on the fixed exponent set. For a fixed extracted row count with `log|K_j|=o(w)`, the error is `o(w)`.

## 4. What the fixed determinant actually contributes

Let `H_U=sqrt(|z|^2+|w|^2)` and `H_j=sqrt(r_j^2+s_j^2)`. Use Euclidean chordal distance at the complex place, and let `lambda=-log(distance)`. Equations (1) give

```text
lambda_infinity(alpha,beta)=2log H_U-log 2,
lambda_infinity(alpha,v_j)
 =lambda_infinity(beta,v_j)=log H_U+log H_j-log|P_j|.  (9)
```

For two primitive Gaussian rows `V,W`, define their finite contact sum by

```text
N_fin(V,W)=(1/2)sum_(Gaussian prime ideals pi)
                       v_pi([V,W]) log Norm(pi).
```

Again using (1),

```text
N_fin(alpha,beta)=log 2,
N_fin(alpha,v_j)=N_fin(beta,v_j)=log|P_j|.              (10)
```

Thus the fixed determinant puts the mutual finite contact of the distinguished pair entirely over two. Its growing mutual proximity is archimedean. Combining (9) and (10) yields the exact height balance

```text
N_fin(alpha,v_j)+N_fin(beta,v_j)
 +lambda_infinity(alpha,v_j)+lambda_infinity(beta,v_j)
    =2log H_U+2log H_j.                               (11)
```

There is no negative correction from `[A,B]=-2i` that can simply be subtracted from every source contact. For example, `U=identity` and `v=(N,1)` keep the distinguished pair fixed, while the total finite contact in (10) is `log(N^2+1)` and tends to infinity. The source height on the right of (11) pays for it exactly. This example only tests that proposed height subtraction; it is not a full-profile or endpoint counterexample.

## 5. Exact scope of the obstruction

Conjugate places carry opposite valuation vectors with equal logarithmic weights. Therefore their vector sums cancel automatically, for every choice of nonnegative cut weights. Linear tropical balancing supplies no further restriction on those weights.

Likewise, (5) proves that projective monomial heights, including the standard heights furnished by nef integral toric polytopes, are already fixed by the old cut masses. Their positivity yields nonnegative combinations of the same masses. Retaining the distinguished pair does not create an additional positive height margin in these calculations.

This conclusion is restricted to the separating-boundary/toric calculation. The full moduli compactification also remembers source-only collision branches, and arithmetic arguments involving them can use the determinants in (2), hence the original primitive residues. The note does not exclude a new inequality using that information together with (9), or a non-toric auxiliary construction. No stronger bound on the endpoint point count has been proved here.

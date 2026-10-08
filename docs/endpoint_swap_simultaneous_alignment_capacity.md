# Simultaneous fixed-parameter overlap has a finite capacity

The four-point Pell family's preservation of endpoint scale for every fixed
positive rational alignment cannot extend in that strong form to the
families with at least five points considered here. A Gaussian resultant
bounds the product of target overlaps over distinct integer parameters.
Combined with the compulsory overlap bound and a five-row restriction,
it allows at most twelve fixed integer parameters to preserve endpoint
scale simultaneously along one family.

This is a limitation theorem. It does not rule out one successful
nonconformal alignment or an adaptive choice, and it does not improve the
general lattice-point count.

## 1. Setup and exact divisibility

Use the notation of the
[content and overlap theorem](endpoint_swap_content_and_overlap.md).
Fix even `x>0`, set `T=x^2+1`, `B=x+i`, and take interior source rows

```text
H_i=d_i+v_i B,     d_i,v_i>0, gcd(d_i,v_i)=1, i=1,...,r.
```

The interior phases are distinct. Include the endpoint rows `1,B`. No
condition `d_i|T` is imposed. Put `Delta=2 arctan(1/x)` and let `N` be
the source's least squared radius. In the `B`-oriented frame, let `Q_i`
be the reduced source denominator after removing its gcd with `B`. Thus

```text
N=T Norm(lcm_G(Q_1,...,Q_r)).
```

For positive integer `k` coprime to `T`, the aligned target uses
`C_i(k)=d_i+k v_i B`. Let `Q_i(k)` be its reduced excess denominator
and define its full overlap

```text
Omega_k=product_i Norm(Q_i(k))/Norm(lcm_G(Q_1(k),...,Q_r(k))).
```

Choose a finite nonempty set `S` of such parameters, and put `K=max S`.
Define Gaussian integers

```text
F_i=gcd_G(B,d_i),
R_ij=B(d_i v_j-d_j v_i)/(F_i F_j).
```

Then there is an exact divisibility of positive ordinary integers:

```text
product_(k in S) Omega_k
  | [(K-1)!]^(r(r-1)) product_(i<j) Norm(R_ij).          (1)
```

Every `R_ij` is nonzero because the interior phases are distinct.

## 2. Primewise proof, including contents and repeated powers

The Gaussian integral linear polynomial

```text
L_i(X)=(d_i+X v_i B)/F_i
```

has coprime coefficients: `gcd_G(d_i,v_i B)=F_i`. The target
`Q_i(k)` divides `L_i(k)`. Indeed ordinary image content is
`gamma_i=gcd(d_i,k)`, which is coprime to `T`; ramified reduction is
also coprime to `F_i`, since `T` is odd. The anchor gcd after those
reductions is still `F_i`, up to a unit. Thus this divisibility retains
all ordinary contents and the factor above two.

At any Gaussian prime, the product/lcm surplus `sum a_i-max a_i` is
nondecreasing in each input valuation. Consequently `Omega_k` divides
the norm of the product/lcm surplus of the `L_i(k)`. That Gaussian
surplus itself divides

```text
product_(i<j) gcd_G(L_i(k),L_j(k)).
```

The determinant of the coefficient pairs of `L_i,L_j` is `R_ij`, so
`R_ij` is integral. Fix a Gaussian prime `pi`, and put

```text
b_k=v_pi(gcd_G(L_i(k),L_j(k))).
```

These valuations satisfy `b_k<=v_pi(R_ij)`. If any is positive, the
linear coefficient of `L_i` is a `pi`-unit, since its two coefficients
are coprime. Choose `k_0` attaining the largest `b_k`. Subtracting
`L_i(k)-L_i(k_0)` gives, for every other `k`,

```text
b_k<=v_pi(k-k_0).
```

Moreover `product_(k in S, k!=k_0) |k-k_0|` divides
`(k_0-1)!(K-k_0)!`, which divides `(K-1)!`. Hence

```text
sum_(k in S) b_k <= v_pi(R_ij)+v_pi((K-1)!),
product_(k in S) gcd_G(L_i(k),L_j(k)) | R_ij (K-1)!.
```

This includes split, inert and ramified primes, with every multiplicity.
Multiplying over pairs and taking norms proves (1).

## 3. A source-height upper bound

Write `h_i=Norm(F_i)=gcd(d_i,T)` and `epsilon_i in {1,2}` for the
source ramified reduction. All source half-angle arguments lie between
zero and `arctan(1/x)`. Their sine difference therefore gives

```text
T(d_i v_j-d_j v_i)^2 <= Norm(H_i) Norm(H_j).
```

It follows that

```text
Norm(R_ij) <= epsilon_i epsilon_j Norm(Q_i) Norm(Q_j)
           <=4 Norm(Q_i) Norm(Q_j).
```

Combining this with (1) yields

```text
product_(k in S) Omega_k
 <= [2(K-1)!]^(r(r-1)) [product_i Norm(Q_i)]^(r-1)
 <= [2(K-1)!]^(r(r-1)) (N/T)^(r(r-1)).                 (2)
```

The last step uses that each `Q_i` divides their Gaussian lcm. No
squarefree, disjoint-support, or small-content assumption is used.

## 4. Consequence for fixed-size families

Fix `r>=3` and a source family with `x` tending to infinity and bounded
normalized width `Delta N^(1/4)`. Suppose `s` distinct fixed positive
integer parameters are coprime to `T` throughout the family and each
has bounded target normalized width along the whole family. Then

```text
s(r-2) <= 2r(r-1),
s <= floor(2r(r-1)/(r-2)).                             (3)
```

Indeed `Delta>=1/x` and `T>=x^2` give `N/T=O(x^2)`, with a constant
depending only on the source width bound. Equation (2) has upper growth
`O(x^(2r(r-1)))`. Equation (10) of the content and overlap theorem gives
`Omega_k >= c_k x^(r-2)` for each successful parameter, where `c_k>0`
is fixed by its target width bound. Multiplication gives the lower growth
`x^(s(r-2))`, proving (3).

This count includes the conformal parameter `k=1` if selected. For
example, all thirteen fixed integer alignments cannot preserve endpoint
scale along a five-point family. Choosing powers of two makes their
coprimality with the odd integer `T` automatic. In contrast, `r=2`
provides no restriction, consistently with the
[four-point Pell family](endpoint_swap_arithmetic_sensitivity.md).

The statement concerns simultaneous boundedness along a family, not the
number of successful parameters at any one finite source. Failure of
simultaneous boundedness also need not mean that one target width tends
to infinity along the entire family. An adaptive parameter or one fixed
successful nonconformal parameter remains possible.

## 5. Five-row restriction gives an absolute capacity of twelve

The number of source rows need not remain fixed. Suppose every source in
a family has at least five points, the two endpoints stay in the stated
`1,B` form, `x` tends to infinity, and source normalized widths are
uniformly bounded. Then at most **twelve distinct fixed positive integer
alignments** coprime to `T` throughout the family can each have uniformly
bounded target normalized width. The target bound may depend on the
fixed parameter. This count includes `k=1` when selected.

To prove this, retain the two endpoints and any three distinct interior
rows from each source. The chosen rows may vary with the source. The
Gaussian denominator lcm of a subset divides that of the full tuple;
therefore both its least source squared radius and each least target
squared radius divide the corresponding full squared radius. The common
angular width `Delta` is unchanged because both endpoints are retained.
Thus the source and target normalized widths of the selected five-row
tuples are bounded by those of the full tuples. Apply (3) with `r=3`:

```text
s <= 2*3*(3-1)/(3-2) = 12.                            (4)
```

In particular, any thirteen fixed powers-of-two alignments cannot all
preserve endpoint scale along such a family, regardless of how its
point count varies. This bounds simultaneous map preservation, not the
source point count or the existence of one useful nonconformal map.

## Verification and scope

The [exact checker](check_endpoint_swap_simultaneous_alignment_capacity.py)
passes 804 configurations. It checks both the Gaussian pair divisibility
and the full integer divisibility (1), independent primitive source and
target radii, the bounds in (2), nonconsecutive parameter sets, nontrivial
ordinary target contents, and deliberate shared split-prime powers
`13^3` at both a source and a target parameter. For every tested tuple
and parameter it also verifies that all endpoint-preserving five-row
subset squared radii divide their full source or target squared radius.
These finite checks supplement the primewise proof. Astra supplied the
resultant argument, and root and Sol independently audited it. Root
derived the five-row strengthening, independently checked by both agents.

A lower bound guaranteeing adequate overlap for some nonconformal
parameter is still missing. Neither (1)--(4) nor the content formulas
supply that existence result or the desired uniform endpoint count.

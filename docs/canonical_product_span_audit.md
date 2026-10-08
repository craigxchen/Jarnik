# A nonsaturated canonical-product span

The canonical weighted products from
[higher_degree_gram_jet_barrier.md](higher_degree_gram_jet_barrier.md) need
not generate the full contact-jet lattice.  This is a lattice-generation
issue, separate from the product-height barrier.

Use the identity frame and three singleton contacts with primitive source
rows
\[
v_0=(1,2),\qquad v_1=(4,5),\qquad v_2=(7,8).
\]
Their pairwise determinants are (-3,-6,-3), while the contact norms are

\[
5=1^2+2^2,\qquad 41=4^2+5^2,\qquad 113=7^2+8^2.
\]

The contact norms are odd and pairwise coprime, and each contact vector is
the corresponding source row.  For degree (d) and jet order (h), form
all canonical products
\[
P_n=\prod_i d_i^{\max(h-n_i,0)}L_i^{n_i},
\qquad \sum_i n_i=d.
\]
Exact integer minor computations for (h=1,2) and (d=1,2,3) give
\[
[\Gamma_{d,h}:\langle P_n\rangle]
 =3^{d(d+1)/2}.
\]
Thus the canonical products are nonsaturated.  The obstruction is the
common factor (3) in the source determinants: modulo (3), all three
linear forms are proportional, producing the triangular exponent in
degree (d).  Their evaluations are nonzero, so the products themselves
do not all lie in the Gram kernel.  In this example the evaluation ideal
is also nonsaturated for ((h,d)=(2,1)), while it has the expected
contact index for the other tested pairs; image saturation and polynomial
saturation are separate questions.
The exceptional pair `h=2,d=1` lies outside the hypothesis `d>=h`
of the exact Gaussian jet-image theorem, so it does not contradict
that theorem. The relative polynomial index here is explained for
every prime away from the contact norms by the
[Smith-invariant formula](prime_to_D_smith_saturation_audit.md).

The checker
[check_canonical_product_span.py](check_canonical_product_span.py) uses
exact integer determinants, rational solves, and Gaussian Euclidean gcds;
it does not infer saturation from floating-point rank or from numerical
evaluation.

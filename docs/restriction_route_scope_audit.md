# Restriction-theory scope audit — 2026-09-05

[Demeter–Langowski, Section 13, Theorem 13.2](https://academic.oup.com/imrn/article/2023/2/1292/6390751)
proves a radius-independent fourth-moment restriction estimate assuming
a uniform short-arc bound at some radius exponent `gamma>1/2`.
Its notation is `R=sqrt(N)`: arc length `N^(gamma/2)` means `R^gamma`.
The theorem is conditional and cannot supply the endpoint assertion
being pursued here.

The proof controls dyadic scales `2^j<=R^(2 gamma-1)` using that
hypothesis. At `gamma=1/2`, this range collapses to `2^j<=1`.
Substituting the endpoint into the stated argument therefore removes
the growing scale range needed to absorb its divisor-bound loss.

A separate calculation shows the limitation of straightforward
endpoint partitioning: an arc of length `sqrt(2^j R)` needs
`O(2^(j/2))` endpoint pieces. The resulting pair-count bound is
`O(2^j)`, so summing against `2^(-j beta)` yields a uniform bound
only when `beta>1`. This does not rule out a better analytic argument.

The apparent disproof in [arXiv:2107.09991](https://arxiv.org/abs/2107.09991)
remains withdrawn because of an error in Lemma 2. The search found
no applicable unconditional replacement; search failure alone is
not proof that none exists. The endpoint goal remains unproved.

## September 21 update

[Zhang–Zhu, June 2026](https://arxiv.org/html/2606.08650v1), Introduction
and Theorem 2.5, gives a lossless comparison between totally geodesic
restriction norms and lattice-point counts in unit-width spherical
bands. Its new higher-dimensional restriction bounds do not settle the
two-dimensional endpoint: the introduction explicitly retains the
square-root arc-count problem as open. This distinguishes a useful
analytic reformulation from an applicable uniform estimate. In
particular, one cannot import a bound for curved segments or closed
rational geodesics as a bound uniform over all straight-segment
directions.

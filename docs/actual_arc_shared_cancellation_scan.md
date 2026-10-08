# Actual short-arc scan: shared cancellation and gap unevenness

This note records a targeted finite scan of actual Gaussian lattice points.
It tests whether short consecutive circle windows force either a large
shared-prime cancellation deficit in coupled products or very uneven angular
gaps. The scan finds neither implication at the tested scale. It is finite
evidence only and contains no endpoint counterexample or uniform theorem.

## Exact quantities

For a primitive equal-norm tuple z_1,...,z_m, let N=|z_i|^2. For each split
prime p=pi*bar(pi) dividing N, write

~~~
e = v_pi(z_i)+v_barpi(z_i),
a_i = v_pi(z_i).
~~~

The raw source weight of an unordered pair i,j at p is e log p. After
conjugate cancellation, its directional factor has weight

~~~
|a_i-a_j| log p.
~~~

Thus the exact deficit is (e-|a_i-a_j|) log p. The reported aggregate
cancellation fraction is the total deficit divided by the total raw pair
weight. This is the pair specialization of the cancellation formula in
[coupled_block_roth_packing.md](coupled_block_roth_packing.md), evaluated
from actual Gaussian valuations rather than a fair valuation model.

For a consecutive angular window, the scan reports

~~~
C_arc = N^(1/4) * angular_width,
gap_ratio = largest consecutive angular gap / smallest consecutive gap.
~~~

The actual common Gaussian gcd of each window is removed before computing the
intrinsic N and C_arc. Angles are used only for exploratory ranking; Gaussian
gcds, divisions, and prime valuations are exact integers.

## Scan output

The standard-library scanner
[check_actual_arc_shared_cancellation.py](check_actual_arc_shared_cancellation.py)
enumerates all lattice points with |x|,|y|<=100, and examines windows of five
through eight points. For windows with intrinsic C_arc<=20, the output was:

~~~
m=5: minimum-C window 3.71014; C<=4 max cancellation 0.400000,
     minimum gap ratio 2.535011; C<=8 max cancellation 0.533333,
     minimum gap ratio 1.064087.
m=6: minimum-C window 4.25882; C<=8 max cancellation 0.577778,
     minimum gap ratio 1.151161.
m=7: minimum-C window 5.54795; C<=8 max cancellation 0.523810,
     minimum gap ratio 1.441016.
m=8: minimum-C window 5.91477; C<=8 max cancellation 0.478262,
     minimum gap ratio 2.141064.
~~~

The fair disjoint binary benchmark has cancellation fractions

~~~
(2^(m-2)-1)/(2^(m-1)-1)
~~~

equal to 0.466667, 0.483871, 0.492063, 0.496063 for m=5,6,7,8. Actual
windows can lie on either side of this benchmark. The minimum gap ratios near
one at larger C_arc show that cancellation and gap unevenness are not
interchangeable finite proxies.

## An exact eight-point fixture

The scan's minimum-C eight-point window is

~~~
N = 5525 = 5^2 * 13 * 17,
z = (70-25i, 71-22i, 73-14i, 74-7i,
     74+7i, 73+14i, 71+22i, 70+25i).
~~~

Its common Gaussian gcd is one. Choosing pi_5=1+2i, pi_13=2+3i, and
pi_17=1+4i, the exact allocation vectors are

~~~
5^2:  (1,0,2,2,0,0,2,1),
13:   (1,0,1,0,1,0,1,0),
17:   (0,1,1,0,1,0,0,1).
~~~

Summed over the 28 unordered pairs, the raw and reduced exponents and their
deficits are

~~~
prime   raw   reduced   deficit
5^2       56      30        26
13        28      16        12
17        28      16        12
~~~

Therefore the exact aggregate deficit is

~~~
26 log(5) + 12 log(13) + 12 log(17),
~~~

over raw weight 56 log(5)+28 log(13)+28 log(17), giving numerical fraction
0.441912. The angular window is symmetric with width 2 atan(5/14), hence
C_arc=5.91477... . Consecutive half-angle tangents are

~~~
1/47, 1/18, 1/21, 7/74, 1/21, 1/18, 1/47,
~~~

so the tangent-gap ratio is exactly 329/74; the angle-gap ratio is
4.433425... . This is an actual primitive tuple, but its C_arc is well above
the endpoint constant 1/2, and it is not a counterexample to the uniform
endpoint claim.

The scan's only robust conclusion is negative: these bounded actual records
do not justify replacing simultaneous phase compatibility by a rule based
solely on shared-prime cancellation or local gap ratios. A larger search or
a proof would need to use the actual Gaussian arguments and their joint
short-arc constraint.

## Verification boundary

Run:

~~~
python3 -B docs/check_actual_arc_shared_cancellation.py
~~~

The default bound is 500 and is intended for targeted exploratory use. The
reported table above uses bound=100, max_rows=8, and cutoff=20. The search is
not a proof and does not assert that any fair valuation model is an actual
endpoint configuration.

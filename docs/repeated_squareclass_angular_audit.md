# Repeated pair-squareclass angular audit

The exact three-partition proof in
[the independent audit](angular_pair_squareclass_independent_audit.md)
shows that every pair-norm squareclass is distinct in a fixed-unit C<=2
cluster. This extends the direct constant range of the older four-wise
baseline in [squareclass_baseline_and_power_groups.md](squareclass_baseline_and_power_groups.md),
Section 3. It does not improve the general growth rate.

The initial pairwise estimate `Im(AB)=ad+bc>=a+c`, for positive integer
imaginary parts, applies whether or not the squareclasses repeat. It gives
no exclusion on its own. The three balanced partitions and their joint
norm budget are essential to the proved result.

## Bounded numerical checks

[The allocation checker](check_repeated_squareclass_angular_audit.py)
audited 3,135 systems; 2,226 exhibited disjoint repeated labels with distinct
ordinary norms, but none of its tested repeated-label four-windows had
C<=2. The actual Pell fixture N=358853785 has a four-window with
C approximately 1.91156867551 and six distinct pair squareclasses.

[The direct circle checker](check_repeated_squareclass_circle_audit.py)
scanned odd squared radii N<=100000. No tested four-window with C<=2 had
repeated labels. The closest repeated-label window had
N=61625, C approximately 2.61990584516, and points
`(236,77),(235,80),(224,107),(220,115)`. Its repeated squareclass
`(5,17)` occurs on disjoint edges with ordinary norms 2125 and 85.
The scan does not impose a common literal Gaussian-unit class, which must
be retained when comparing it to the proved fixed-unit theorem.

Point coordinates and norms are exact integers. The window ordering and
reported constants use floating point, so these searches are bounded
numerical evidence, not a proof. The general proof and its normalization
conditions are in the independent audit linked above.

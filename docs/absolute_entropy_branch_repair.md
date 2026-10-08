# Absolute-entropy prime reveal and branch-repair renormalization

## Status

This note strengthens the prime-reveal framework by replacing the rare-event
information cost `log(R_N/M)` with an exact identity equal to `log M`, the
absolute logarithmic cluster size. It then packages the remaining arithmetic
problem as a branch-repair hypergraph. This is not yet a proof of the uniform
endpoint bound.

## 1. The successful-leaf process

Assume first that the split support is squarefree. Let

    A subset { -1,+1 }^r

be the set of orientation vectors whose angles lie in one endpoint interval,
and let `|A| = M`. Sample one successful vector `Epsilon` uniformly from A.
Choose an ordering of the prime coordinates and reveal them successively.

For a revealed prefix `omega` at time k, let

    N_k(omega) = number of successful completions of omega in A.

Thus `N_0 = M` and `N_r = 1` along every successful leaf. If the two children
of a prefix have completion counts `n_+` and `n_-`, the conditional law of the
next sign under the uniform successful-leaf measure is

    p_+ = n_+/(n_+ + n_-),
    p_- = n_-/(n_+ + n_-).

The logarithmic completion count telescopes:

    log M = sum_{k=1}^r log(N_{k-1}/N_k).

Taking expectation along the uniformly sampled successful leaf yields the
exact chain-rule identity

    log M = sum_{k=1}^r E h(p_k),

where `h(p) = -p log p -(1-p) log(1-p)` is binary entropy and `p_k` is the
conditional sign probability at the revealed prefix.

This is the correct quantity for the uniform-bound problem. It detects the
absolute number M of successful atoms rather than merely the rare-event
probability M/2^r.

## 2. Backward profiles and renormalization

Let the unrevealed tail-angle multiset after time k be represented by the
integer-valued completion profile

    F_k(x) = number of tail sign assignments whose tail angle lies in I-x.

For a squarefree prime block of angle theta_k,

    F_{k-1}(x) = F_k(x+theta_k) + F_k(x-theta_k).

After normalization by total mass, this is the same backward convolution
operator as in the Doob-transform formulation. The ratio

    F_k(x+theta_k) / (F_k(x+theta_k)+F_k(x-theta_k))

is exactly the conditioned transition probability p_k.

Thus the local branching entropy of the reveal tree is a local functional of
the renormalization profile. Large M means that a total amount log M of
branching entropy survives under repeated convolution and rescaling.

## 3. Every genuine branch requires a repair

At a prefix where both signs of prime p have successful descendants, choose
one successful completion from each child. Their final angles both lie in I,
so their difference satisfies

    2 theta_p + sum_{q in T} c_q theta_q = O(delta) mod 2 pi,

where `c_q in { -2,0,2 }` records the difference of the two tail assignments
and T is contained in the unrevealed coordinates.

Equivalently, the current prime flip is repaired by a tail support T. After
halving coefficients, this gives a Gaussian quotient

    U = (pi_p / conjugate(pi_p)) * product_{q in T}
        (pi_q / conjugate(pi_q))^{c_q/2}

with `|arg U| = O(delta)`.

The elementary Gaussian separation estimate implies that every nontrivial
repair relation has total conductor width at least

    (1/2) log N - O_C(1).

Therefore each information-bearing branching node carries a repair witness of
half-conductor scale.

## 4. Minimal repairs and the branch-repair hypergraph

For every branching node v, choose a minimal tail support T_v that repairs the
revealed sign split. Associate the hyperedge

    R_v = {current coordinate} union T_v.

The successful-leaf process gives probability weight to branching nodes. The
total entropy weight of all nodes is exactly log M.

A uniform bound would follow from an amortized packing theorem of the form

    sum_v Prob(reach v) h(p_v) <= O_C(1),

proved by charging each entropy contribution to the conductor mass of its
minimal repair hyperedge with bounded total multiplicity.

The half-conductor lower bound alone is insufficient because repair supports
can overlap heavily. The correct combinatorial object is therefore the
weighted hypergraph of minimal repairs.

## 5. Matching-versus-core strategy

The proof pattern suggested by the mathematical walkthroughs is to separate
independent witnesses from a concentrated obstruction.

If the repair hypergraph contains many essentially disjoint entropy-bearing
hyperedges, then each costs about half the total conductor. Only constantly
many can be disjoint.

If instead the repair hyperedges have a small transversal/core, then most
branching information must be routed through a bounded collection of prime
coordinates. Conditioning on those core coordinates leaves only bounded
residual entropy unless the core itself supports a new endpoint cluster. The
latter should permit descent or a bounded-dimensional arithmetic argument.

The missing statement must be quantitative and probability-weighted, because
one should charge branching nodes according to the probability that the
successful-leaf process reaches them. A purely unweighted sunflower or
matching theorem loses the entropy distribution on the tree.

## 6. Relation to the prime-reveal Doob martingale

There are now two exact entropy identities:

1. Under the original product law conditioned on landing in I,

       log(2^r/M)

   is the sum of relative-entropy costs of biasing the independent prime
   reveals.

2. Under the uniform law on the M successful leaves,

       log M

   is the sum of conditional branching entropies.

Adding the identities gives `r log 2`, the total entropy of the full
orientation cube. The first measures selection cost; the second measures
surviving multiplicity. The endpoint theorem concerns the second.

This duality should be preserved by any renormalization potential.

## 7. A candidate multiscale potential

For a completion profile F and a scale parameter u>0, define the resolvent
coordinate

    Psi_u(F)(x) = F(x)/(F(x)+u).

A scale-integrated potential can combine:

- local branching entropy between the two translated children;
- overlap of their resolvent profiles;
- a rapidly decaying family of enlarged target intervals.

The enlarged intervals serve as a remote correction: a repair that escapes
the endpoint cell must still register at a larger scale. The desired one-step
inequality is

    expected potential after reveal
      <= potential before reveal
         - c * branching entropy
         + arithmetic repair error.

Summing over reveals would bound log M if the arithmetic errors have bounded
total mass.

## 8. Concrete obstruction still to resolve

The remaining obstruction is heavy overlap of minimal repair supports. A
large decision tree can in principle route many branching events through the
same pool of future primes. Abstract binary systems can do this.

The Gaussian arithmetic must show that repeatedly using the same repair core
forces one of the following concrete outcomes:

- an exact multiplicative relation among Gaussian prime orientations;
- a smaller-conductor endpoint instance after conditioning on the core;
- or a quantitative loss in the multiscale survival profile.

Unlike the earlier approximate-grid statement, this formulation is sensitive
to the absolute cluster size and identifies the object that must be controlled:
the probability-weighted hypergraph of branch repairs.

## 9. Immediate next theorem

The next theorem to prove is a probability-weighted matching/core lemma tailored
to the reveal tree:

> If the total branching entropy is L, then either the minimal repair
> hypergraph contains Omega(L) entropy weight supported on essentially
> disjoint repair sets, or there is a conductor-small core of coordinates
> conditioning on which the residual successful-leaf entropy drops by a fixed
> fraction.

The disjoint branch is ruled out by the half-conductor arithmetic cost. An
iterable core branch would give a genuine renormalization contraction.

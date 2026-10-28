# Outreach drafts (not sent)

Three short drafts for putting the open questions to experts. Fill in the bracketed fields. Attach
`cg_questions.tex` (compiled) or link to `docs/claude_review/expert_note/expert_note.md` on the branch
`claude/upbeat-einstein-p0peju` of `craigxchen/Jarnik`. Nothing here has been sent.

---

## 1. Approximation of rational points on toric varieties (Question B)

**Subject:** Canonical boundedness at the identity of a twisted toric 4-fold, and lattice points on short arcs

Dear [Name],

I am working on the Cilleruelo–Granville question of whether an arc of length $C\sqrt R$ on the circle
$x^2+y^2=R^2$ contains a bounded number of lattice points, uniformly in $R$. A long, refereed computation
reduces the uniform bound to a single statement in the language of your work on approximation constants.

Let $X'$ be the Gorenstein toric Fano 4-fold of the $A_4$ root polytope, twisted by $-1$ over $\mathbb Q(i)$, so
that $X'(\mathbb Q)=T_5(\mathbb Q)$ for the anisotropic norm-one torus $T_5$. The question is whether $X'$ is
*canonically bounded* at the identity $e$, in your sense: does every Zariski-dense sequence of rational points
tending to $e$ satisfy $\alpha(-K)\ge 4$? It would be enough to show that the essential constant exceeds $15/4$,
even only along sequences whose 15 toric heights are asymptotically equal. A yes implies the uniform bound.

Some facts we proved along the way may be useful to you:
* $\varepsilon_e(-K)=2$;
* $\tau_e(-K)=10/3$ exactly, for both the split and the twisted form, so every section or Liouville argument stops
  there;
* the split twin of the question implies Erdős Problem #887;
* Galois-aware pointwise arguments give no gain over the split ones.

We found no unconditional theorem of this kind, with $\alpha^{\rm ess}>\tau_x$ at a rational point, for any
variety. We would be very grateful to know whether you think this is known, tractable, or a natural test case.

A short note with precise statements and pointers is attached.

Best regards,
Craig Chen

---

## 2. Zero estimates and the determinant method (Questions A/route 5)

**Subject:** A zero estimate at one rational point of a torus, from a lattice-point problem

Dear [Name],

The determinant method with forced divisibility proves that every large cluster of lattice points on a short
arc of $x^2+y^2=R^2$ satisfies an exact algebraic relation. It is LP-supercritical in the relevant regime.
Converting this into a growth bound $M=O(\log\log R)$ needs one missing lemma: a zero estimate at the actual point
$Z\in T_M(\mathbb Q)$ of the cluster torus, whose coordinates are ratios of Gaussian integers of one norm.

We showed that Philippon-type zero estimates fed with the certified data are met by the trivial subgroup.
The method certifies vanishing only along progressions of length 1–4, against roughly the dimension needed.
"Fake" configurations with the same combinatorics also satisfy every averaged version.

Do you know of a zero or multiplicity estimate that uses the arithmetic of a single rational point in this
way, or would you regard the statement as out of reach? A short note with the precise statement is attached.

Best regards,
Craig Chen

---

## 3. Divisors in short intervals (split calibration)

**Subject:** A five-divisor construction for Erdős #887 and a geometric reformulation

Dear [Name],

In work on lattice points on short arcs of circles, we found that a "split" analogue of our key open
question is a Zariski-closure form of Erdős Problem #887. The analogue concerns five divisors of $n$ within a
factor $1+\mathrm{lcm}^{-(1+2\delta)/4}$ of each other. We prove that it implies #887 in its absolute form, and
part of Ruzsa's #886.

The 2026 five-divisor construction ($n=(x^2-1)(x^2-4)(x^2-16)(x^2-49)$, $z^2=56x^2+65$) turns out to be a
rational curve at exactly the critical (Liouville) exponent $15/4$. It lies in a 2-parameter family that sweeps a
hypersurface. A finiteness statement for such "Pell-balanced" families would be the split analogue of what the
circle problem needs.

Do you know whether the Liouville-extremal five-divisor patterns lie in finitely many such families? Details are
in the attached note.

Best regards,
Craig Chen

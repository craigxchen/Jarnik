# Round 15 (interrupted): Rig^SL, the single-angle truncation of P3-all*

**Status: incomplete.** The workflow stopped when the account hit its weekly usage limit. What finished:

* `rigSLtheory.md` (theory; restricted theorems, **unrefereed**: its referee was cut off, see the partial
  `referee_rigSLtheory.md`). Main claims: closed form of the per-monomial box criterion (Lemma J0); the
  C/A0/A1 classification of affordable extra collisions at one modulus, which explains why q = 83 and
  q = 163 are rigid once one angle vector is required; angle pinning by additive triples (Theorem B), which
  turns one-modulus SL rigidity into an integer criterion; Theorem S, every actual cluster (and every
  P3-all* member with X >= 2M) has r > M/4 distinct split primes; and an outcome table showing that no
  outcome of Rig^SL changes the growth rate of actual circles by itself (rigidity routes back to the
  genuine-prime inputs H^gen + L1; full strictness only makes P3-all* exact).
* `rigSL.md` (the exact test on loose prime sets) is a **partial draft**: the agent was stopped before
  finishing, and nothing in it is refereed.

No growth improvement and no uniform bound. Only files under 1 MiB are committed; caches (2.5 GB of
generated tables in `rigSL_checks/base/cache/`), pickles, arrays and compiled binaries are not. The
theory agent's result is cached in the workflow journal, so the round can be resumed.

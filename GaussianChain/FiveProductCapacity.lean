import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

namespace GaussianChain
namespace FiveProductCapacity

/-!
# Five-product pair capacity

This file records the exact scalar obstruction to obtaining fourth-order
angular decay by multiplying reduced pair determinants of five signed Gaussian
products.

For each of the ten `2 | 3` sign cuts, a one-use conductor budget bounds the
sum of the six determinant multiplicities crossing that cut by one.  Averaging
the ten inequalities counts each pair six times, so the total determinant
multiplicity is at most `5 / 3`.

The theorem is deliberately abstract.  It does not construct Gaussian
products or assert that an arithmetic configuration satisfies the ten cut
hypotheses.
-/

/--
The ten balanced signed cuts on five vertices force total pair multiplicity at
most `5 / 3`.  The conclusion is written without division as
`3 * total ≤ 5`.
-/
theorem signed_pair_capacity_le
    {x01 x02 x03 x04 x12 x13 x14 x23 x24 x34 : ℝ}
    (h01 : x02 + x03 + x04 + x12 + x13 + x14 ≤ 1)
    (h02 : x01 + x03 + x04 + x12 + x23 + x24 ≤ 1)
    (h03 : x01 + x02 + x04 + x13 + x23 + x34 ≤ 1)
    (h04 : x01 + x02 + x03 + x14 + x24 + x34 ≤ 1)
    (h12 : x01 + x02 + x13 + x14 + x23 + x24 ≤ 1)
    (h13 : x01 + x03 + x12 + x14 + x23 + x34 ≤ 1)
    (h14 : x01 + x04 + x12 + x13 + x24 + x34 ≤ 1)
    (h23 : x02 + x03 + x12 + x13 + x24 + x34 ≤ 1)
    (h24 : x02 + x04 + x12 + x14 + x23 + x34 ≤ 1)
    (h34 : x03 + x04 + x13 + x14 + x23 + x24 ≤ 1) :
    3 * (x01 + x02 + x03 + x04 + x12 + x13 + x14 + x23 + x24 + x34) ≤ 5 := by
  linarith

/-- Uniform multiplicity `1 / 6` on every pair attains the signed capacity. -/
theorem signed_uniform_attains_capacity :
    3 * (10 * ((1 : ℝ) / 6)) = 5 := by
  norm_num

/--
For rooted binary occurrence patterns a crossing costs one half of a conductor
use.  The analogous ten cut inequalities therefore have right side two and
give the exact fractional upper bound `10 / 3`.
-/
theorem rooted_binary_pair_capacity_le
    {x01 x02 x03 x04 x12 x13 x14 x23 x24 x34 : ℝ}
    (h01 : x02 + x03 + x04 + x12 + x13 + x14 ≤ 2)
    (h02 : x01 + x03 + x04 + x12 + x23 + x24 ≤ 2)
    (h03 : x01 + x02 + x04 + x13 + x23 + x34 ≤ 2)
    (h04 : x01 + x02 + x03 + x14 + x24 + x34 ≤ 2)
    (h12 : x01 + x02 + x13 + x14 + x23 + x24 ≤ 2)
    (h13 : x01 + x03 + x12 + x14 + x23 + x34 ≤ 2)
    (h14 : x01 + x04 + x12 + x13 + x24 + x34 ≤ 2)
    (h23 : x02 + x03 + x12 + x13 + x24 + x34 ≤ 2)
    (h24 : x02 + x04 + x12 + x14 + x23 + x34 ≤ 2)
    (h34 : x03 + x04 + x13 + x14 + x23 + x24 ≤ 2) :
    3 * (x01 + x02 + x03 + x04 + x12 + x13 + x14 + x23 + x24 + x34) ≤ 10 := by
  linarith

/-- Uniform multiplicity `1 / 3` on every pair attains the binary capacity. -/
theorem rooted_binary_uniform_attains_capacity :
    3 * (10 * ((1 : ℝ) / 3)) = 10 := by
  norm_num

end FiveProductCapacity
end GaussianChain

"""Finite certificate for the explicit twelve-row pencil propagation note.

All ranks are exact over F_101. The companion proof supplies the
characteristic-zero upper bounds and the all-parameter conclusion.
"""

import numpy as np

import check_generic_terminal_coherent_syzygies as base


def main():
    directions = list(base.LABELS)
    determinants = base.triple_determinants(base.test_vectors())
    constant_columns = []
    for six in base.SIX_SETS:
        if 0 in six:
            continue
        complement = tuple(i for i in base.LABELS if i not in six)
        matching = base.matching_values(six, directions)
        local_relation = base.gradient_values(six, determinants) @ matching % base.PRIME
        constant_columns.append(
            local_relation[:, None] * base.gradient_values(complement, determinants)
            % base.PRIME
        )
    constant = np.concatenate(constant_columns, axis=1)
    assert constant.shape == (462, 2310)
    assert base.modular_rank(constant.T) == 286

    full = base.local_generator_evaluations(directions)
    assert full.shape == (462, 4620)
    assert base.modular_rank(full.T) == 341

    at_zero = base.coherent_matrix(directions)
    changed = list(directions)
    changed[0] = 13
    slope = (base.coherent_matrix(changed) - at_zero) * pow(13, -1, base.PRIME)
    slope %= base.PRIME
    assert base.modular_rank(at_zero) == 121
    assert base.modular_rank(np.concatenate((at_zero, slope))) == 176
    print("F_101 ranks: constant local family 286, full local family 341,")
    print("coherent fiber 121, stacked pencil 176.")


if __name__ == "__main__":
    main()

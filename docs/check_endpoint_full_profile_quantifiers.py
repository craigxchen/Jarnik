"""Exact threshold/gcd certificate for the sharpened pair residual bound."""
from itertools import combinations, product


def main():
    tuples = pairs = 0
    for k in range(3, 6):
        for e in range(1, 9):
            for allocation in product(range(e + 1), repeat=k):
                interior = [min(a, e - a) for a in allocation]
                radius = max(interior)
                retained = e - 2 * radius
                labels = [2 * a >= e for a in allocation]
                differences = [a - allocation[0] for a in allocation]
                for i in range(1, k):
                    core = retained if labels[i] != labels[0] else 0
                    assert 0 <= abs(differences[i]) - core <= 2 * radius
                for i, j in combinations(range(1, k), 2):
                    di, dj = differences[i], differences[j]
                    gcd_depth = min(abs(di), abs(dj)) if di * dj > 0 else 0
                    assert 2 * gcd_depth == (abs(di) + abs(dj)
                                              - abs(allocation[i] - allocation[j]))
                    common_core = (retained if labels[i] != labels[0]
                                   and labels[j] != labels[0] else 0)
                    original_mass = core_mass = 0
                    for threshold in range(1, e + 1):
                        anchor = allocation[0] >= threshold
                        both = ((allocation[i] >= threshold) != anchor
                                and (allocation[j] >= threshold) != anchor)
                        original_mass += both
                        if radius < threshold <= e - radius:
                            core_mass += both
                    assert original_mass == gcd_depth
                    assert core_mass == common_core
                    assert 0 <= gcd_depth - common_core <= 2 * radius
                    pairs += 1
                tuples += 1
    print(f'Gcd and retained-layer identities verified on {tuples} tuples, {pairs} pairs.')


if __name__ == '__main__':
    main()

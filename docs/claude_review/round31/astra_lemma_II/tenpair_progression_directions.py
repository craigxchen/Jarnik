"""Exact check for five affine points whose ten directions have slopes 0..9."""
from fractions import Fraction
from itertools import combinations, permutations

# Translate and scale the two ends of the unique horizontal edge to (0, 0), (1, 0).
# Each other point has slopes a and b from these two anchors.
slopes = frozenset(range(1, 10))
possible = []
for a, b in permutations(sorted(slopes), 2):
    x = Fraction(b, b-a)
    possible.append((a, b, x, a*x))

tested = 0
solutions = 0
for triangle in combinations(possible, 3):
    used = {direction for point in triangle for direction in point[:2]}
    if len(used) != 6:
        continue
    tested += 1
    for first, second in combinations(triangle, 2):
        _, _, x, y = first
        _, _, u, v = second
        if x == u:
            break
        direction = (v-y)/(u-x)
        if direction not in slopes or direction in used:
            break
        used.add(direction)
    else:
        solutions += 1

print("Tested", tested, "triangles with distinct anchor slopes; solutions =", solutions)
assert tested == 10080
assert solutions == 0

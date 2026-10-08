"""Exact first-moment barycentric witness for the even-q2 compressed automaton."""

from check_critical_moment_compressed_automaton import Automaton


BASE = (128, 64, 206, 205, 35, 19, 240, 8, 4, 60, 195)

# Each pair is an active path terminal -> axial height and an active return.
EXCURSIONS = {
    (1, -1): ((195, 60, 253), (2, 62, 193)),
    (1, 1): ((193, 2, 62), (60, 253, 195)),
    (2, -1): ((199, 56, 251), (4, 60, 195)),
    (2, 1): ((195, 4, 60), (56, 251, 199)),
    (3, -1): ((207, 48, 247), (8, 56, 199)),
    (3, 1): ((199, 8, 56), (48, 247, 207)),
    (4, -1): ((223, 32, 239), (16, 48, 207)),
    (4, 1): ((207, 16, 48), (32, 239, 223)),
    (5, -1): ((223,), (32, 32, 223)),
    (5, 1): ((223, 32, 32), (223,)),
    (6, -1): ((64, 191, 191), (64,)),
    (6, 1): ((64,), (191, 191, 64)),
    (7, -1): ((192, 63, 127), (128, 191, 64)),
    (7, 1): ((64, 128, 191), (63, 127, 192)),
}


def run():
    machine = Automaton(0, 2)
    delta = machine.delta
    state = machine.initial
    prefix = []
    columns = []

    def step(mask):
        nonlocal state
        edges = {m: nxt for m, nxt, _ in machine.transitions(state)}
        assert mask in edges and mask not in (0, 255)
        prefix.append(state[0])
        columns.append(mask)
        state = edges[mask]

    for mask in BASE:
        step(mask)
    assert state == machine.terminal

    for (i, sign), (out, back) in EXCURSIONS.items():
        assert len(out) <= 3 and len(back) <= 3
        for mask in out:
            step(mask)
        expected = tuple(delta[j] + sign * (i == j) for j in range(8))
        assert state[0] == expected
        for mask in back:
            step(mask)
        assert state == machine.terminal

    n = len(columns)
    assert n == 87 and state[0] == delta
    assert all(abs(r[i] - delta[i]) <= 1 for r in prefix for i in range(8))
    numerator = [2] * n  # denominator 16*n
    v = [sum(delta[i] - r[i] for r in prefix) for i in range(8)]
    for i in range(1, 8):
        assert abs(v[i]) <= n
        for sign in (-1, 1):
            expected = tuple(delta[j] + sign * (i == j) for j in range(8))
            t = prefix.index(expected)
            numerator[t] += n + sign * v[i]

    assert min(numerator) >= 2
    assert sum(numerator) == 16 * n
    assert all(sum(numerator[t] * prefix[t][i] for t in range(n))
               == 16 * n * delta[i] for i in range(8))

    magnitudes = []
    running = 0
    for increment in numerator:
        running += increment
        magnitudes.append(running)
    assert magnitudes[0] >= 2 and magnitudes[-1] == 16 * n
    assert all(x < y for x, y in zip(magnitudes, magnitudes[1:]))
    moments = [sum((1 if mask >> i & 1 else -1) * a
                   for mask, a in zip(columns, magnitudes)) for i in range(8)]
    assert len(set(moments)) == 1
    print(f"PASS: {n} active columns, 14 axial heights, integer magnitudes 2..{16*n}, common first moment {moments[0]}.")


if __name__ == "__main__":
    run()

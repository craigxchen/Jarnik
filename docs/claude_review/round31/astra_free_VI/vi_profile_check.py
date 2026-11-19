"""Exact weighted packing obstruction for the ordinary Paley double core."""
from math import comb

from ii_paley import double_core


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            leading = row.bit_length()
            if leading not in pivots:
                pivots[leading] = row
                break
            row ^= pivots[leading]
    return len(pivots)


def check(q):
    first, second = double_core(q)
    columns = [[v for v in col[1:]] for col in first + second]
    assert all([x ^ y for x, y in zip(a, b)] ==
               [int(row == index) for row in range(q)]
               for index, (a, b) in enumerate(zip(columns, columns[q:])))
    anchor = [column[0] for column in columns]
    rows = [sum((column[i] ^ bit) << j
                for j, (column, bit) in enumerate(zip(columns, anchor)))
            for i in range(1, q)]
    rank = binary_rank(rows)
    n = len(columns)
    j = n // 8 - 1
    # Markov gives ball size >= (1-8*j/n)*binom(n,j), for all weights.
    numerator = (n - 8*j) * comb(n, j)
    denominator = n * (1 << (n - rank))
    print(f"rows={q}, columns={n}, rank={rank}, subset_size={j}, "
          f"packed_volume/space > 1: {numerator > denominator}")
    assert rank == q - 1 and numerator > denominator


def check_two_cores(d, modulus):
    # modulus is a supplied irreducible binary polynomial of degree d.
    m = 1 << d
    def multiply(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            if a & m:
                a ^= modulus
            b >>= 1
        return result

    cubes = [multiply(multiply(x, x), x) for x in range(m)]
    assert len(set(cubes)) == m
    assert binary_rank([x | (y << d) for x, y in enumerate(cubes)]) == 2*d
    expressions = [
        {v for v in range(1, m) if (v & x).bit_count() % 2}
        for x in range(1, m)
    ]
    assert all(len(s) == m//2 for s in expressions)
    assert all(len(s & t) == m//4 for i, s in enumerate(expressions)
               for t in expressions[:i])
    print(f"two Walsh cores: m={m}, rank={2*d}, first allowed even degree={4*d-2}")


if __name__ == "__main__":
    for prime in (167, 191, 251):
        check(prime)
    check_two_cores(3, 0b1011)
    check_two_cores(5, 0b100101)

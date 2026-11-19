"""Independently audit the witnesses in vii_squarefree_table.out using python-flint."""
from pathlib import Path
from re import match

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

prime = fmpz_mod_poly_ctx(3)
field = fq_default_ctx(modulus=prime([1, 0, 1]), var="i")
ring = fq_default_poly_ctx(field)


def decode(index):
    coefficients = []
    while index:
        index, digit = divmod(index, 9)
        b, a = divmod(digit, 3)
        coefficients.append(a + field.gen() * b)
    return ring(coefficients)


def audit(line):
    record = match(r"D=(\d+) M=(\d+) base9_marks=([\d,]+)", line)
    degree, size = map(int, record.groups()[:2])
    points = [decode(int(s)) for s in record.group(3).split(",") if s]
    assert len(points) == size
    norms = [z * ring([c**3 for c in z.coeffs()]) for z in points]
    parent = norms[0]
    assert all(z.degree() == degree for z in points)
    assert all(norm == parent for norm in norms)
    assert parent.degree() == 2*degree
    assert all(c == c**3 for c in parent.coeffs())
    assert parent.gcd(parent.derivative()) == 1
    gaps = [(u-v).degree() for j, u in enumerate(points) for v in points[:j]]
    assert min(gaps) >= 0 and 2*max(gaps) <= degree
    print(f"D={degree} marks={size} max_gap={max(gaps)} squarefree=yes")


if __name__ == "__main__":
    for row in Path(__file__).with_name("vii_squarefree_table.out").read_text().splitlines():
        audit(row)

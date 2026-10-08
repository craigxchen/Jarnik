"""Exact SymPy irreducibility and real-root certificates for two witnesses.

Run with SymPy available, or in the existing cached environment:
UV_CACHE_DIR=/tmp/uv-cache uv run --offline --with sympy python \
    docs/check_six_point_branch_arithmetic.py
"""

import sympy as sp

from check_six_point_isotropic_circle_cover import probe


def main():
    parameter = sp.Symbol('t')
    for parameters in ((0,1,2,3,4,5), (0,1,2,4,7,11)):
        weights, _, _, _, coefficients = probe(parameters)
        polynomial = sp.Poly.from_dict(
            {(i,): sp.Rational(c.numerator, c.denominator)
             for i, c in enumerate(coefficients)}, parameter, domain=sp.QQ,
        )
        assert polynomial.degree() == 12
        _, factors = polynomial.factor_list()
        assert len(factors) == 1
        factor, multiplicity = factors[0]
        assert factor.degree() == 12 and multiplicity == 1
        real_roots = polynomial.count_roots(-sp.oo, sp.oo)
        assert real_roots == 2
        print(f'PASS: weights {weights}: irreducible degree 12 over Q; '
              'exactly two real roots.')


if __name__ == '__main__':
    main()

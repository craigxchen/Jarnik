"""Exact algebra for an entire one-group family; not finite-template feasibility."""

from fractions import Fraction as F


ZERO = (0, 0, 0, 0, 0)  # Exponents of B,C,z,t_i,t_j.


def variable(index):
    exponent = list(ZERO)
    exponent[index] = 1
    return {tuple(exponent): F(1)}


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            out[exponent] = out.get(exponent, F(0)) + coefficient
    return {e: c for e, c in out.items() if c}


def negate(polynomial):
    return {e: -c for e, c in polynomial.items()}


def multiply(*polynomials):
    out = {ZERO: F(1)}
    for polynomial in polynomials:
        updated = {}
        for left, a in out.items():
            for right, b in polynomial.items():
                exponent = tuple(x + y for x, y in zip(left, right))
                updated[exponent] = updated.get(exponent, F(0)) + a * b
        out = {e: c for e, c in updated.items() if c}
    return out


def product(values):
    result = F(1)
    for value in values:
        result *= value
    return result


def run():
    b, c, z, ti, tj = [variable(i) for i in range(5)]
    ai = add(b, multiply(ti, z, c))
    ai_reflected = add(c, negate(multiply(ti, z, b)))
    aj = add(b, multiply(tj, z, c))
    aj_reflected = add(c, negate(multiply(tj, z, b)))
    squares = add(multiply(b, b), multiply(c, c))
    assert add(multiply(ai_reflected, c), multiply(ai, b)) == squares
    assert add(multiply(ai_reflected, aj), negate(multiply(ai, aj_reflected))) == (
        multiply(add(tj, negate(ti)), z, squares))
    determinant = add({ZERO: F(1)}, multiply(ti, tj, z, z))
    assert determinant == add({ZERO: F(1)}, negate(multiply(ti, z, negate(multiply(tj, z)))))

    # Rational values of B,C satisfying B^2+C^2=2; all identities are
    # algebraic in these variables after the trigonometric identity is used.
    b_value, c_value, z_value = F(7, 5), F(1, 5), F(2, 7)
    assert b_value ** 2 + c_value ** 2 == 2
    for m in [1, 3, 4, 7, 12]:
        parameters = [F(i + 1) for i in range(m)]
        direct = [b_value + t * z_value * c_value for t in parameters]
        reflected = [c_value - t * z_value * b_value for t in parameters]
        common_norm = b_value * c_value * product(
            left * right for left, right in zip(direct, reflected))
        for i in range(m):
            row = reflected[i] * c_value * product(
                direct[j] for j in range(m) if j != i)
            reflected_row = direct[i] * b_value * product(
                reflected[j] for j in range(m) if j != i)
            assert row * reflected_row == common_norm
            assert reflected[i] * c_value + direct[i] * b_value == 2
            for j in range(m):
                assert reflected[i] * direct[j] - direct[i] * reflected[j] == (
                    2 * (parameters[j] - parameters[i]) * z_value)
        last_row = b_value * product(direct)
        last_reflected = c_value * product(reflected)
        assert last_row * last_reflected == common_norm
        total_derivative = sum(1 + t for t in parameters)
        beta = [3 + 2 * t - total_derivative for t in parameters]
        beta_last = -1 - total_derivative
        assert len(set(beta + [beta_last])) == m + 1
        assert all(value - beta_last == 4 + 2 * t
                   for value, t in zip(beta, parameters))
        assert all(1 + ti * tj * z_value ** 2 > 0
                   for ti in parameters for tj in parameters)
    print("PASS: symbolic cross/within residuals and opposite-root determinant")
    print("PASS: reflection products and distinct derivative labels for m=1,3,4,7,12")


if __name__ == "__main__":
    run()

"""Constants of Theorem 5.1 (50-digit Decimal)."""
from decimal import Decimal as D, getcontext
getcontext().prec = 50
PI = D('3.14159265358979323846264338327950288419716939937510')
pre = 4 * (PI + 1) * (D(1) / 16).exp()
print('4(pi+1)e^(1/16) =', pre, '<= 17.64:', pre <= D('17.64'))
x = D('0.1')
s = x * x * (2 - x) / (1 - x) ** 2
print('sum_{n>=2} n x^n at x=0.1 =', s, ' ratio to x^2 =', s / (x * x), '<= 2.346:', s / (x * x) <= D('2.346'))
print('17.64 * 2.346 / 100 =', D('17.64') * D('2.346') / 100, '<= 0.414')
# arcsin(t)/t at t = 5^(-1/2) (arcsin(t)/t is increasing)
t = D(5).sqrt() ** -1
# arcsin via series / Newton on sin
y = t
for _ in range(60):
    # Newton: f(y) = sin y - t
    s_, c_ = D(0), D(0)
    term = y; k = 1; s_ = D(0)
    for n in range(40):
        s_ += term
        term *= -y * y / ((2 * n + 2) * (2 * n + 3))
    term = D(1); c_ = D(0)
    for n in range(40):
        c_ += term
        term *= -y * y / ((2 * n + 1) * (2 * n + 2))
    y = y - (s_ - t) / c_
print('arcsin(5^-1/2)/5^-1/2 =', y / t, '<= 1.0368:', y / t <= D('1.0368'))
print('log 160000 =', D(160000).ln(), ' + log 2 =', D(160000).ln() + D(2).ln())

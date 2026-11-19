# F_27 = F_3[a]/(a^3 - a - 1); elements encoded as ints 0..26 = c0 + 3 c1 + 9 c2
import itertools
def dec(x): return [x % 3, (x // 3) % 3, x // 9]
def enc(c): return c[0] % 3 + 3 * (c[1] % 3) + 9 * (c[2] % 3)
def fadd(x, y): a, b = dec(x), dec(y); return enc([a[i] + b[i] for i in range(3)])
def fneg(x): return enc([-c for c in dec(x)])
def fmul(x, y):
    a, b = dec(x), dec(y)
    r = [0] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] += a[i] * b[j]
    # a^4 = a^2 + a, a^3 = a + 1
    r[2] += r[4]; r[1] += r[4]
    r[1] += r[3]; r[0] += r[3]
    return enc(r[:3])
def fpow(x, e):
    r = 1
    while e:
        if e & 1: r = fmul(r, x)
        x = fmul(x, x); e >>= 1
    return r
def frob(x): return fpow(x, 3)
def chi(x):
    if x == 0: return 0
    v = fpow(x, 13)
    assert v in (1, 2)
    return 1 if v == 1 else -1
def from_f3(c): return c % 3
# sanity: field
assert all(fmul(x, fpow(x, 25)) == 1 for x in range(1, 27))
assert chi(2) == -1

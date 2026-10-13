# Exhaustive: (tau3,tau4,tau5) >= 0 with expected dim 15 - t3 - 2t4 - 3t5 >= 3 and (B+8K) = 6t3+11t4+15t5-70 >= 0
sol = [(a, b, c, 6*a + 11*b + 15*c - 70) for a in range(20) for b in range(10) for c in range(6)
       if 15 - a - 2*b - 3*c >= 3 and 6*a + 11*b + 15*c - 70 >= 0]
print(sol)

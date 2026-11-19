from sympy import Poly, symbols
t=symbols('t')
def P(s): return Poly([int(c) for c in s],t,modulus=3)
her=Poly(t**24+t**20+2*t**18+t**16+2*t**10+2*t**2+2,t,modulus=3)
ws=["1020000020200020200000201","1002121010020212012021002","1201021110200002001121202"]
def orbit(p):
    out=set()
    for a in range(3):
        for s in (1,2):
            q=Poly(p.as_expr().subs(t,s*t+a),t,modulus=3)
            q=q*pow(int(q.LC())%3,1,3)  # make monic: LC in {1,2}, inverse = itself
            out.add(tuple(int(c)%3 for c in q.all_coeffs()))
    return out
herc=tuple(int(c)%3 for c in her.all_coeffs())
for w in ws:
    p=P(w); o=orbit(p)
    print(w, "orbit size",len(o),"contains Astra N:",herc in o, "even:", all(int(c)%3==0 for c in p.all_coeffs()[1::2]))
    print("   factors:",[(str(f.as_expr()).replace('**','^'),e) for f,e in p.factor_list()[1]])

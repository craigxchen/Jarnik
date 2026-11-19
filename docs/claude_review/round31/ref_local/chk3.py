exec(open('chk2.py').read().split('tab={}')[0])
def v3(n):
    k=0
    while n%3==0: n//=3;k+=1
    return k,n
print("twisted at 3: m mod 9 -> (A mod 81, B mod 81)")
for r in range(9): print(r, sorted({(A(x)%81,B(x)%81) for x in range(r,3**6,9)})[:4])
print("split at 2: odd m mod 8 -> (-A mod 256, -B mod 256)")
for r in (1,3,5,7): print(r,{((-A(x))%256,(-B(x))%256) for x in range(r,4096,8)})

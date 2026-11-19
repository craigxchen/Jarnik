import pickle, sympy as sp
exec(open('canon_helper.py').read())
exec(open('filter_ref.py').read().split('out=[]')[0].replace("reps=pickle.load(open('reps6.pkl','rb'))",""))
R8=pickle.load(open('reps8_rank4.pkl','rb'))
A=[(15,51,60,85,170,195,204,240),(15,51,85,102,153,170,204,240),(15,51,85,105,150,170,204,240)]
cA=[canon_rows(list(t)) for t in A]
print('mine',R8); print('astra canon',cA)
for S in R8:
    rk,M,ns=kernel(S); H=sp.Matrix.hstack(*ns)
    print(S,'rank',rk,'H rows',[tuple(H.row(j)) for j in range(8)])

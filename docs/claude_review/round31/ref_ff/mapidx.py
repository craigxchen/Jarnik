import pickle
exec(open('canon_helper.py').read())
J=pickle.load(open('/home/user/Jarnik/docs/claude_review/round30/judge_checks/classes_D8_c0.pkl','rb'))
surv=pickle.load(open('survivors.pkl','rb'))
pos={s:i for i,s in enumerate(surv)}
m=[pos[canon_signs(v)] for v in J.values()]
print('judge idx -> my idx',m)
for g in (7,0,3,5):
    r=pickle.load(open(f'nonlin_res_{g}.pkl','rb'))
    print(g,[r[m[i]][0] for i in range(21)])

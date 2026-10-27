import json
for s in ['first4','rand4','bal4','last4']:
    mine=[];last=0
    for l in open('G_%s_K5.err'%s):
        t=l.split()
        if t[0]=='q':
            q=int(t[1]); last=q
            if t[2]=='forced': mine.append(q)
    a=[q for q in json.load(open('Fa_%s_K5.json'%s))['F'] if q<=last]
    print(s,'checked up to',last,'agree' if a==mine else 'DIFFER', sorted(set(a)^set(mine)), 'author max',max(json.load(open('Fa_%s_K5.json'%s))['F']))

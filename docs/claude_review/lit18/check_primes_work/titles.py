import sys, json, re
db = json.load(open(sys.argv[1]))
used = set(open(sys.argv[2]).read().split())
pats = [re.compile(p, re.I|re.S) for p in sys.argv[3:]]
n=0
for k, r in sorted(db.items(), key=lambda kv: kv[1]['first']):
    if not r['cat'].startswith('math.'): continue
    txt = r['title'] + ' ' + r['abstract']
    if all(p.search(txt) for p in pats):
        n+=1
        flag = 'USED' if r['id'] in used else 'new '
        print(f"{flag} {r['id']} [{r['first'][:7]} {r['cat']}] {r['title'][:110]} | {r['authors'][:60]}")
print('matches', n)

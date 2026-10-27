import sys, json, re
db = json.load(open(sys.argv[1]))
pats = [re.compile(p, re.I) for p in sys.argv[2:]]
n=0
for k, r in sorted(db.items(), key=lambda kv: kv[1]['first']):
    txt = r['title'] + ' ' + r['abstract']
    if all(p.search(txt) for p in pats):
        n+=1
        print(f"### {r['id']} [{r['first']}..{r['date']} v{r['ver']} {r['cat']}] {r['title']}\n   AUTH: {r['authors'][:150]}\n   {r['abstract'][:900]}\n")
print('matches', n)

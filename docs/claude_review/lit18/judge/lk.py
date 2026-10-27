import sys, json
db = json.load(open(sys.argv[1]))
for i in sys.argv[2:]:
    r = db.get(i)
    if not r: print('MISSING', i); continue
    print(f"### {r['id']} [{r['first']}..{r['date']} v{r['ver']} {r['cat']}] {r['title']}\n   AUTH: {r['authors'][:200]}\n   {r['abstract'][:1400]}\n")

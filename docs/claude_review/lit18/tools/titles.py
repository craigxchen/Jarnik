import json,re,sys
d=json.load(open(sys.argv[1]))
kw=sys.argv[2:] 
print('total',d.get('total_count'))
seen=set()
for it in d['items']:
    for tm in it.get('text_matches',[]):
        f=tm['fragment']
        items=re.findall(r'<title>(.*?)</title>\s*<link>(.*?)</link>\s*<description>(.*?)(?:</description>|$)',f,re.S)
        if not items:
            print('RAW', it['path'].split('/')[-1][:10], f[:300].replace('\n',' '))
        for t in items:
            if t[1] not in seen:
                seen.add(t[1]); print('*', it['path'].split('/')[-1][:10], t[1].split('/')[-1], t[0][:150]); print('    ', t[2][:600].replace('\n',' '))

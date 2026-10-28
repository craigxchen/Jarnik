import re,subprocess,sys
p='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round17/calib.md'
t=subprocess.run(['python3','tables.py'],capture_output=True,text=True,timeout=120).stdout.strip()
s=open(p).read()
blk='<!--TABLE-->\n'+t+'\n<!--/TABLE-->'
if 'TABLE_MAIN' in s: s=s.replace('TABLE_MAIN',blk)
else: s=re.sub(r'<!--TABLE-->.*?<!--/TABLE-->',lambda m:blk,s,flags=re.S)
open(p,'w').write(s); print('ok',len(t.splitlines()))

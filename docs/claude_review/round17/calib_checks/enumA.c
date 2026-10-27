/* Enumerate the T=1 slab box A_L of sdc_q: k in {-1,0,1}^q, |k(P)| <= L for the 2q patterns
   P_a = {a} u {x : chi(a-x) = 1}, P'_a = {a} u {x : chi(a-x) = -1}.  Writes canonical translation-orbit
   representatives (lexicographically smallest rotation) as int8 rows if outfile given; prints counts.
   usage: enumA q L [outfile] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int q, L, np; static int inc[128][256]; static int ninc[128]; static int rem[129][256];
static int sums[256]; static int smax, supp=0; static signed char k[128]; static long long cnt=0, nfree=0, nfixed=0; static FILE *fo=NULL;
static int chi(int a){ a%=q; if(a<0)a+=q; if(a==0) return 0; long long r=1,b=a; int e=(q-1)/2; while(e){ if(e&1) r=r*b%q; b=b*b%q; e>>=1;} return r==1?1:-1; }
static int pat[256][128];
static void visit(void){
  /* is k the lexicographically smallest rotation? (compare as sequences over -1<0<1) */
  int s, t, minimal=1, fixed=1;
  for(s=1;s<q;s++){ int c=0; for(t=0;t<q;t++){ int a=k[(t+s)%q], b=k[t]; if(a!=b){ c=(a<b)?-1:1; break;} }
    if(c==0){ /* periodic: q prime => constant */ } else fixed=0;
    if(c<0){ minimal=0; break; } }
  cnt++;
  if(!minimal) return;
  if(fixed) nfixed++; else nfree++;
  if(fo) fwrite(k,1,q,fo);
}
static void rec(int x){
  if(x==q){ visit(); return; }
  for(int v=-1; v<=1; v++){
    int ok=1;
    for(int i=0;i<ninc[x];i++) sums[inc[x][i]]+=v;
    int left = smax - supp - (v!=0);
    if(left<0) ok=0;
    if(ok) for(int j=0;j<np;j++){ int r=rem[x+1][j]; if(r>left) r=left; if(abs(sums[j]) - r > L){ ok=0; break; } }
    if(ok && x==q/2 && v==1){ static long long ctr=0; if((++ctr & 1023)==0){ fprintf(stderr,"progress k0..k2=%d %d %d cnt=%lld\n",k[0],k[1],k[2],cnt); } }
    if(ok){ k[x]=v; if(v) supp++; rec(x+1); if(v) supp--; }
    for(int i=0;i<ninc[x];i++) sums[inc[x][i]]-=v;
  }
  k[x]=0;
}
int main(int argc,char**argv){
  q=atoi(argv[1]); L=atoi(argv[2]); if(argc>3) fo=fopen(argv[3],"wb");
  np=2*q;
  /* Parseval: sum_a (chi*k)(a)^2 = q|supp k| - S^2 and |S+k_a|+|chi*k(a)| <= 2L, |S| < 2L  => |supp k| <= smax */
  { int best=0; for(int s=0;s<=q;s++){ int ok=0; for(int S=-(2*L-1);S<=2*L-1;S++){ /* max over k_a pattern: crude bound sum_a (2L)^2 */ if((long)q*s - (long)S*S <= (long)q*4*L*L) ok=1;} if(ok) best=s;} smax=best; }
  fprintf(stderr,"smax=%d\n",smax);
  for(int a=0;a<q;a++) for(int x=0;x<q;x++){ pat[a][x]=(x==a)||chi(a-x)==1; pat[q+a][x]=(x==a)||chi(a-x)==-1; }
  for(int x=0;x<q;x++){ ninc[x]=0; for(int j=0;j<np;j++) if(pat[j][x]) inc[x][ninc[x]++]=j; }
  for(int x=0;x<=q;x++) for(int j=0;j<np;j++){ int c=0; for(int y=x;y<q;y++) c+=pat[j][y]; rem[x][j]=c; }
  rec(0);
  printf("q=%d L=%d |A|=%lld free_orbits=%lld fixed=%lld\n",q,L,cnt,nfree,nfixed);
  if(fo) fclose(fo);
  return 0;
}

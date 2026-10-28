/* Referee: DFS enumeration of A_L(q) = {k in {-1,0,1}^q : |k(P)| <= L for the 2q sdc patterns}. Writes int8 rows. */
#include <stdio.h>
#include <stdlib.h>
int q, L, np; int pat[128][128]; int rem[128][128]; int s[128]; signed char k[128]; long cnt=0; FILE *fo;
int chi(int a){ a%=q; if(a<0)a+=q; if(a==0) return 0; long r=1,b=a; int e=(q-1)/2; while(e){ if(e&1) r=r*b%q; b=b*b%q; e>>=1;} return r==1?1:-1; }
void rec(int x){
  if(x==q){ if(fo) fwrite(k,1,q,fo); cnt++; if((cnt&((1L<<24)-1))==0){fprintf(stderr,"cnt=%ld\n",cnt);} return; }
  for(int v=-1; v<=1; v++){
    int ok=1;
    for(int p=0;p<np;p++){ int ns=s[p]+v*pat[p][x]; if(ns-rem[p][x+1]>L || ns+rem[p][x+1]< -L){ok=0;break;} }
    if(!ok) continue;
    for(int p=0;p<np;p++) s[p]+=v*pat[p][x];
    k[x]=v; rec(x+1);
    for(int p=0;p<np;p++) s[p]-=v*pat[p][x];
  }
}
int main(int argc,char**argv){ q=atoi(argv[1]); L=atoi(argv[2]); fo=(argv[3][0]==0x2d)?NULL:fopen(argv[3],"wb"); np=2*q;
  for(int a=0;a<q;a++) for(int x=0;x<q;x++){ pat[2*a][x]=(x==a)||chi(a-x)==1; pat[2*a+1][x]=(x==a)||chi(a-x)==-1; }
  for(int p=0;p<np;p++){ rem[p][q]=0; for(int x=q-1;x>=0;x--) rem[p][x]=rem[p][x+1]+pat[p][x]; }
  rec(0); if(fo) fclose(fo); printf("q=%d L=%d |A|=%ld\n",q,L,cnt); return 0; }

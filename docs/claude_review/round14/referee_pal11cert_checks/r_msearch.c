/* Referee: generic exhaustive multiset search: m generators from a pool (n<=12 rows, base row 0), counts only. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int ND,dd; int D[80000][12]; int tau[80000]; int P; int Z[5000][12]; int nt[5000]; int *tl[5000];
int cnt[80000]; int m; int cur[16]; long long leaves=0; FILE*out;
void rec(int start,int k){
  if(k==m){ leaves++; for(int i=0;i<m;i++) fprintf(out,"%d ",cur[i]); fprintf(out,"\n"); return; }
  for(int i=start;i<P;i++){ int ok=1; for(int j=0;j<nt[i];j++){int v=tl[i][j]; if(++cnt[v]>tau[v]) ok=0;}
    if(ok){cur[k]=i; rec(i,k+1);} for(int j=0;j<nt[i];j++) cnt[tl[i][j]]--; }
}
int main(int argc,char**argv){
  dd=atoi(argv[1]); FILE*f=fopen(argv[2],"r"); if(fscanf(f,"%d",&ND)!=1) return 1;
  for(int i=0;i<ND;i++){ for(int x=1;x<=dd;x++) if(fscanf(f,"%d",&D[i][x])!=1) return 1; if(fscanf(f,"%d",&tau[i])!=1) return 1; } fclose(f);
  f=fopen(argv[3],"r"); if(fscanf(f,"%d",&P)!=1) return 1;
  for(int i=0;i<P;i++){ for(int x=0;x<=dd;x++) if(fscanf(f,"%d",&Z[i][x])!=1) return 1;
    int buf[80000],c=0; for(int j=0;j<ND;j++){ long s=0; for(int x=1;x<=dd;x++) s+=(long)Z[i][x]*D[j][x]; if(s) buf[c++]=j; }
    tl[i]=malloc(sizeof(int)*(c+1)); memcpy(tl[i],buf,sizeof(int)*c); nt[i]=c; }
  m=atoi(argv[4]); out=fopen(argv[5],"w"); rec(0,0); printf("count-feasible %d-multisets: %lld\n",m,leaves); return 0; }

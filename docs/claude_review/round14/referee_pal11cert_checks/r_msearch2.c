/* Referee: exhaustive multiset search (m generators from a pool, base row 0) under width counts.
   Candidate lists are filtered by saturation: a candidate touching a direction v with cnt[v] = tau[v] would exceed
   the count, so dropping it is a consequence of the counts only. Every leaf is re-checked in Python. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define W 40
#define MAXD 2560
int ND,dd; int D[MAXD][12]; int tau[MAXD]; int P; int Z[5000][12]; int nt[5000]; int *tl[5000]; u64 bs[5000][W];
int cnt[MAXD]; int m; int cur[16]; long long leaves=0; FILE*out;
void rec(int *cand,int nc,int k){
  if(k==m){ leaves++; for(int i=0;i<m;i++) fprintf(out,"%d ",cur[i]); fprintf(out,"\n"); return; }
  int *nc2=malloc(sizeof(int)*(nc+1));
  for(int a=0;a<nc;a++){ int i=cand[a];
    for(int j=0;j<nt[i];j++) cnt[tl[i][j]]++;
    u64 sat[W]; memset(sat,0,sizeof sat); for(int v=0;v<ND;v++) if(cnt[v]>=tau[v]) sat[v>>6]|=1ULL<<(v&63);
    int n2=0; for(int b=a;b<nc;b++){ int c=cand[b],hit=0; for(int w=0;w<W;w++) if(sat[w]&bs[c][w]){hit=1;break;}
      if(!hit) nc2[n2++]=c; }
    cur[k]=i;
    if(k+1==m || n2>0) rec(nc2,n2,k+1);
    for(int j=0;j<nt[i];j++) cnt[tl[i][j]]--;
  }
  free(nc2);
}
int main(int argc,char**argv){
  dd=atoi(argv[1]); FILE*f=fopen(argv[2],"r"); if(fscanf(f,"%d",&ND)!=1) return 1;
  if(ND>MAXD){fprintf(stderr,"too many dirs\n");return 1;}
  for(int i=0;i<ND;i++){ for(int x=1;x<=dd;x++) if(fscanf(f,"%d",&D[i][x])!=1) return 1; if(fscanf(f,"%d",&tau[i])!=1) return 1; } fclose(f);
  f=fopen(argv[3],"r"); if(fscanf(f,"%d",&P)!=1) return 1;
  for(int i=0;i<P;i++){ for(int x=0;x<=dd;x++) if(fscanf(f,"%d",&Z[i][x])!=1) return 1;
    int buf[MAXD],c=0; for(int j=0;j<ND;j++){ long s=0; for(int x=1;x<=dd;x++) s+=(long)Z[i][x]*D[j][x]; if(s) buf[c++]=j; }
    tl[i]=malloc(sizeof(int)*(c+1)); memcpy(tl[i],buf,sizeof(int)*c); nt[i]=c; memset(bs[i],0,sizeof bs[i]);
    for(int j=0;j<c;j++) bs[i][buf[j]>>6]|=1ULL<<(buf[j]&63); }
  m=atoi(argv[4]); out=fopen(argv[5],"w");
  /* initial candidates: those whose own touches do not exceed (tau>=1 always) */
  int *cand=malloc(sizeof(int)*P); for(int i=0;i<P;i++) cand[i]=i;
  rec(cand,P,0); fclose(out); printf("count-feasible %d-multisets: %lld\n",m,leaves); return 0; }

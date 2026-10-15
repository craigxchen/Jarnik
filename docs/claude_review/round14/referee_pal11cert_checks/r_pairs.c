/* Referee: count-feasible multisets of k pairs e_a-e_b (rows 0..11, 0=base) w.r.t. the 1519 width counts. Pure count pruning. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int ND; int D[2000][12]; int tau[2000];
int E[66][2]; int nt[66]; int tl[66][2000];
int cnt[2000]; long long found[9]; int K;
FILE*out;
int cur[8];
void rec(int start,int k){
  found[k]++;
  if(k==K){ for(int i=0;i<k;i++) fprintf(out,"%d ",cur[i]); fprintf(out,"\n"); return; }
  for(int e=start;e<66;e++){
    int ok=1,j;
    for(j=0;j<nt[e];j++){ int v=tl[e][j]; if(++cnt[v]>tau[v]) ok=0; }
    if(ok){ cur[k]=e; rec(e,k+1); }
    for(j=0;j<nt[e];j++) cnt[tl[e][j]]--;
  }
}
int main(int argc,char**argv){
  K=atoi(argv[1]); out=fopen(argv[2],"w");
  FILE*f=fopen("dirs1519.txt","r"); fscanf(f,"%d",&ND);
  for(int i=0;i<ND;i++){ D[i][0]=0; for(int x=1;x<12;x++) fscanf(f,"%d",&D[i][x]); fscanf(f,"%d",&tau[i]); }
  int n=0; for(int a=0;a<12;a++) for(int b=a+1;b<12;b++){E[n][0]=a;E[n][1]=b;n++;}
  for(int e=0;e<66;e++){ nt[e]=0; for(int i=0;i<ND;i++) if(D[i][E[e][0]]!=D[i][E[e][1]]) tl[e][nt[e]++]=i; }
  memset(cnt,0,sizeof cnt); rec(0,0);
  for(int k=0;k<=K;k++) printf("k=%d feasible multisets %lld\n",k,found[k]);
  return 0;
}

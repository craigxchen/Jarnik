/* Referee: for every count-feasible multiset of k pairs (file), extend by nonpair slots drawn from given pools
   (multiset over each pool; slots of the same pool nondecreasing), keeping all 1519 width counts. Pure count pruning.
   Output each feasible full multiset: pair edge indices | pool:index ... */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXD 1600
int ND; int D[MAXD][12]; int tau[MAXD];
typedef unsigned long long u64;
#define W 25
typedef struct {int n; int (*z)[12]; int *nt; int **tl; u64 (*bs)[W];} Pool;
Pool pools[8]; int npools; int slotpool[8]; int nslots;
int cnt[MAXD]; FILE*out; long long nfeas=0, nnodes=0;
int pe[8]; int K; int sel[8];
void maketouch(int z[12], int *nt, int **tl){
  int buf[MAXD], c=0;
  for(int i=0;i<ND;i++){ long s=0; for(int x=1;x<12;x++) s+= (long)z[x]*D[i][x]; if(s) buf[c++]=i; }
  *tl=malloc(sizeof(int)*(c+1)); memcpy(*tl,buf,sizeof(int)*c); *nt=c;
}
void loadpool(Pool*p,const char*fn){
  FILE*f=fopen(fn,"r"); if(fscanf(f,"%d",&p->n)!=1) exit(1);
  p->z=malloc(sizeof(int[12])*p->n); p->nt=malloc(sizeof(int)*p->n); p->tl=malloc(sizeof(int*)*p->n);
  p->bs=calloc(p->n,sizeof(u64[W]));
  for(int i=0;i<p->n;i++){ for(int x=0;x<12;x++) if(fscanf(f,"%d",&p->z[i][x])!=1) exit(1); maketouch(p->z[i],&p->nt[i],&p->tl[i]);
    for(int j=0;j<p->nt[i];j++) p->bs[i][p->tl[i][j]>>6]|=1ULL<<(p->tl[i][j]&63); }
  fclose(f);
}
int add(Pool*p,int i){ int ok=1; for(int j=0;j<p->nt[i];j++){int v=p->tl[i][j]; if(++cnt[v]>tau[v]) ok=0;} return ok; }
void rem(Pool*p,int i){ for(int j=0;j<p->nt[i];j++) cnt[p->tl[i][j]]--; }
void rec(int s,int start){
  nnodes++;
  if(s==nslots){ nfeas++; for(int i=0;i<K;i++) fprintf(out,"%d ",pe[i]); fprintf(out,"|"); for(int i=0;i<nslots;i++) fprintf(out," %d:%d",slotpool[i],sel[i]); fprintf(out,"\n"); return; }
  Pool*p=&pools[slotpool[s]];
  int st = (s>0 && slotpool[s-1]==slotpool[s]) ? start : 0;
  u64 sat[W]; memset(sat,0,sizeof sat); for(int v=0;v<ND;v++) if(cnt[v]>=tau[v]) sat[v>>6]|=1ULL<<(v&63);
  for(int i=st;i<p->n;i++){ int hit=0; for(int w=0;w<W;w++) if(sat[w]&p->bs[i][w]){hit=1;break;} if(hit) continue;
    if(add(p,i)){ sel[s]=i; rec(s+1,i);} rem(p,i); }
}
int main(int argc,char**argv){
  /* argv: dirs pairsfile k outfile npools pool1 ... nslots slotpool... */
  int a=1; FILE*f=fopen(argv[a++],"r"); if(fscanf(f,"%d",&ND)!=1) return 1;
  for(int i=0;i<ND;i++){ D[i][0]=0; for(int x=1;x<12;x++) if(fscanf(f,"%d",&D[i][x])!=1) return 1; if(fscanf(f,"%d",&tau[i])!=1) return 1;} fclose(f);
  char*pf=argv[a++]; K=atoi(argv[a++]); out=fopen(argv[a++],"w");
  npools=atoi(argv[a++]); for(int i=0;i<npools;i++) loadpool(&pools[i],argv[a++]);
  nslots=atoi(argv[a++]); for(int i=0;i<nslots;i++) slotpool[i]=atoi(argv[a++]);
  /* pairs pool */
  Pool pp; pp.n=66; pp.z=malloc(sizeof(int[12])*66); pp.nt=malloc(sizeof(int)*66); pp.tl=malloc(sizeof(int*)*66);
  int n=0; for(int x=0;x<12;x++) for(int y=x+1;y<12;y++){ memset(pp.z[n],0,sizeof(int[12])); pp.z[n][x]=1; pp.z[n][y]=-1; maketouch(pp.z[n],&pp.nt[n],&pp.tl[n]); n++; }
  FILE*g=fopen(pf,"r"); long long npm=0;
  while(1){
    int ok=1; for(int i=0;i<K;i++) if(fscanf(g,"%d",&pe[i])!=1){ok=0;break;}
    if(!ok) break; npm++;
    memset(cnt,0,sizeof cnt); int good=1;
    for(int i=0;i<K;i++) if(!add(&pp,pe[i])) good=0;
    if(!good){ fprintf(stderr,"infeasible pair multiset in input!\n"); return 2; }
    rec(0,0);
  }
  printf("pair multisets read %lld, nodes %lld, feasible full multisets %lld\n",npm,nnodes,nfeas);
  return 0;
}

// smoothscan: exhaustive scan of all n <= X whose prime factors are all = 1 mod 4 and <= B.
// For each n: the d(n) lattice-point classes mod units are angles <2a-e,theta> mod pi/2, stored as u64
// (2^64 = pi/2). For k = 3..KMAX: C_k(n) = n^(1/4) * least span of k cyclically consecutive angles
// (same measurement as growth/families.md section 6).
// Output: per k and per largest-prime index j (P_+(n) = prime j): best C_k and its n;
//         per k, C0 in list: least n with C_k(n) <= C0, per P_+ index.
// Usage: smoothscan primes.txt B X [shard nshard]   (shard splits on the smallest prime of n's DFS root)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define KMAX 24
#define NPMAX 64
#define DMAX (1<<21)
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static int np; static u64 P[NPMAX]; static u64 TH[NPMAX];
static u64 X; static int shard=0,nshard=1;
static long double bestC[KMAX+1][NPMAX]; static u64 bestN[KMAX+1][NPMAX];
#define NC0 7
static const double C0s[NC0]={0.5,1.0,1.41421356237,2.0,4.0,6.0,8.0};
static u64 firstN[NC0][KMAX+1][NPMAX];
static u64 *buf[64]; static u64 *tmp;
static u64 count=0, sumd=0;
static void rsort(u64 *a, u64 n){
  static u64 cnt[2048];
  u64 *src=a,*dst=tmp;
  for(int pass=0;pass<6;pass++){
    int sh=pass*11; memset(cnt,0,sizeof cnt);
    for(u64 i=0;i<n;i++) cnt[(src[i]>>sh)&2047]++;
    u64 s=0; for(int b=0;b<2048;b++){u64 c=cnt[b];cnt[b]=s;s+=c;}
    for(u64 i=0;i<n;i++) dst[cnt[(src[i]>>sh)&2047]++]=src[i];
    u64 *t=src;src=dst;dst=t;
  }
  // 6 passes: even number, result back in a
}
static const long double HALFPI=1.570796326794896619231321691639751442L;
static void measure(u64 n, u64 *ang, u64 d, int jmax){
  count++; sumd+=d;
  if(d<3) return;
  long double q=powl((long double)n,0.25L)*HALFPI/18446744073709551616.0L;
  int kk = d<KMAX? (int)d : KMAX;
  for(int k=3;k<=kk;k++){
    u64 best=~0ULL;
    for(u64 i=0;i<d;i++){ u64 j=i+k-1; if(j>=d) j-=d; u64 sp=ang[j]-ang[i]; if(sp<best) best=sp; }
    if(k==(int)d) { /* all points: span is full circle minus largest gap; fine as computed */ }
    long double C=q*(long double)best;
    if(C<bestC[k][jmax]){bestC[k][jmax]=C;bestN[k][jmax]=n;}
    for(int c=0;c<NC0;c++) if(C<=C0s[c]){ if(firstN[c][k][jmax]==0||n<firstN[c][k][jmax]) firstN[c][k][jmax]=n; }
  }
}
// dfs: current n with angle list buf[lev] of size d; next primes from index i
static void dfs(int lev,int i,u64 n,u64 d,int jmax){
  if(lev>0) measure(n,buf[lev],d,jmax);
  for(int j=i;j<np;j++){
    if(lev==0 && (j%nshard)!=shard) continue;
    u64 p=P[j]; if((u128)n*p>X) break;
    u64 m=n; int e=0;
    while((u128)m*p<=X){
      m*=p; e++;
      u64 nd=d*(e+1); if(nd>DMAX){fprintf(stderr,"DMAX\n");exit(1);}
      u64 *src=buf[lev],*dst=buf[lev+1]; u64 t=0;
      for(int a=0;a<=e;a++){ u64 sh=(u64)(2*a-e)*TH[j]; for(u64 x=0;x<d;x++) dst[t++]=src[x]+sh; }
      rsort(dst,nd);
      dfs(lev+1,j+1,m,nd,j);
    }
  }
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); u64 B=strtoull(argv[2],0,10); X=(u64)strtold(argv[3],0);
  if(argc>5){shard=atoi(argv[4]);nshard=atoi(argv[5]);}
  int tot; if(fscanf(f,"%d",&tot)!=1) return 1; np=0;
  for(int t=0;t<tot;t++){u64 p,a,b,u; if(fscanf(f,"%llu %llu %llu %llu",&p,&a,&b,&u)!=4) return 1; if(p<=B){P[np]=p;TH[np]=u;np++;}}
  for(int l=0;l<64;l++) buf[l]=malloc(sizeof(u64)*DMAX); tmp=malloc(sizeof(u64)*DMAX);
  for(int k=0;k<=KMAX;k++) for(int j=0;j<NPMAX;j++){bestC[k][j]=1e30L;bestN[k][j]=0;}
  buf[0][0]=0;
  dfs(0,0,1,1,-1);
  printf("# smoothscan B=%llu X=%llu shard=%d/%d count=%llu sumd=%llu\n",B,X,shard,nshard,count,sumd);
  printf("# best: k jmax P_+ bestC n\n");
  for(int k=3;k<=KMAX;k++) for(int j=0;j<np;j++) if(bestN[k][j]) printf("best %d %d %llu %.9Lf %llu\n",k,j,P[j],bestC[k][j],bestN[k][j]);
  printf("# first: C0 k jmax P_+ n\n");
  for(int c=0;c<NC0;c++) for(int k=3;k<=KMAX;k++) for(int j=0;j<np;j++) if(firstN[c][k][j]) printf("first %.4f %d %d %llu %llu\n",C0s[c],k,j,P[j],firstN[c][k][j]);
  return 0;
}

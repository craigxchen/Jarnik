// Independent referee scan: all n <= X with all prime factors = 1 mod 4 and <= B.
// For each n: point classes mod i have angles sum_j (2a_j - e_j) theta_j mod pi/2, a_j in [0,e_j].
// C_k(n) = n^(1/4) * (least span of k cyclically consecutive angles mod pi/2).
// Records best C_k per (decade, k) overall and restricted to P_+(n) <= k^2 and P_+(n) <= k^1.75.
// Usage: rscan B X KMAX countonly
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
typedef long double ld;
static int np; static long long pr[200]; static ld th[200];
static long double X; static int KMAX; static int countonly;
static ld *buf[80]; static long long cap[80];
static long long nnorms=0, nclasses=0;
static double best[20][30][3]; static long long bestn[20][30][3];
static const ld HALFPI = 1.5707963267948966192313216916397514L;
static int cmpld(const void*a,const void*b){ ld x=*(const ld*)a,y=*(const ld*)b; return (x<y)?-1:(x>y);}
static ld *tmp; static long long tmpcap=0;
static void process(long double n, int depth, long long m, long long Pplus){
  nnorms++; nclasses+=m;
  if(countonly) return;
  if(m<2) return;
  if(tmpcap<m){ tmpcap=2*m; tmp=realloc(tmp,tmpcap*sizeof(ld)); }
  memcpy(tmp,buf[depth],m*sizeof(ld));
  qsort(tmp,m,sizeof(ld),cmpld);
  int dec=(int)floorl(log10l(n)+1e-12L);
  ld q=powl(n,0.25L);
  for(int k=2;k<=KMAX && k<=m;k++){
    ld mn=1e9L;
    for(long long i=0;i<m;i++){ long long j=i+k-1; ld d = (j<m)? tmp[j]-tmp[i] : tmp[j-m]+HALFPI-tmp[i]; if(d<mn) mn=d; }
    double C=(double)(q*mn);
    for(int r=0;r<3;r++){
      int ok = (r==0) || (r==1 && Pplus <= (long long)k*k) || (r==2 && (double)Pplus <= pow((double)k,1.75));
      if(ok && (best[dec][k][r]==0 || C<best[dec][k][r])){ best[dec][k][r]=C; bestn[dec][k][r]=(long long)n; }
    }
  }
}
static void dfs(int idx, long double n, int depth, long long m, long long Pplus){
  process(n,depth,m,Pplus);
  for(int i=idx;i<np;i++){
    long double p=pr[i];
    if(n*p > X) break;
    long double nn=n; int e=0;
    while(nn*p <= X){
      nn*=p; e++;
      long long m2=m*(e+1);
      if(!countonly){
        if(cap[depth+1]<m2){ cap[depth+1]=2*m2; buf[depth+1]=realloc(buf[depth+1],cap[depth+1]*sizeof(ld)); }
        long long t=0;
        for(int a=0;a<=e;a++){ ld sh=(2*a-e)*th[i];
          for(long long u=0;u<m;u++){ ld v=buf[depth][u]+sh; v=fmodl(v,HALFPI); if(v<0) v+=HALFPI; buf[depth+1][t++]=v; } }
      }
      dfs(i+1,nn,depth+1,m2,pr[i]);
    }
  }
}
int main(int argc,char**argv){
  int B=atoi(argv[1]); X=strtold(argv[2],0); KMAX=atoi(argv[3]); countonly=atoi(argv[4]);
  for(long long p=5;p<=B;p++){ if(p%4!=1) continue; int isp=1; for(long long d=2;d*d<=p;d++) if(p%d==0){isp=0;break;} if(!isp) continue;
    long long a=0,b=0; for(a=1;a*a<p;a++){ long long r=p-a*a; long long s=(long long)llroundl(sqrtl((ld)r)); if(s*s==r){b=s;break;} }
    pr[np]=p; th[np]=atan2l((ld)b,(ld)a); np++; }
  for(int i=0;i<80;i++){cap[i]=0;buf[i]=0;}
  cap[0]=1; buf[0]=malloc(sizeof(ld)); buf[0][0]=0;
  dfs(0,1.0L,0,1,1);
  printf("# B=%d X=%Lg primes=%d norms(incl n=1)=%lld classes=%lld\n",B,X,np,nnorms,nclasses);
  if(!countonly){
    for(int dec=0;dec<20;dec++) for(int k=2;k<=KMAX;k++){ if(best[dec][k][0]==0) continue;
      printf("dec %d k %d all %.4f n %lld | Pk2 %.4f n %lld | Pk175 %.4f n %lld\n",dec,k,best[dec][k][0],bestn[dec][k][0],best[dec][k][1],bestn[dec][k][1],best[dec][k][2],bestn[dec][k][2]); }
  }
  return 0;
}

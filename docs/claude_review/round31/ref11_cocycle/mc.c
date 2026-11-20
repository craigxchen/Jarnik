// Unbiased Monte Carlo over random pairs (i,j): count k such that {i,j,k} is an
// integral-cofactor configuration at constant C, i.e. diam <= C q^{1/4} gcd(T_ij,T_ik,T_jk),
// |T_ab| = gcd(dx,dy)/2.  Hits are printed and re-verified exactly in Python.
// usage: mc marks.txt C npairs seed
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
typedef __int128 i128; typedef unsigned __int128 u128;
static i128 parse(const char *s){int neg=0; if(*s=='-'){neg=1;s++;} i128 v=0; while(*s>='0'&&*s<='9'){v=v*10+(*s-'0');s++;} return neg?-v:v;}
static u128 ab(i128 x){return x<0?(u128)(-x):(u128)x;}
static int ctz128(u128 x){unsigned long long lo=(unsigned long long)x; if(lo) return __builtin_ctzll(lo); return 64+__builtin_ctzll((unsigned long long)(x>>64));}
static u128 gcd(u128 a,u128 b){ if(!a) return b; if(!b) return a; int sh=ctz128(a|b); a>>=ctz128(a);
  do{ b>>=ctz128(b); if(a>b){u128 t=a;a=b;b=t;} b-=a; }while(b); return a<<sh; }
static unsigned long long rng=88172645463325252ULL; static unsigned long long xs(){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;return rng;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); double C=atof(argv[2]); long npairs=atol(argv[3]); rng^=atol(argv[4])*2654435761ULL;
  char buf[256],b2[128]; fgets(buf,256,f); long double q=strtold(buf,0); long double q4=powl(q,0.25L);
  int N=0; i128 *X=malloc(sizeof(i128)*70000),*Y=malloc(sizeof(i128)*70000);
  while(fscanf(f,"%s %s",buf,b2)==2){X[N]=parse(buf);Y[N]=parse(b2);N++;}
  long hits=0;
  for(long t=0;t<npairs;t++){
    int i=xs()%N, j=xs()%N; if(argc>6){i=atoi(argv[5]);j=atoi(argv[6]);} if(i==j){t--;continue;}
    u128 Tij=gcd(ab(X[j]-X[i]),ab(Y[j]-Y[i]))/2;
    long double dij=hypotl((long double)(X[j]-X[i]),(long double)(Y[j]-Y[i]));
    if((long double)Tij*C*q4 < dij) continue;        // even G = T_ij too small
    for(int k=0;k<N;k++){ if(k==i||k==j) continue;
      u128 g=gcd(Tij, ab(X[k]-X[i]));
      long double dik=hypotl((long double)(X[k]-X[i]),(long double)(Y[k]-Y[i]));
      long double djk=hypotl((long double)(X[k]-X[j]),(long double)(Y[k]-Y[j]));
      long double D=dij; if(dik>D)D=dik; if(djk>D)D=djk;
      if((long double)g*C*q4*1.0001L < D) continue;
      u128 Tik=gcd(ab(X[k]-X[i]),ab(Y[k]-Y[i]))/2, Tjk=gcd(ab(X[k]-X[j]),ab(Y[k]-Y[j]))/2;
      u128 G=gcd(gcd(Tij,Tik),Tjk);
      if((long double)G*C*q4*1.0001L < D) continue;
      hits++; printf("HIT %d %d %d %.6Lg\n",i,j,k,D/((long double)G*q4)); fflush(stdout);
    }
  }
  printf("DONE pairs=%ld hits=%ld N=%d\n",npairs,hits,N);
}

/* Exact minimum of sigma over 4-term characters f = (+1 at 0, +-1 at 1, two more points), F = 0, for sdc_q.
   By AGL(1,q) (sharply 2-transitive, preserves patterns) and negation, every F=0 ternary f with |supp|=4 is
   equivalent to one with f(0)=+1, f(1)!=0. sigma = sum_a |c(a)| - q, c = f*chi. Exact integer arithmetic. */
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char**argv){
  long q=atol(argv[1]); int *ch=malloc(sizeof(int)*q); long m=(q-1)/2;
  for(long a=0;a<q;a++){ if(a==0){ch[a]=0;continue;} long r=1,b=a,e=m; while(e){ if(e&1) r=r*b%q; b=b*b%q; e>>=1;} ch[a]=(r==1)?1:-1; }
  int *s=malloc(sizeof(int)*2*q); /* s01(a) for f1=-1 and f1=+1, duplicated index helpers */
  int *chx=malloc(sizeof(int)*2*q); for(long a=0;a<2*q;a++) chx[a]=ch[a%q];
  long best=1L<<60; long bx=-1,by=-1; int bcfg=-1; long cnt=0;
  for(int cfg=0;cfg<3;cfg++){
    int f1 = (cfg==0)? 1 : -1; int fx,fy;
    if(cfg==0){fx=-1;fy=-1;} else if(cfg==1){fx=1;fy=-1;} else {fx=-1;fy=1;}
    for(long a=0;a<q;a++) s[a]=ch[a]+f1*ch[(a-1+q)%q];
    for(long x=2;x<q;x++){
      long ystart = (cfg==0)? x+1 : 2;
      for(long y=ystart;y<q;y++){ if(y==x) continue;
        /* c(a)=s[a]+fx*ch[a-x]+fy*ch[a-y]; ch[a-x] = chx[a-x+q] */
        long N=0; const int *px=chx+q-x, *py=chx+q-y;
        for(long a=0;a<q;a++){ int v=s[a]+fx*px[a]+fy*py[a]; N+= v<0?-v:v; }
        cnt++;
        if(N<best){best=N;bx=x;by=y;bcfg=cfg;}
      }
    }
  }
  printf("q=%ld configs=%ld min sigma(n=4,F=0)=%ld sigma/r=%.6f witness: f(0)=1 cfg=%d x=%ld y=%ld\n",q,cnt,best-q,(double)(best-q)/(2.0*q),bcfg,bx,by);
  return 0;
}

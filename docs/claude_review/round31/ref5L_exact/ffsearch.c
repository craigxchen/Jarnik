/* Exhaustive function-field test of the Part V rigidity theorem (and of its peeling extension).
   Model: k=F_p, form I(X,Y)=XY (split; the proof of the theorem is form-agnostic), columns
   V_i=(z_i,z'_i) with z_i,z'_i monic quadratics in F_p[T], blocks J_e = T - a_e (ten distinct a_e),
   det(V_i,V_j)=z_i z'_j - z'_i z_j = t_ij (T - a_ij), t_ij in F_p^* (constant cofactors),
   n_i = prod_{j!=i}(T-a_ij), delta_i = z_i z'_i - n_i (deg<=3 automatically).  "P" = T, so
   P^{s-eps} <-> deg delta <= s-1.  We enumerate ALL such tuples with deg delta_i <= 2 for all i
   (the IV.5 range), up to translation T->T+c (beta_1=0) and scaling T->lambda T (sigma in {0,1}).
   Parallel leading parts force beta'_i=beta_i+sigma, gamma'_i=gamma_i+sigma*beta_i+tau.          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int p;
static int md(long x){x%=p; if(x<0)x+=p; return (int)x;}
static int inv(int x){int r=1,b=x,e=p-2; while(e){if(e&1)r=md((long)r*b); b=md((long)b*b); e>>=1;} return r;}
static int B[6],G[6],Bp[6],Gp[6]; static int sig,tau;
static int ev(int b,int g,int x){return md((long)x*x+(long)b*x+g);}
/* det coefficients: z_i z'_j - z'_i z_j ; returns L (T coeff), K (const) ; asserts higher vanish */
static void detc(int i,int j,int*L,int*K){
  /* z=T^2+bT+g */
  long t3=(Bp[j]+B[i])-(B[j]+Bp[i]);
  long t2=(Gp[j]+(long)B[i]*Bp[j]+G[i])-(G[j]+(long)Bp[i]*B[j]+Gp[i]);
  if(md(t3)||md(t2)){fprintf(stderr,"structure error\n");exit(1);}
  *L=md((long)B[i]*Gp[j]+(long)G[i]*Bp[j]-(long)Bp[i]*G[j]-(long)Gp[i]*B[j]);
  *K=md((long)G[i]*Gp[j]-(long)Gp[i]*G[j]);
}
static long cnt[4][11][4]; /* [defect class][#iso edges][graph class] */
static long viol[3]; static long total=0, nonprim=0, propfail=0;
static int E[10][2]={{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
static int eid(int i,int j){for(int e=0;e<10;e++) if((E[e][0]==i&&E[e][1]==j)||(E[e][0]==j&&E[e][1]==i)) return e; return -1;}
static int degS(int *iso,int v,int r1,int r2){int d=0; for(int e=0;e<10;e++) if(iso[e]&&(E[e][0]==v||E[e][1]==v)){int o=E[e][0]==v?E[e][1]:E[e][0]; if(o!=r1&&o!=r2) d++;} return d;}
static int peel(int *iso,int s){ for(int a=0;a<5;a++)for(int j=0;j<5;j++)for(int k=0;k<5;k++){ if(a==j||a==k||j==k)continue;
   if(degS(iso,a,-1,-1)>=s && degS(iso,j,a,-1)>=s-1 && degS(iso,k,a,j)>=s-2) return 1;} return 0;}
static FILE *ex;
static void survivor(void){
  int a[10],L[10];
  for(int e=0;e<10;e++){int K; detc(E[e][0],E[e][1],&L[e],&K); if(!L[e]) return; a[e]=md((long)-K*inv(L[e]));}
  for(int e=0;e<10;e++)for(int f=e+1;f<10;f++) if(a[e]==a[f]) return;
  /* primitivity */
  for(int i=0;i<5;i++){ if(B[i]==Bp[i]&&G[i]==Gp[i]){nonprim++;return;}
    for(int x=0;x<p;x++) if(ev(B[i],G[i],x)==0&&ev(Bp[i],Gp[i],x)==0){nonprim++;return;} }
  total++;
  int dmax=-1; int delta[5][5];
  for(int i=0;i<5;i++){
    long n[5]={1,0,0,0,0}; /* product of (T-a) */
    for(int e=0;e<10;e++) if(E[e][0]==i||E[e][1]==i){ long m[5]={0}; for(int d=0;d<4;d++){ m[d+1]+=n[d]; m[d]-= (long)a[e]*n[d]; } for(int d=0;d<5;d++) n[d]=md(m[d]); }
    long zz[5]={0}; int zc[3]={G[i],B[i],1}, zpc[3]={Gp[i],Bp[i],1};
    for(int u=0;u<3;u++)for(int v=0;v<3;v++) zz[u+v]+=(long)zc[u]*zpc[v];
    for(int d=0;d<5;d++){delta[i][d]=md(zz[d]-n[d]); if(delta[i][d]&&d>dmax) dmax=d;}
  }
  if(dmax>2){fprintf(stderr,"defect degree error\n");exit(1);}
  int iso[10]; int niso=0;
  for(int e=0;e<10;e++){ int i=E[e][0],j=E[e][1];
    int ii=(ev(B[i],G[i],a[e])==0)||(ev(Bp[i],Gp[i],a[e])==0);
    int jj=(ev(B[j],G[j],a[e])==0)||(ev(Bp[j],Gp[j],a[e])==0);
    if(ii!=jj) propfail++; iso[e]=ii&&jj; niso+=iso[e]; }
  int cls = dmax<0?0: dmax==0?1: dmax==1?2:3;   /* 0 exact, 1 const, 2 deg<=1, 3 deg<=2 */
  int g = peel(iso,3)?3: peel(iso,2)?2: peel(iso,1)?1:0;
  cnt[cls][niso][g]++;
  if(cls==1 && g>=1) viol[0]++;
  if(cls>=1 && cls<=2 && g>=2) viol[1]++;
  if(cls>=1 && g>=3) viol[2]++;
  if(cls>=1 && niso>=2 && ex){ fprintf(ex,"cls=%d niso=%d graph=%d sig=%d tau=%d |",cls,niso,g,sig,tau);
     for(int i=0;i<5;i++) fprintf(ex," (%d,%d;%d,%d)",B[i],G[i],Bp[i],Gp[i]); fprintf(ex," | iso:");
     for(int e=0;e<10;e++) if(iso[e]) fprintf(ex," %d%d",E[e][0]+1,E[e][1]+1); fprintf(ex,"\n"); }
}
int main(int argc,char**argv){
  p=atoi(argv[1]); ex=fopen(argv[2],"w");
  for(sig=0;sig<2;sig++) for(tau=0;tau<p;tau++){
   for(int g1=0;g1<p;g1++){ B[0]=0;G[0]=g1;
    for(long m2=0;m2<(long)p*p;m2++){ B[1]=m2/p;G[1]=m2%p;
     for(long m3=0;m3<(long)p*p;m3++){ B[2]=m3/p;G[2]=m3%p;
      for(long m4=0;m4<(long)p*p;m4++){ B[3]=m4/p;G[3]=m4%p;
        for(int i=0;i<4;i++){Bp[i]=md(B[i]+sig);Gp[i]=md(G[i]+(long)sig*B[i]+tau);}
        int a6[4][4]; int ok=1;
        for(int i=0;i<4&&ok;i++)for(int j=i+1;j<4;j++){int L,K; detc(i,j,&L,&K); if(!L){ok=0;break;} a6[i][j]=a6[j][i]=md((long)-K*inv(L));}
        if(!ok) continue;
        int a5[4];
        for(int i=0;i<4;i++){ long s=-(2L*B[i]+sig); for(int j=0;j<4;j++) if(j!=i) s-=a6[i][j]; a5[i]=md(s);}
        /* mark 5: need z_i(a) z'_5(a) - z'_i(a) z_5(a) = 0 at a=a5[i], i=0..3 ; linear in (b5,g5) */
        for(int b5=0;b5<p;b5++){
          int x=a5[0]; int zi=ev(B[0],G[0],x), zpi=ev(Bp[0],Gp[0],x);
          /* f = (zi-zpi)*z5(x) + zi*(sig*x+sig*b5+tau) ; z5(x)=x^2+b5 x+g5 */
          int cg=md(zi-zpi); int lo,hi;
          int gfix=-1;
          if(cg){ long rest=md((long)cg*md((long)x*x+(long)b5*x)+(long)zi*md((long)sig*x+(long)sig*b5+tau)); gfix=md(-rest*inv(cg)); lo=gfix;hi=gfix;} else {lo=0;hi=p-1;}
          for(int g5=lo;g5<=hi;g5++){
            B[4]=b5;G[4]=g5;Bp[4]=md(b5+sig);Gp[4]=md(g5+(long)sig*b5+tau);
            int good=1;
            for(int i=0;i<4&&good;i++){int L,K; detc(i,4,&L,&K); if(!L||md((long)L*a5[i]+K)) good=0;}
            if(!good) continue;
            long s=-(long)sig; for(int i=0;i<4;i++) s-=a5[i]; if(md(2L*b5-s)) continue; /* T^3 at mark 5 */
            survivor();
          }
        }
      }}}}}
  printf("p=%d tuples (10 distinct blocks, primitive, deg delta<=2): %ld ; nonprimitive skipped %ld ; propagation failures %ld\n",p,total,nonprim,propfail);
  const char*cn[4]={"exact","deg0 ","deg1 ","deg2 "};
  for(int c=0;c<4;c++){ printf("%s:",cn[c]); for(int k=0;k<=10;k++){ long s=0; for(int g=0;g<4;g++) s+=cnt[c][k][g]; if(s) printf("  %dedges:%ld[g0=%ld,g1=%ld,g2=%ld,g3=%ld]",k,s,cnt[c][k][0],cnt[c][k][1],cnt[c][k][2],cnt[c][k][3]);} printf("\n"); }
  printf("violations: T1(deg0 & >=1 edge)=%ld  T2(deg<=1 & nonstar triple)=%ld  T3(deg<=2 & 3-2-1)=%ld\n",viol[0],viol[1],viol[2]);
  return 0;
}

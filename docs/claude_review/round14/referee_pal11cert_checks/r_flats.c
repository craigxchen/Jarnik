/* Referee: brute force over supports R (|R|>=2) and subsets Zc of the 11 core Hadamard columns.
   Matrix M = rows {1|_R} u {beta_c|_R : c in Zc}.  Generic element of null(M) (within R^R) has full support iff no e_x in rowspace.
   t(generic) = #{c : beta_c|_R not in rowspace(M)}. Output closed flats (|Zc| = 11-t) with full-support complement.
   All ranks mod p = 2^31-1; minors of 0/+-1 matrices of size <=13 are < 13^6.5 < p, so ranks = rational ranks. */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
typedef int64_t i64; const i64 P=2147483647LL;
static i64 md(i64 a){a%=P; return a<0?a+P:a;}
static i64 pw(i64 a,i64 e){i64 r=1;a=md(a);while(e){if(e&1)r=r*a%P;a=a*a%P;e>>=1;}return r;}
int B[12][12];
int main(){
  FILE*f=fopen("hadB.txt","r"); for(int i=0;i<12;i++)for(int j=0;j<12;j++)fscanf(f,"%d",&B[i][j]); fclose(f);
  /* B[0]=all ones, B[1..11]=core betas */
  for(int R=0;R<4096;R++){
    int cols[12],s=0; for(int x=0;x<12;x++) if(R>>x&1) cols[s++]=x;
    if(s<2) continue;
    for(int Z=0;Z<2048;Z++){
      i64 m[13][12]; int nr=0;
      for(int j=0;j<s;j++) m[nr][j]=B[0][cols[j]]; nr++;
      for(int c=0;c<11;c++) if(Z>>c&1){ for(int j=0;j<s;j++) m[nr][j]=md(B[c+1][cols[j]]); nr++; }
      /* RREF */
      int piv[12]; int r=0;
      for(int j=0;j<s && r<nr;j++){
        int pr=-1; for(int i=r;i<nr;i++) if(m[i][j]){pr=i;break;}
        if(pr<0) continue;
        if(pr!=r) for(int k=0;k<s;k++){i64 t=m[pr][k];m[pr][k]=m[r][k];m[r][k]=t;}
        i64 inv=pw(m[r][j],P-2); for(int k=0;k<s;k++) m[r][k]=m[r][k]*inv%P;
        for(int i=0;i<nr;i++) if(i!=r && m[i][j]){ i64 fct=m[i][j]; for(int k=0;k<s;k++) m[i][k]=md(m[i][k]-fct*m[r][k]); }
        piv[r]=j; r++;
      }
      int dimnull=s-r; if(dimnull==0) continue;
      /* e_x in rowspace iff column x is pivot column whose row is zero on all nonpivot columns */
      int ispiv[12]={0}; int prow[12]; for(int i=0;i<r;i++){ispiv[piv[i]]=1;prow[piv[i]]=i;}
      int full=1;
      for(int j=0;j<s && full;j++) if(ispiv[j]){ int i=prow[j]; int zero=1; for(int k=0;k<s;k++) if(!ispiv[k] && m[i][k]) {zero=0;break;} if(zero) full=0; }
      if(!full) continue;
      /* t: beta_c in rowspace? reduce */
      int t=0, closedcount=0;
      for(int c=0;c<11;c++){
        i64 w[12]; for(int j=0;j<s;j++) w[j]=md(B[c+1][cols[j]]);
        for(int i=0;i<r;i++){ i64 fct=w[piv[i]]; if(fct) for(int k=0;k<s;k++) w[k]=md(w[k]-fct*m[i][k]); }
        int nz=0; for(int k=0;k<s;k++) if(w[k]) {nz=1;break;}
        if(nz) t++; else closedcount++;
      }
      if(closedcount!=__builtin_popcount(Z)) continue; /* not closed */
      /* also print null vector if dimnull==1: free column fc, x_fc=1, x_piv = -m[i][fc] */
      printf("%d %d %d %d %d",R,Z,s,t,dimnull);
      if(dimnull==1){ int fc=-1; for(int k=0;k<s;k++) if(!ispiv[k]) fc=k;
        i64 v[12]={0}; v[cols[fc]]=1; for(int i=0;i<r;i++) v[cols[piv[i]]]=md(-m[i][fc]);
        for(int x=0;x<12;x++){ i64 a=v[x]; if(a>P/2) a-=P; printf(" %lld",(long long)a);} }
      printf("\n");
    }
  }
  return 0;
}

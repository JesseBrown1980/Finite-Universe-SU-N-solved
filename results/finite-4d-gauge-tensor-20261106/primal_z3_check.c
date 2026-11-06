/* Independent integer-only PRIMAL link-flux verification of the 2x2x2 Z3
 * spatial magnetic partition function. No Fourier dual, Python import, or
 * transfer implementation. Enumerates every DISTINCT plaquette-flux vector:
 * 3^14 rather than the 3^24 raw links or 3^17 gauge orbit states.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { N=24, D=14 };
static unsigned char field[N][N];

static unsigned mod3(int x) { int v=x%3;return (unsigned)(v<0?v+3:v); }
static unsigned inv3(unsigned x) { if(x==1)return 1;if(x==2)return 2;abort(); }
static int edge(int x,int mu) {return x*3+mu;}
static void populate(void) {
  int face=0;
  for(int x=0;x<8;x++)for(int a=0;a<3;a++)for(int b=a+1;b<3;b++){
    int e[4]={edge(x,a),edge(x^(1<<a),b),edge(x^(1<<b),a),edge(x,b)};
    int s[4]={1,1,-1,-1};
    for(int k=0;k<4;k++)field[face][e[k]]=(unsigned char)mod3(field[face][e[k]]+s[k]);
    face++;
  }
  if(face!=N)abort();
}
static int pivots[D];
static int find_rank(void) {
  unsigned char a[N][N];
  for(int i=0;i<N;i++)for(int j=0;j<N;j++)a[i][j]=field[i][j];
  int rank=0;
  for(int col=0;col<N;col++){
    int k=rank;
    while(k<N && a[k][col]==0)k++;
    if(k==N)continue;
    for(int j=0;j<N;j++){unsigned char t=a[rank][j];a[rank][j]=a[k][j];a[k][j]=t;}
    unsigned invert=inv3(a[rank][col]);
    for(int j=0;j<N;j++)a[rank][j]=(unsigned char)mod3((int)a[rank][j]*(int)invert);
    for(int i=0;i<N;i++)if(i!=rank && a[i][col]){
      unsigned mul=a[i][col];
      for(int j=0;j<N;j++)a[i][j]=(unsigned char)mod3((int)a[i][j]-(int)mul*(int)a[rank][j]);
    }
    if(rank>=D)abort();
    pivots[rank]=col;rank++;
  }
  return rank;
}
static void print128(unsigned __int128 v){
  char s[64];int i=0;
  do{s[i++]=(char)('0'+v%10);v/=10;}while(v);
  while(i)putchar(s[--i]);
}
int main(void){
  populate();
  if(find_rank()!=D)abort();
  unsigned char digits[D]={0};
  unsigned char flux[N]={0};
  uint64_t n=1;
  for(int i=0;i<D;i++)n*=3;
  if(n!=4782969ULL)abort();
  uint64_t sum=0;
  int64_t wilson1_twice=0;
  int64_t wilson2_twice=0;
  uint64_t cnt=0;
  for(uint64_t code=0;code<n;code++){
    unsigned zeros=0;
    for(int p=0;p<N;p++)zeros+=(unsigned)(flux[p]==0);
    uint64_t w=((uint64_t)1)<<(2*zeros);
    sum+=w;
    wilson1_twice+=(flux[0]==0 ? (int64_t)(2*w) : -(int64_t)w);
    unsigned combined=(flux[0]+flux[1])%3;
    wilson2_twice+=(combined==0 ? (int64_t)(2*w) : -(int64_t)w);
    cnt++;
    if(code+1==n)break;
    for(int k=0;k<D;k++){
      digits[k]=(unsigned char)((digits[k]+1)%3);
      for(int p=0;p<N;p++)flux[p]=(unsigned char)((flux[p]+field[p][pivots[k]])%3);
      if(digits[k]!=0)break;
    }
  }
  if(cnt!=4782969ULL || sum!=474407207856378ULL || wilson1_twice!=803118340327086LL || wilson2_twice!=708536247566832LL)abort();
  unsigned __int128 raw=(unsigned __int128)sum;
  for(int i=0;i<10;i++)raw*=3;
  printf("PRIMAL_Z3|rank=14|flux_states=%llu|sum_flux=%llu|raw_M=",
        (unsigned long long)cnt,(unsigned long long)sum);
  print128(raw);
  printf("|twice_W1_sum=%lld|twice_W2_sum=%lld|match_frozen=1|json=0\n",
      (long long)wilson1_twice,(long long)wilson2_twice);
  return 0;
}

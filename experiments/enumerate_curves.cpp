// Exact genus-two F_p-isomorphism census, p > 3.
// Normalized squarefree binary sextics, modulo PGL_2 and square scalars.
// Compile: g++ -O3 -std=c++17 enumerate_curves.cpp -o enumerate_curves
// Run: ./enumerate_curves 7 output_prefix
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <string>
using namespace std;
int p,nu; long long pd[8]; vector<int> invs,chi;
int mod(long long n){n%=p;return int(n<0?n+p:n);}
using Poly=array<int,7>;
int deg(const vector<int>&f){int d=int(f.size())-1;while(d>=0&&!f[d])--d;return d;}
bool smooth(const Poly&f){
 vector<int>a(f.begin(),f.end()),b(6);for(int i=1;i<=6;i++)b[i-1]=mod(i*f[i]);
 while(deg(b)>=0){int db=deg(b);while(deg(a)>=db){int da=deg(a),c=mod(a[da]*invs[b[db]]);for(int j=0;j<=db;j++)a[j+da-db]=mod(a[j+da-db]-c*b[j]);}a.swap(b);}
 return deg(a)==0;
}
Poly decode(long long id){
 Poly f{}; long long block=pd[6]+pd[5];int lead=id>=block?nu:1;id%=block;
 int d=6;if(id>=pd[6]){d=5;id-=pd[6];}for(int i=0;i<d;i++){f[i]=id%p;id/=p;}f[d]=lead;return f;
}
long long encode(Poly f){
 int d=f[6]?6:5;assert(f[d]);int lead=chi[f[d]]==1?1:nu;int c=mod(lead*invs[f[d]]);
 long long id=0;for(int i=0;i<d;i++)id+=mod(f[i]*c)*pd[i];
 if(d==5)id+=pd[6];if(lead==nu)id+=pd[6]+pd[5];return id;
}
using Action=array<int,49>;
vector<Action> actions(){
 vector<Action> out;
 for(int a=0;a<p;a++)for(int b=0;b<p;b++)for(int c=0;c<p;c++)for(int d=0;d<p;d++){
  if(!mod(a*d-b*c))continue;
  int first=a?a:(b?b:(c?c:d));if(first!=1)continue;
  int powers[2][7][7]={};powers[0][0][0]=powers[1][0][0]=1;
  for(int i=1;i<=6;i++)for(int j=0;j<=i;j++){
   powers[0][i][j]=mod((j? a*powers[0][i-1][j-1]:0)+b*powers[0][i-1][j]);
   powers[1][i][j]=mod((j? c*powers[1][i-1][j-1]:0)+d*powers[1][i-1][j]);
  }
  Action t{};for(int i=0;i<=6;i++)for(int j=0;j<=i;j++)for(int k=0;k<=6-i;k++)t[(j+k)*7+i]=mod(t[(j+k)*7+i]+powers[0][i][j]*powers[1][6-i][k]);out.push_back(t);
 }return out;
}
Poly transform_form(const Action&t,const Poly&f){Poly g{};for(int j=0;j<=6;j++){int s=0;for(int i=0;i<=6;i++)s+=t[j*7+i]*f[i];g[j]=mod(s);}return g;}
array<int,4> cm(const Poly&f){
 vector<int> h(1,1);for(int k=0;k<(p-1)/2;k++){vector<int> n(h.size()+6);for(size_t i=0;i<h.size();i++)for(int j=0;j<=6;j++)n[i+j]=mod(n[i+j]+h[i]*f[j]);h.swap(n);}
 return {h[p-1],h[p-2],h[2*p-1],h[2*p-2]};
}
int rank2(array<int,4>m){if(mod(m[0]*m[3]-m[1]*m[2]))return 2;return any_of(m.begin(),m.end(),[](int a){return a!=0;})?1:0;}
pair<int,int> pointcounts(const Poly&f){
 int d=f[6]?6:5;int n1=p+(d==5?1:1+chi[f[6]]),n2=p*p+(d==5?1:2);
 for(int x=0;x<p;x++){int y=0;for(int i=d;i>=0;i--)y=mod(y*x+f[i]);n1+=chi[y];}
 for(int a=0;a<p;a++)for(int b=0;b<p;b++){
  int u=0,v=0;for(int i=d;i>=0;i--){int un=mod(u*a+nu*v*b+f[i]);v=mod(u*b+v*a);u=un;}n2+=chi[mod(u*u-nu*v*v)];
 }return {n1,n2};
}
int main(int argc,char**argv){
 if(argc!=3){cerr<<"usage: enumerate_curves p output_prefix\n";return 2;}p=stoi(argv[1]);assert(p>3&&p<100);
 pd[0]=1;for(int i=1;i<8;i++)pd[i]=pd[i-1]*p;invs.resize(p);chi.assign(p,-1);chi[0]=0;
 for(int i=1;i<p;i++){chi[i*i%p]=1;for(int j=1;j<p;j++)if(i*j%p==1)invs[i]=j;}
 nu=2;while(chi[nu]!=-1)nu++;
 auto ts=actions();assert(ts.size()==size_t(p*(p*p-1)));
 long long size=2*(pd[6]+pd[5]);vector<int> orbit(size,-1);
 string prefix=argv[2];ofstream out(prefix+"_curves.csv");out<<"class_id,c0,c1,c2,c3,c4,c5,c6,p_rank,a_number,aut_order,orbit_size,n1,n2,weil_a1,weil_a2\n";
 int classes=0,counts[3][3]={};long long n_smooth=0;vector<long long> reps;
 for(long long id=0;id<size;id++){
  if(orbit[id]!=-1)continue;Poly f=decode(id);if(!smooth(f)){orbit[id]=-2;continue;}
  int stab=0,norb=0;
  for(auto&t:ts){long long j=encode(transform_form(t,f));assert(j>=0&&j<size);if(j==id)stab++;if(orbit[j]==-1){orbit[j]=classes;norb++;}else assert(orbit[j]==classes);}
  assert(stab>0&&norb*stab==int(ts.size()));n_smooth+=norb;
  auto m=cm(f);array<int,4>m2={mod(m[0]*m[0]+m[1]*m[2]),mod(m[0]*m[1]+m[1]*m[3]),mod(m[2]*m[0]+m[3]*m[2]),mod(m[2]*m[1]+m[3]*m[3])};
  int prank=rank2(m2),anum=2-rank2(m);counts[prank][anum]++;
  out<<classes;for(int c:f)out<<','<<c;out<<','<<prank<<','<<anum<<','<<2*stab<<','<<norb;
  if(prank==0){auto[n1,n2]=pointcounts(f);int s1=p+1-n1,s2=p*p+1-n2;out<<','<<n1<<','<<n2<<','<<-s1<<','<<(s1*s1-s2)/2;}else out<<",,,,";out<<'\n';
  reps.push_back(id);classes++;
 }
 out.close();ofstream mapout(prefix+"_orbit.bin",ios::binary);mapout.write((char*)orbit.data(),orbit.size()*sizeof(int));
 ofstream summary(prefix+"_summary.txt");summary<<"p="<<p<<"\nnon_square="<<nu<<"\nnormalized_models="<<size<<"\nsmooth_models="<<n_smooth<<"\nisomorphism_classes="<<classes<<"\n";
 for(int r=0;r<=2;r++)for(int a=0;a<=2;a++)if(counts[r][a])summary<<"p_rank="<<r<<",a_number="<<a<<",count="<<counts[r][a]<<"\n";
 cout<<"p="<<p<<" classes="<<classes<<" smooth_models="<<n_smooth<<" supersingular="<<counts[0][1]+counts[0][2]<<" superspecial="<<counts[0][2]<<endl;
 return 0;
}

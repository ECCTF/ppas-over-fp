// Exhaustive maximal genus-two models, independent of isogeny connectivity.
// Maximality over Fp2 makes all six branch points Fp2-rational. Cover the
// four Frobenius cycle types using PGL2(Fp) normalizations; retain both twists.
// No graph input, class-group input, or connectedness theorem is used.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int p,nu; vector<int> chi,iv; long long tried=0,traceok=0,cmok=0,maximal=0;
ofstream out;
int mod(long long x){int v=x%p;return v<0?v+p:v;}
using Poly=vector<int>;
Poly mul(const Poly&a,const Poly&b,int cap=10000){
 int n=min(int(a.size()+b.size()-1),cap);vector<long long> raw(n,0);
 for(int i=0;i<int(a.size());++i)for(int j=0;j<int(b.size())&&i+j<n;++j)raw[i+j]+=a[i]*b[j];
 Poly c(n);for(int i=0;i<n;++i)c[i]=raw[i]%p;return c;
}
struct F{int a,b;};
F add(F x,F y){return {mod(x.a+y.a),mod(x.b+y.b)};}
F mult(F x,F y){return {mod(x.a*y.a+nu*x.b*y.b),mod(x.a*y.b+x.b*y.a)};}
F divide(F x,F y){int z=iv[mod(y.a*y.a-nu*y.b*y.b)];return mult(x,{mod(y.a*z),mod(-y.b*z)});}
struct Quad{int a,b;Poly f;};vector<Quad> qs;
int qindex(F x){if(x.b>p/2)x.b=p-x.b;assert(x.b);return x.a*((p-1)/2)+x.b-1;}
void check(const Poly&f,int type){
 ++tried;int d=f.size()-1,sum=0;
 for(int x=0;x<p;++x){int v=0;for(int i=d;i>=0;--i)v=(v*x+f[i])%p;sum+=chi[v];}
 if(sum!=(d==5?0:-1))return;
 ++traceok;Poly h={1},b=f;int n=(p-1)/2;
 while(n){if(n&1)h=mul(h,b,2*p);n>>=1;if(n)b=mul(b,b,2*p);}
 if(h[p-1]||h[p-2]||h[2*p-1]||h[2*p-2])return;
 ++cmok;int n2=p*p+(d==5?1:2);
 for(int a=0;a<p;++a)for(int b=0;b<=(p-1)/2;++b){
  F v={0,0},x={a,b};for(int i=d;i>=0;--i)v=add(mult(v,x),{f[i],0});
  n2+=(b?2:1)*chi[mod(v.a*v.a-nu*v.b*v.b)];
 }
 if(n2!=p*p+1+4*p)return;
 for(int scalar:{1,nu}){++maximal;out<<type;for(int i=0;i<7;++i)out<<','<<(i<=d?f[i]*scalar%p:0);out<<'\n';}
}
void status(int type){cerr<<"p="<<p<<" type="<<type<<" models="<<tried<<" trace_zero="<<traceok<<" cm_zero="<<cmok<<" maximal_models_with_twists="<<maximal<<endl;}
int main(int argc,char**argv){
 assert(argc==3);p=stoi(argv[1]);assert(p>=7&&p<1000);
 chi.assign(p,-1);chi[0]=0;iv.resize(p);
 for(int x=1;x<p;++x){chi[x*x%p]=1;for(int y=1;y<p;++y)if(x*y%p==1)iv[x]=y;}
 nu=2;while(chi[nu]!=-1)++nu;
 for(int a=0;a<p;++a)for(int b=1;b<=(p-1)/2;++b)qs.push_back({a,b,{mod(a*a-nu*b*b),mod(-2*a),1}});
 out.open(argv[2]);out<<"rational_branch_points,c0,c1,c2,c3,c4,c5,c6\n";
 auto start=chrono::steady_clock::now();
 // Six rational points: normalize three to infinity,0,1.
 for(int a=2;a<p;++a)for(int b=a+1;b<p;++b)for(int c=b+1;c<p;++c)
  check(mul(mul(mul({0,p-1,1},{p-a,1}),{p-b,1}),{p-c,1}),6);
 status(6);
 // Four rational points: normalize three to infinity,0,1.
 for(int a=2;a<p;++a){Poly base=mul({0,p-1,1},{p-a,1});for(auto&q:qs)check(mul(base,q.f),4);}
 status(4);
 // Two rational points: normalize to infinity,0. Scale a chosen nonrational
 // root a+b*u to real part 1 (a!=0), or to u (a==0).
 vector<int> first={0};for(int b=1;b<=(p-1)/2;++b)first.push_back(qindex({1,b}));
 for(int i:first)for(int j=0;j<int(qs.size());++j)if(i!=j)check(mul({0,1},mul(qs[i].f,qs[j].f)),2);
 status(2);
 // No rational points: PGL2(Fp) is transitive on conjugate pairs. Fix {u,-u}.
 // Its stabilizer is x->(a*x+nu)/(x+a), x->x, together with x->-x.
 // Choose a representative for the second pair under that stabilizer.
 vector<bool> seen(qs.size(),false);seen[0]=true;vector<int> reps;
 for(int i=1;i<int(qs.size());++i)if(!seen[i]){
  reps.push_back(i);F x={qs[i].a,qs[i].b};
  for(int a=-1;a<p;++a){F y=a==-1?x:divide(add(mult({a,0},x),{nu,0}),add(x,{a,0}));
   seen[qindex(y)]=true;seen[qindex({mod(-y.a),mod(-y.b)})]=true;}
 }
 for(int i:reps)for(int j=1;j<int(qs.size());++j)if(i!=j)check(mul(mul(qs[0].f,qs[i].f),qs[j].f),0);
 status(0);out.close();
 cout<<"p="<<p<<" nonsquare="<<nu<<" second_pair_orbits="<<reps.size()<<" candidates="<<tried<<" maximal_models="<<maximal
     <<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
}

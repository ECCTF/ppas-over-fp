"""Exact 10-bit experiments without optional CAS dependencies.

Compute the complete closed Weil-restriction (3,3) sector for p=1 mod 3.
Every edge is an individual geometric elliptic 3-kernel, NOT a Galois
orbit of kernels.  Each vertex denotes the restriction of the unique
maximal elliptic form attached to its Frobenius orbit of j-invariants.
Auxiliary Weierstrass models are used only for the geometric quotient j.
For generic j their possible quadratic twist does not change these data.

Usage: python ten_bit_experiments.py [p ...]
Only Python's standard library is required.  Plotting is separate.
"""
import csv, itertools, json, math, random, sys, time
from collections import Counter, deque
from pathlib import Path

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'

class Field:
    def __init__(self,p):
        self.p=p; self.q=p*p
        self.nu=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
    def add(self,x,y):
        p=self.p
        return (x%p+y%p)%p+p*((x//p+y//p)%p)
    def neg(self,x):
        p=self.p
        return (-x%p)%p+p*((-(x//p))%p)
    def sub(self,x,y):return self.add(x,self.neg(y))
    def mul(self,x,y):
        p=self.p; a,b=x%p,x//p; c,d=y%p,y//p
        return (a*c+self.nu*b*d)%p+p*((a*d+b*c)%p)
    def power(self,x,n):
        r=1
        while n:
            if n&1:r=self.mul(r,x)
            x=self.mul(x,x);n>>=1
        return r
    def inv(self,x):
        p=self.p;a,b=x%p,x//p
        t=pow((a*a-self.nu*b*b)%p,-1,p)
        return (a*t)%p+p*((-b*t)%p)
    def div(self,x,y):return self.mul(x,self.inv(y))
    def frob(self,x):return x%self.p+self.p*((-(x//self.p))%self.p)
    def coeff(self,x):return [x%self.p,x//self.p]

def trim(f):
    while len(f)>1 and not f[-1]:f.pop()
    return f

class Polys:
    def __init__(self,F):self.F=F;self.rng=random.Random(F.p)
    def add(self,f,g):
        return trim([self.F.add(f[i] if i<len(f) else 0,g[i] if i<len(g) else 0) for i in range(max(len(f),len(g)))])
    def sub(self,f,g):return self.add(f,[self.F.neg(c) for c in g])
    def mul(self,f,g):
        h=[0]*(len(f)+len(g)-1)
        for i,a in enumerate(f):
            for j,b in enumerate(g):h[i+j]=self.F.add(h[i+j],self.F.mul(a,b))
        return trim(h)
    def rem(self,f,g):
        f=list(f);iv=self.F.inv(g[-1])
        while f!=[0] and len(f)>=len(g):
            c=self.F.mul(f[-1],iv);d=len(f)-len(g)
            for i in range(len(g)):f[d+i]=self.F.sub(f[d+i],self.F.mul(c,g[i]))
            trim(f)
        return f
    def monic(self,f):return [self.F.div(c,f[-1]) for c in f]
    def gcd(self,f,g):
        while g!=[0]:f,g=g,self.rem(f,g)
        return self.monic(f)
    def power(self,f,n,m):
        r=[1]
        while n:
            if n&1:r=self.rem(self.mul(r,f),m)
            f=self.rem(self.mul(f,f),m);n>>=1
        return r
    def quotient(self,f,g):
        f=list(f);r=[0]*(len(f)-len(g)+1);iv=self.F.inv(g[-1])
        while len(f)>=len(g):
            c=self.F.mul(f[-1],iv);d=len(f)-len(g);r[d]=c
            for i in range(len(g)):f[d+i]=self.F.sub(f[d+i],self.F.mul(c,g[i]))
            trim(f)
            if f==[0]:break
        assert f==[0]
        return trim(r)
    def roots(self,f):
        f=self.monic(trim(list(f)))
        assert self.rem(self.sub(self.power([0,1],self.F.q,f),[0,1]),f)==[0], 'Polynomial does not split over Fp2'
        pending=[f];out=[]
        while pending:
            g=pending.pop()
            if len(g)==2:out.append(self.F.neg(g[0]));continue
            while True:
                h=[self.rng.randrange(self.F.q) for _ in range(len(g)-1)]
                d=self.gcd(g,self.sub(self.power(h,(self.F.q-1)//2,g),[1]))
                if 1<len(d)<len(g):break
            pending.extend([d,self.quotient(g,d)])
        assert len(out)==len(f)-1 and len(set(out))==len(out)
        return sorted(out)

def class_number(D):
    assert D<0 and D%4 in (0,1)
    forms=[]
    for a in range(1,math.isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D)%(4*a):continue
            c=(b*b-D)//(4*a)
            if a>c or math.gcd(a,math.gcd(b,c))!=1:continue
            if (abs(b)==a or a==c) and b<0:continue
            forms.append((a,b,c))
    return len(forms)

def trace_fp(p,a,b):
    squares={i*i%p for i in range(1,p)}
    return -sum(0 if (v:=(x*x*x+a*x+b)%p)==0 else 1 if v in squares else -1 for x in range(p))

def model(F,j):
    p=F.p
    if j==0:return 0,1
    if j==1728%p:return 1,0
    c=F.div(j,F.sub(1728%p,j))
    return F.mul(3,c),F.mul(2,c)

def jinvariant(F,a,b):
    a3=F.power(a,3);b2=F.mul(b,b)
    return F.div(F.mul(1728%F.p,F.mul(4,a3)),F.add(F.mul(4,a3),F.mul(27%F.p,b2)))

def seed(F):
    if F.p%4==3:return 1728%F.p
    for j in range(1,F.p):
        if j==1728%F.p:continue
        a,b=model(F,j)
        if trace_fp(F.p,a,b)==0:return j
    raise AssertionError('No rational supersingular j found')

def quotient3(F,a,b,r):
    # Velu's t=2(3r^2+a) and w=2(5r^3+3ar+2b).
    a1=F.sub(F.mul((-9)%F.p,a),F.mul(30%F.p,F.power(r,2)))
    b1=F.sub(F.sub(F.mul((-27)%F.p,b),F.mul(70%F.p,F.power(r,3))),F.mul(42%F.p,F.mul(a,r)))
    return a1,b1,jinvariant(F,a1,b1)

def verify_velu(P,a,b,r,aa,bb):
    F=P.F;p=F.p
    scale=lambda f,c:[F.mul(x,c) for x in f]
    deriv=lambda f:trim([F.mul(i%p,f[i]) for i in range(1,len(f))])
    d=[F.mul(r,r),F.mul((-2)%p,r),1]
    n=[F.add(F.add(F.mul((-2)%p,F.power(r,3)),F.mul(2,F.mul(a,r))),F.mul(4,b)),
       F.add(F.mul(7,F.mul(r,r)),F.mul(2,a)),F.mul((-2)%p,r),1]
    q=P.sub(P.mul(deriv(n),d),P.mul(n,deriv(d)))
    # x'=n/d, y'=y (n/d)' is an exact normalized isogeny identity.
    lhs=P.mul([b,a,0,1],P.mul(q,q))
    d2=P.mul(d,d);d3=P.mul(d2,d)
    rhs=P.mul(P.add(P.add(P.mul(P.mul(n,n),n),scale(P.mul(n,d2),aa)),scale(d3,bb)),d)
    assert lhs==rhs

def components(B):
    unseen=set(range(len(B)));out=[]
    while unseen:
        c={min(unseen)};stack=list(c)
        while stack:
            i=stack.pop()
            for j in unseen-c:
                if B[i][j] or B[j][i]:c.add(j);stack.append(j)
        unseen-=c;out.append(sorted(c))
    return sorted(out,key=lambda c:(-len(c),c))

def split_kernel_count2(n1,n2):
    # Direct enumeration in F_2^4 with symplectic pairs (e1,f1),(e2,f2).
    pair=lambda a,b:((a&1)*((b>>1)&1)+((a>>1)&1)*(b&1)+((a>>2)&1)*((b>>3)&1)+((a>>3)&1)*((b>>2)&1))%2
    planes={tuple(sorted((0,a,b,a^b))) for a in range(1,16) for b in range(a+1,16) if pair(a,b)==0}
    assert len(planes)==15
    def frob(v):
        return v ^ (((v>>1)&1) if n1==1 else 0) ^ ((4*((v>>3)&1)) if n2==1 else 0)
    stable=[L for L in planes if {frob(v) for v in L}==set(L)]
    product=[L for L in stable if len(set(L)&{0,1,2,3})==2 and len(set(L)&{0,4,8,12})==2]
    assert len(product)==n1*n2
    return len(stable),len(product)

def run(p):
    start=time.time()
    assert 512<=p<1024 and p%3==1
    assert all(p%d for d in range(2,math.isqrt(p)+1))
    F=Field(p);P=Polys(F);j0=seed(F);seen={j0};queue=deque([j0]);elliptic={}
    while queue:
        j=queue.popleft();a,b=model(F,j)
        # psi_3(x)=3x^4+6a x^2+12b x-a^2.
        poly=[F.neg(F.mul(a,a)),F.mul(12,b),F.mul(6,a),0,3]
        roots=P.roots(poly);assert len(roots)==4
        edges=[]
        for r in roots:
            aa,bb,jj=quotient3(F,a,b,r)
            verify_velu(P,a,b,r,aa,bb)
            edges.append({'kernel_x':r,'target_j':jj,'quotient_a':aa,'quotient_b':bb})
            if jj not in seen:seen.add(jj);queue.append(jj)
        elliptic[j]={'j':j,'a':a,'b':b,'kernels':edges}
    t=p//12+{1:0,5:1,7:1,11:2}[p%12]
    assert len(seen)==t,(p,len(seen),t)
    assert {F.frob(j) for j in seen}==seen
    rational=sorted(j for j in seen if j<p);s=len(rational)
    for j in rational:
        a,b=model(F,j);assert trace_fp(p,a,b)==0
    h=class_number(-p if p%4==3 else -4*p)
    predicted_s=h if p%8==7 else 2*h if p%8==3 else h//2
    assert s==predicted_s
    orbits=sorted({tuple(sorted({j,F.frob(j)})) for j in seen})
    oi={j:i for i,o in enumerate(orbits) for j in o}
    n=len(orbits);assert n==(s+t)//2
    B=[[0]*n for _ in range(n)]
    vertices=[]
    for i,o in enumerate(orbits):
        aorder=2 if len(o)==2 else 8 if o[0]==1728%p else 4
        vertices.append({'id':i,'j_orbit':list(o),'automorphism_order':aorder})
        for e in elliptic[o[0]]['kernels']:B[i][oi[e['target_j']]]+=1
    assert all(sum(r)==4 for r in B)
    for i in range(n):
        for j in range(n):assert B[i][j]*vertices[j]['automorphism_order']==B[j][i]*vertices[i]['automorphism_order'],(p,i,j)
    comps=components(B)
    # Enumerate actual elliptic Fp-forms and the split pairs, retaining twists.
    forms=[]
    for j in rational:
        a,b=model(F,j)
        candidates=[(a,b),(F.nu*a%p,0)] if j==1728%p else [(a,b),(a*F.nu**2%p,b*F.nu**3%p)]
        assert len(set(candidates))==2
        for aa,bb in candidates:
            assert trace_fp(p,aa,bb)==0
            roots2=[x for x in range(p) if (x**3+aa*x+bb)%p==0]
            assert len(roots2) in [1,3]
            forms.append({'j':j,'a':aa,'b':bb,'rational_2_kernels':len(roots2)})
    pairs=list(itertools.combinations_with_replacement(range(len(forms)),2))
    split_degrees=Counter();product2=0
    local2={(i,j):split_kernel_count2(i,j) for i in (1,3) for j in (1,3)}
    for i,j in pairs:
        ni,nj=forms[i]['rational_2_kernels'],forms[j]['rational_2_kernels']
        actual_degree,actual_product=local2[ni,nj]
        product2+=actual_product
        split_degrees[actual_degree]+=1
    assert len(pairs)==s*(2*s+1)
    data={'p':p,'bits':p.bit_length(),'p_mod_8':p%8,'field_relation':f'u^2={F.nu}',
          'field_encoding':'a+p*b denotes a+b*u, 0<=a,b<p','ell':3,
          'scope':'Complete closed Weil-restriction sector in the maximal block; not the full PPAS graph',
          'model_convention':'Auxiliary geometric models compute quotient j. Vertices represent unique maximal elliptic forms; generic quadratic twisting leaves the geometric kernel/codomain-j multiset unchanged.',
          'h_quadratic_field':h,'s':s,'t':t,'split_vertices':len(pairs),'restriction_vertices':n,
          'elliptic_seed_j':j0,'elliptic_models':list(elliptic.values()),'elliptic_fp_forms':forms,
          'vertices':vertices,'adjacency':B,'components':comps,
          'directed_kernels':4*n,'loops':sum(B[i][i] for i in range(n)),
          'split_ell2_degree_histogram':dict(split_degrees),'split_ell2_product_kernels':product2,
          'jacobian_asymptotic_main_term':f'{p*s}/3','jacobian_exact_count':None,
          'checks':['Geometric quotient calculation by all four roots of psi_3','Vélu quotient j and rational-map polynomial identity','Closure and geometric elliptic census','Frobenius orbit census','Trace-zero verification for every rational elliptic form','All four kernels per restriction vertex','Weighted polarized duality','Direct enumeration of all 15 symplectic planes for each elliptic two-torsion action type'],
          'seconds':round(time.time()-start,3)}
    (OUT/f'p{p}_ten_bit_restriction_ell3.json').write_text(json.dumps(data,indent=2),encoding='utf8')
    with (OUT/f'p{p}_ten_bit_restriction_ell3_adjacency.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.writer(stream);writer.writerow(['source/target']+list(range(n)))
        for i,row in enumerate(B):writer.writerow([i]+row)
    print(p, 's/t',s,t,'S/R',len(pairs),n,'edges',4*n,'components',list(map(len,comps)),'seconds',data['seconds'],flush=True)
    return {k:data[k] for k in ['p','bits','p_mod_8','h_quadratic_field','s','t','split_vertices','restriction_vertices','directed_kernels','loops','split_ell2_degree_histogram','split_ell2_product_kernels','jacobian_asymptotic_main_term','seconds']}|{'component_sizes':list(map(len,comps))}

if __name__=='__main__':
    primes=list(map(int,sys.argv[1:])) or [523,541,547,613,619,631,661,673,751,823,991]
    rows=[run(p) for p in primes]
    (OUT/'ten_bit_restriction_summary.json').write_text(json.dumps(rows,indent=2),encoding='utf8')

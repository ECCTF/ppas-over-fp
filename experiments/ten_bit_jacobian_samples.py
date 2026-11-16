"""Actual rational Richelot neighborhoods at 10-bit primes.

The curves are constructed with two trace-zero elliptic quotients.
Their Frobenius polynomial and every smooth Richelot quotient are checked
by direct Fp/Fp2 point counts and the Cartier--Manin matrix.  Vertices are
deduplicated by an exact rational binary-sextic isomorphism test, including
the square-scalar condition. This is a complete neighborhood of each
specified source, NOT an enumeration of every genus-two curve.

Requires NumPy and the adjacent ten_bit_experiments.py finite-field helpers.
Run with the bundled Python or Sage Python: python ten_bit_jacobian_samples.py.
"""
from pathlib import Path
import itertools, json, time
import numpy as np
from ten_bit_experiments import Field, Polys, trace_fp, jinvariant

OUT=Path(__file__).resolve().parent/'results'

def matchings(items):
    if not items:
        yield [];return
    a=items[0]
    for b in items[1:]:
        for rest in matchings([i for i in items if i not in (a,b)]):
            yield [(a,b)]+rest

class Experiment:
    def __init__(self,p):
        self.p=p;self.F=Field(p);self.P=Polys(self.F)
        self.chi=[-1]*p;self.chi[0]=0
        for i in range(1,p):self.chi[i*i%p]=1
        self.cache={}
        self.root_cache={}
    def scale(self,f,c):return [self.F.mul(x,c) for x in f]
    def derivative(self,f):return [self.F.mul(i%self.p,f[i]) for i in range(1,len(f))] or [0]
    def value(self,f,x):
        y=0
        for c in reversed(f):y=self.F.add(self.F.mul(y,x),c)
        return y
    def points(self,f):
        p=self.p;F=self.F;chi=self.chi
        assert all(0<=c<p for c in f)
        n1=p+1+(chi[f[-1]] if len(f)==7 else 0)
        n1+=sum(chi[sum(c*pow(x,i,p) for i,c in enumerate(f))%p] for x in range(p))
        # chi_(p^2)(z)=chi_p(Norm(z)); pair conjugate x values.
        total=sum(0 if self.value(f,x)==0 else 1 for x in range(p))
        for b in range(1,(p+1)//2):
            for a in range(p):
                y=self.value(f,a+p*b);u,v=y%p,y//p
                total+=2*chi[(u*u-F.nu*v*v)%p]
        n2=p*p+1+(1 if len(f)==7 else 0)+total
        return n1,n2
    def cm(self,f):
        # All integer convolution coefficients fit int64: at most 2p terms
        # each bounded by (p-1)^2, for the 10-bit parameters used here.
        p=self.p;n=(p-1)//2;f=np.array(f,dtype=np.int64);r=np.array([1],dtype=np.int64)
        while n:
            if n&1:r=np.convolve(r,f)[:2*p] % p
            f=np.convolve(f,f)[:2*p] % p;n>>=1
        return [int(r[i]) if i<len(r) else 0 for i in [p-1,p-2,2*p-1,2*p-2]]
    def verify(self,f):
        key=tuple(f)
        if key not in self.cache:
            assert self.P.gcd(f,self.derivative(f))==[1]
            m=self.cm(f);assert m==[0,0,0,0],(self.p,f,m)
            n1,n2=self.points(f)
            assert (n1,n2)==(self.p+1,self.p*self.p+1+4*self.p),(self.p,f,n1,n2)
            self.cache[key]={'n1':n1,'n2':n2,'cartier_manin':m,'weil_polynomial':[self.p*self.p,0,2*self.p,0,1]}
        return self.cache[key]
    def roots(self,f):
        if tuple(f) in self.root_cache:return self.root_cache[tuple(f)]
        out=self.P.roots(f)
        if len(f)==6:out.append(None)
        assert len(out)==6
        self.root_cache[tuple(f)]=out
        return out
    def triple(self,a,b,c):
        F=self.F
        def point(z):return (1,0) if z is None else (z,1)
        aa,bb,cc=map(point,(a,b,c))
        la=(aa[1],F.neg(aa[0]));lb=(bb[1],F.neg(bb[0]))
        va=F.add(F.mul(la[0],cc[0]),F.mul(la[1],cc[1]))
        vb=F.add(F.mul(lb[0],cc[0]),F.mul(lb[1],cc[1]))
        return (F.mul(lb[0],va),F.mul(lb[1],va),F.mul(la[0],vb),F.mul(la[1],vb))
    def matmul(self,A,B):
        F=self.F
        return tuple(F.add(F.mul(A[2*i],B[j]),F.mul(A[2*i+1],B[2+j])) for i in range(2) for j in range(2))
    def inverse(self,A):
        F=self.F;a,b,c,d=A;z=F.inv(F.sub(F.mul(a,d),F.mul(b,c)))
        return tuple(F.mul(z,x) for x in (d,F.neg(b),F.neg(c),a))
    def transform(self,g,M):
        a,b,c,d=M;P=self.P
        aa=[[1]];bb=[[1]]
        for _ in range(6):aa.append(P.mul(aa[-1],[b,a]));bb.append(P.mul(bb[-1],[d,c]))
        h=[0]
        for i,v in enumerate(g):h=P.add(h,self.scale(P.mul(aa[i],bb[6-i]),v))
        return h
    def isomorphisms(self,f,g,stop=False):
        F=self.F;p=self.p;source=self.triple(*self.roots(f)[:3]);count=0
        for target in itertools.permutations(self.roots(g),3):
            M=self.matmul(self.inverse(self.triple(*target)),source)
            M=tuple(F.div(x,next(y for y in M if y)) for x in M)
            if any(x>=p for x in M):continue
            h=self.transform(g,M)
            if len(h)!=len(f):continue
            c=F.div(h[-1],f[-1])
            if c>=p or self.chi[c]!=1:continue
            if h==self.scale(f,c):
                count+=2
                if stop:return count
        return count
    def short(self,g):
        F=self.F;p=self.p;d,c,b,a=g
        A=F.sub(F.mul(a,c),F.div(F.mul(b,b),3))
        B=F.add(F.sub(F.mul(F.mul(a,a),d),F.div(F.mul(F.mul(a,b),c),3)),F.div(F.mul(2,F.power(b,3)),27%p))
        return A,B
    def canonical_elliptic(self,a,b):
        p=self.p
        return min(((a*pow(u,4,p))%p,(b*pow(u,6,p))%p) for u in range(1,p))
    def source(self):
        p=self.p
        # C_a : y²=(x²-a)^3+(x²-a), first elliptic quotient y²=t³+t.
        for a in range(1,p):
            if self.chi[a]!=-1:continue
            f=[(-a*a*a-a)%p,0,(3*a*a+1)%p,0,(-3*a)%p,0,1]
            if f[0]==0:continue
            A,B=self.short([1,f[4],f[2],f[0]])
            if trace_fp(p,A,B)==0:
                assert trace_fp(p,1,0)==0
                return a,f,(A,B)
        raise AssertionError('No source found')
    def det3(self,rows):
        F=self.F;r=0
        for perm in itertools.permutations(range(3)):
            z=1
            for i,j in enumerate(perm):z=F.mul(z,rows[i][j])
            odd=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))%2
            r=F.sub(r,z) if odd else F.add(r,z)
        return r
    def count_reverse(self,g,f):
        F=self.F;P=self.P;roots=self.roots(g);count=0
        rho=[roots.index(None if x is None else F.frob(x)) for x in roots]
        source_fixed=sum(x is None or x<self.p for x in self.roots(f))
        for m in matchings(list(range(6))):
            if {frozenset((rho[i],rho[j])) for i,j in m}!={frozenset(q) for q in m}:continue
            qs=[]
            for pair in m:
                q=[1]
                for i in pair:
                    if roots[i] is not None:q=P.mul(q,[F.neg(roots[i]),1])
                qs.append(q)
            qs[0]=self.scale(qs[0],g[-1])
            delta=self.det3([[q[i] if i<len(q) else 0 for i in (2,1,0)] for q in qs])
            if not delta:continue
            def bracket(u,v):return P.sub(P.mul(self.derivative(u),v),P.mul(u,self.derivative(v)))
            h=self.scale(P.mul(P.mul(bracket(qs[1],qs[2]),bracket(qs[2],qs[0])),bracket(qs[0],qs[1])),F.neg(F.inv(delta)))
            assert all(c<self.p for c in h)
            if sum(x is None or x<self.p for x in self.roots(h))!=source_fixed:continue
            if self.isomorphisms(f,h,stop=True):count+=1
        return count
    def decompose(self,qs):
        F=self.F;P=self.P;p=self.p;A,B,C=qs
        coeff=lambda q,i:q[i] if i<len(q) else 0
        disc=P.sub(P.mul([coeff(A,1),coeff(B,1)],[coeff(A,1),coeff(B,1)]),
                   self.scale(P.mul([coeff(A,2),coeff(B,2)],[coeff(A,0),coeff(B,0)]),4))
        ts=P.roots(disc)
        fixed=[]
        for t in ts:
            den=F.add(coeff(A,2),F.mul(t,coeff(B,2)))
            fixed.append(F.div(F.neg(F.add(coeff(A,1),F.mul(t,coeff(B,1)))),F.mul(2,den)) if den else None)
        if len(disc)<3:fixed.append(F.div(F.neg(coeff(B,1)),F.mul(2,coeff(B,2))) if coeff(B,2) else None)
        assert len(fixed)==2 and fixed[0]!=fixed[1]
        ls=[[1] if z is None else [F.neg(z),1] for z in fixed]
        ss=[P.mul(l,l)+[0]*3 for l in ls]
        ij=next((i,j) for i,j in itertools.combinations(range(3),2) if F.sub(F.mul(ss[0][i],ss[1][j]),F.mul(ss[1][i],ss[0][j])))
        i,j=ij;M=(ss[0][i],ss[1][i],ss[0][j],ss[1][j]);Mi=self.inverse(M)
        uv=[]
        for q in qs:
            u=F.add(F.mul(Mi[0],coeff(q,i)),F.mul(Mi[1],coeff(q,j)))
            v=F.add(F.mul(Mi[2],coeff(q,i)),F.mul(Mi[3],coeff(q,j)))
            assert all(F.add(F.mul(u,ss[0][k]),F.mul(v,ss[1][k]))==coeff(q,k) for k in range(3))
            uv.append((u,v))
        g1=g2=[1]
        for u,v in uv:g1=P.mul(g1,[v,u]);g2=P.mul(g2,[u,v])
        if all(z is None or z<p for z in fixed):
            assert all(c<p for c in g1+g2)
            es=[self.short(g) for g in (g1,g2)]
            assert all(trace_fp(p,a,b)==0 for a,b in es)
            pair=sorted(self.canonical_elliptic(a,b) for a,b in es)
            return {'type':'split','elliptic_models':pair,'elliptic_traces':[0,0],
                    'aut':8 if pair[0]==pair[1] else 4}
        assert F.frob(fixed[0])==fixed[1]
        a,b=self.short(g1);j=jinvariant(F,a,b)
        # The parent isogeny implies the maximal restriction block. Its
        # elliptic geometric model is recorded, and the j orbit identifies
        # the unique maximal elliptic form used for the restriction vertex.
        aut=2 if j>=p else 12 if j==0 else 8 if j==1728%p else 4
        return {'type':'restriction','j_orbit':sorted({j,F.frob(j)}),'auxiliary_elliptic_model':[a,b],'aut':aut}
    def run(self):
        started=time.time();p=self.p;F=self.F;P=self.P
        a,f,second=self.source();check=self.verify(f);roots=self.roots(f)
        vertices=[{'label':'J0','type':'Jacobian','coefficients':f,**check,'aut':self.isomorphisms(f,f)}]
        permutation=[roots.index(None if x is None else F.frob(x)) for x in roots]
        kernels=[]
        for matching in matchings(list(range(6))):
            if {frozenset((permutation[i],permutation[j])) for i,j in matching}!={frozenset(q) for q in matching}:continue
            qs=[]
            for pair in matching:
                q=[1]
                for i in pair:
                    if roots[i] is not None:q=P.mul(q,[F.neg(roots[i]),1])
                qs.append(q)
            qs[0]=self.scale(qs[0],f[-1])
            delta=self.det3([[q[i] if i<len(q) else 0 for i in (2,1,0)] for q in qs])
            rec={'matching':matching,'delta':delta,'quadratics':qs}
            if delta:
                def bracket(u,v):return P.sub(P.mul(self.derivative(u),v),P.mul(u,self.derivative(v)))
                # Descending coefficient determinant: the rational quotient is
                # -[B,C][C,A][A,B]/delta. The opposite sign twists the codomain.
                g=self.scale(P.mul(P.mul(bracket(qs[1],qs[2]),bracket(qs[2],qs[0])),bracket(qs[0],qs[1])),F.neg(F.inv(delta)))
                assert all(c<p for c in g),(p,g)
                assert len(g) in (6,7)
                verification=self.verify(g)
                candidates=[i for i,v in enumerate(vertices) if v['type']=='Jacobian' and self.isomorphisms(v['coefficients'],g,stop=True)]
                assert len(candidates)<=1
                if candidates:target=candidates[0]
                else:
                    target=len(vertices)
                    vertices.append({'label':f'J{sum(v["type"]=="Jacobian" for v in vertices)}','type':'Jacobian','coefficients':g,**verification,'aut':self.isomorphisms(g,g)})
                rec['quotient_coefficients']=g
            else:
                v=self.decompose(qs)
                same=lambda old:old['type']==v['type'] and (old.get('elliptic_models')==v.get('elliptic_models') if v['type']=='split' else old['j_orbit']==v['j_orbit'])
                ids=[i for i,old in enumerate(vertices) if same(old)]
                if ids:target=ids[0]
                else:
                    target=len(vertices);v['label']=('S' if v['type']=='split' else 'R')+str(sum(old['type']==v['type'] for old in vertices));vertices.append(v)
            rec['source']=0;rec['target']=target;kernels.append(rec)
        assert len(kernels)==7,(p,permutation,len(kernels))
        for vertex in vertices:
            if vertex['type']=='Jacobian':
                rr=self.roots(vertex['coefficients']);rho=[rr.index(None if x is None else F.frob(x)) for x in rr]
                vertex['rational_outdegree']=sum({frozenset((rho[i],rho[j])) for i,j in m}=={frozenset(q) for q in m} for m in matchings(list(range(6))))
        edges=[]
        for target in sorted({k['target'] for k in kernels}):
            count=sum(k['target']==target for k in kernels)
            num=count*vertices[target]['aut'];den=vertices[0]['aut']
            assert num%den==0
            edge={'source':0,'target':target,'forward_kernels':count,'reverse_kernels_by_duality':num//den}
            if vertices[target]['type']=='Jacobian':
                reverse=self.count_reverse(vertices[target]['coefficients'],f)
                assert reverse==num//den
                edge['reverse_kernels_directly_verified']=reverse
            edges.append(edge)
        result={'p':p,'bits':p.bit_length(),'source_parameter_a':a,'source_elliptic_quotients':[[1,0],list(second)],'extension_nonresidue':F.nu,'roots':roots,'frobenius_permutation':permutation,'vertices':vertices,'kernels':kernels,'edges':edges,'source_outdegree':len(kernels),'scope':'All rational (2,2) kernels from one specified Jacobian, with reverse multiplicities from polarized duality; other outgoing rows are not enumerated.','isomorphism_test':'All 120 branch-triple maps, rational matrix and square scalar in Fp','seconds':round(time.time()-started,3)}
        (OUT/f'p{p}_jacobian_neighborhood.json').write_text(json.dumps(result,indent=2)+'\n')
        print(p,'source a=',a,'vertices',[(v['label'],v['type']) for v in vertices],'targets',[k['target'] for k in kernels],'seconds',result['seconds'],flush=True)
        return result

if __name__=='__main__':
    small=Experiment(7)
    f=[0,5,4,1,1,0,1];g=[1,0,4,0,1,0,1]
    assert small.isomorphisms(f,f)==4 and small.isomorphisms(g,g)==8
    assert small.isomorphisms(f,g)==0
    print('Exact rational isomorphism routine agrees with both F7 regression examples.',flush=True)
    for p in (719,751,823):Experiment(p).run()

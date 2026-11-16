"""Exact rational Richelot BFS, retaining rational twists.

Reproducible exact maximal-block census. The canonical hash normalizes every ordered
branch triple to infinity, zero, one. It identifies geometric branch sets.
The ensuing exact binary-sextic comparison retains only Fp matrices and
square multipliers, so the hash NEVER identifies rational twists itself.
"""
from pathlib import Path
import argparse
from array import array
from collections import Counter
import itertools
import json
import sys
import time

BASE = Path(__file__).resolve().parent
OUT = BASE / 'results'
sys.path.insert(0, str(BASE))
from ten_bit_jacobian_samples import Experiment, matchings


class Census(Experiment):
    def __init__(self, p, richelot_sign=-1):
        super().__init__(p)
        self.richelot_sign=richelot_sign
        self.started = time.time()
        self.inv_fp = [0] + [pow(x, -1, p) for x in range(1, p)]
        self.F.inv = self.fast_inv
        self.sqrt = array('i', [-1]) * (p*p)
        for x in range(p*p):
            self.sqrt[self.F.mul(x, x)] = x
        self.old_roots = self.P.roots
        self.P.roots = self.fast_roots
        self.vertices = []
        self.canonical_data = []
        self.buckets = {}
        self.literal = {}
        self.edges = Counter()
        self.decomposed = Counter()
        self.decomposed_records = {}
        self.all_matchings = list(matchings(list(range(6))))
        self.edge_count = 0

    def fast_inv(self, x):
        p = self.p
        a, b = x % p, x // p
        norm = (a*a - self.F.nu*b*b) % p
        assert norm
        inv = self.inv_fp[norm]
        return (a*inv % p) + p*(-b*inv % p)

    def fast_roots(self, f):
        while len(f)>1 and not f[-1]:
            f=f[:-1]
        if len(f)==1:
            assert f[0]
            return []
        F=self.F
        if len(f)==2:
            return [F.div(F.neg(f[0]),f[1])]
        if len(f)==3:
            c,b,a=f
            disc=F.sub(F.mul(b,b),F.mul(4,F.mul(a,c)))
            root=self.sqrt[disc]
            assert root>=0, (self.p, f, disc)
            den=F.inv(F.mul(2,a))
            return sorted({F.mul(F.add(F.neg(b),root),den),
                           F.mul(F.sub(F.neg(b),root),den)})
        return self.old_roots(f)

    def normalize(self, f):
        f=list(f)
        while not f[-1]:f.pop()
        assert len(f) in (6,7) and all(c<self.p for c in f)
        lead=f[-1]
        desired=1 if self.chi[lead]==1 else self.F.nu
        scale=desired*pow(lead,-1,self.p)%self.p
        assert self.chi[scale]==1
        return tuple(c*scale%self.p for c in f)

    def canonical(self, roots):
        """Exact geometric key and all maps attaining the least key."""
        F=self.F
        def det(x,y):
            if x is None:return 1
            if y is None:return F.neg(1)
            return F.sub(x,y)
        ds=[[det(x,y) for y in roots] for x in roots]
        best=None
        minima=[]
        for a in range(6):
            for b in range(6):
                if a==b:continue
                remain=[i for i in range(6) if i!=a and i!=b]
                ratios={i:F.div(ds[i][b],ds[i][a]) for i in remain}
                for c in remain:
                    scale=F.inv(ratios[c])
                    key=tuple(sorted(F.mul(ratios[d],scale) for d in remain if d!=c))
                    if best is None or key<best:
                        best=key;minima=[(a,b,c)]
                    elif key==best:
                        minima.append((a,b,c))
        mats=[self.triple(*(roots[i] for i in triple)) for triple in minima]
        return best,mats

    def rational_isomorphisms(self, f, g, canonical_f, canonical_g, stop=False):
        """Maps x on f to x on g. Count both hyperelliptic lifts."""
        F=self.F
        if canonical_f[0]!=canonical_g[0]:return 0
        source=canonical_f[1][0]
        count=0
        for target in canonical_g[1]:
            M=self.matmul(self.inverse(target),source)
            scalar=F.inv(next(x for x in M if x))
            M=tuple(F.mul(x,scalar) for x in M)
            if any(x>=self.p for x in M):continue
            h=self.transform(g,M)
            if len(h)!=len(f):continue
            scalar=F.div(h[-1],f[-1])
            if scalar>=self.p or self.chi[scalar]!=1:continue
            if h==self.scale(f,scalar):
                count+=2
                if stop:return count
        return count

    def insert(self, f, roots=None):
        f=self.normalize(f)
        if f in self.literal:return self.literal[f]
        if roots is None:roots=self.roots(list(f))
        roots=sorted((x for x in roots if x is not None)) + ([None] if None in roots else [])
        assert len(roots)==6 and len(set(roots))==6
        for root in roots:
            assert root is None or self.value(f,root)==0
        canon=self.canonical(roots)
        for i in self.buckets.get(canon[0],[]):
            if self.rational_isomorphisms(f,self.vertices[i]['coefficients'],canon,self.canonical_data[i],True):
                self.literal[f]=i
                return i
        i=len(self.vertices)
        rho=[roots.index(None if x is None else self.F.frob(x)) for x in roots]
        assert all(rho[rho[j]]==j for j in range(6))
        fixed=sum(rho[j]==j for j in range(6))
        aut=self.rational_isomorphisms(f,f,canon,canon)
        assert aut>=2
        self.vertices.append(dict(id=i,label=f'J{i}',type='Jacobian',coefficients=list(f),roots=roots,
                                  aut=aut,Frobenius2type=f'1^{fixed}2^{(6-fixed)//2}',
                                  frobenius_permutation=rho,geometric_reduced_aut=len(canon[1])))
        self.canonical_data.append(canon)
        self.buckets.setdefault(canon[0],[]).append(i)
        self.literal[f]=i
        return i

    def expand(self,i):
        F=self.F;P=self.P
        v=self.vertices[i];f=v['coefficients'];roots=v['roots'];rho=v['frobenius_permutation']
        degree=0
        for matching in self.all_matchings:
            if {frozenset((rho[a],rho[b])) for a,b in matching}!={frozenset(pair) for pair in matching}:continue
            degree+=1
            qs=[]
            for pair in matching:
                q=[1]
                for j in pair:
                    if roots[j] is not None:q=P.mul(q,[F.neg(roots[j]),1])
                qs.append(q)
            qs[0]=self.scale(qs[0],f[-1])
            delta=self.det3([[q[j] if j<len(q) else 0 for j in (2,1,0)] for q in qs])
            if not delta:
                target=self.decompose(qs)
                key=json.dumps(target,sort_keys=True)
                # Auxiliary elliptic equations are not part of the target label.
                if target['type']=='restriction':
                    key=('restriction',tuple(target['j_orbit']))
                else:key=('split',tuple(map(tuple,target['elliptic_models'])))
                self.decomposed_records[key]=target
                self.decomposed[(i,key)]+=1
                continue
            def bracket(u,w):return P.sub(P.mul(self.derivative(u),w),P.mul(u,self.derivative(w)))
            brackets=[bracket(qs[1],qs[2]),bracket(qs[2],qs[0]),bracket(qs[0],qs[1])]
            h=self.scale(P.mul(P.mul(brackets[0],brackets[1]),brackets[2]),F.mul(self.richelot_sign%self.p,F.inv(delta)))
            newroots=[]
            for q in brackets:
                newroots+=self.fast_roots(q)
                if len(q)==2:newroots.append(None)
            j=self.insert(h,newroots)
            self.edges[(i,j)]+=1
        expected={6:15,4:3,2:3,0:7}[sum(x is None or x<self.p for x in roots)]
        assert degree==expected,(i,degree,expected)
        v['degree']=degree
        self.edge_count+=degree

    def run(self,seeds,output,limit=None):
        seed_indices=[]
        for record in seeds['kernels']:
            seed_indices.append(self.insert(record['coefficients']))
        print(f'p={self.p}: {len(seed_indices)} split/J kernels, {len(self.vertices)} distinct seed Jacobians, t={time.time()-self.started:.2f}s',flush=True)
        done=0
        while done<len(self.vertices):
            self.expand(done)
            done+=1
            if done%100==0 or done==len(self.vertices):
                print(f'p={self.p}: processed={done}, discovered={len(self.vertices)}, geometric={len(self.buckets)}, kernels={self.edge_count}, t={time.time()-self.started:.2f}s',flush=True)
            if limit and done>=limit:break
        if done==len(self.vertices):
            for (i,j),n in self.edges.items():
                assert n*self.vertices[j]['aut']==self.edges[(j,i)]*self.vertices[i]['aut'],(i,j,n,self.edges[(j,i)])
        result=dict(p=self.p,richelot_sign=self.richelot_sign,field={'p':self.p,'nu':self.F.nu,'encoding':'a+p*b with u^2=nu'},
                    complete_jacobian_closure=done==len(self.vertices),processed=done,
                    elapsed_seconds=time.time()-self.started,vertices=self.vertices,
                    edges=[{'source':i,'target':j,'multiplicity':n} for (i,j),n in sorted(self.edges.items())],
                    decomposed_edges=[dict(source=i,multiplicity=n,**self.decomposed_records[key]) for (i,key),n in self.decomposed.items()],
                    seed_indices=seed_indices,geometric_jacobian_classes=len(self.buckets),
                    rational_jacobian_classes=len(self.vertices),jacobian_kernel_count=self.edge_count)
        output.write_text(json.dumps(result,indent=2),encoding='utf8')
        print(f'Saved {output}',flush=True)
        return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('p',type=int)
    parser.add_argument('--limit',type=int)
    args=parser.parse_args()
    path=OUT/f'p{args.p}_split_gluing_seeds.json'
    seeds=json.loads(path.read_text(encoding='utf8'))
    Census(args.p).run(seeds,OUT/f'p{args.p}_rational_census.json',args.limit)

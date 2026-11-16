"""Explicit exchanged-factor descent of HLP gluing, retaining the prime-field quotient form."""
import itertools,json,sys,time
from pathlib import Path
from collections import Counter
BASE=Path(__file__).resolve().parent
OUT=BASE/'results'
sys.path.insert(0,str(BASE))
from ten_bit_jacobian_samples import Experiment
from ten_bit_experiments import jinvariant
from split_gluing_seeds import glue,elliptic_data


def restriction_models(p,E):
    F=E.F;P=E.P
    fn=OUT/f'p{p}_restrictions.json'
    if not fn.exists():
        out=[];used=set()
        for e in elliptic_data(p):
            if e['j'] in used:continue
            used.add(e['j']);out.append((e['a'],e['b'],[e['j']]))
        return out
    d=json.loads(fn.read_text());z=P.roots(d['field_modulus'])[0]
    def element(coeff):
        s=0
        for c in reversed(coeff):s=F.add(F.mul(s,z),c)
        return s
    out=[];seen=set()
    for e in d['elliptic_curves']:
        if e['trace']!=-2*p:continue
        a,b=element(e['a4']),element(e['a6'])
        j=jinvariant(F,a,b);orb=tuple(sorted({j,F.frob(j)}))
        if orb in seen:continue
        seen.add(orb);out.append((a,b,list(orb)))
    return out


def generate_restriction_seeds(p):
    E=Experiment(p);F=E.F;P=E.P;start=time.time()
    models=restriction_models(p,E)
    # In the descent coordinate, x=(z-u)/(z+u), u^p=-u.
    u=p;M=(1,F.neg(u),1,u)
    result=dict(p=p,extension_nonresidue=F.nu,kernels=[],decomposed_graph_kernels=[],models=[])
    for a,b,orb in models:
        alpha=P.roots([b,a,0,1])
        row=dict(elliptic_coefficients=[a,b],j_orbit=orb,Jacobian_kernels=0,decomposed_graph_kernels=0)
        for perm in itertools.permutations(range(3)):
            if any(perm[perm[i]]!=i for i in range(3)):continue
            beta=[F.frob(alpha[perm[i]]) for i in range(3)]
            h=glue(F,P,alpha,beta,rational=False)
            rec=dict(j_orbit=orb,root_permutation=list(perm))
            if h is None:
                result['decomposed_graph_kernels'].append(rec);row['decomposed_graph_kernels']+=1
                continue
            assert [F.frob(c) for c in h]==list(reversed(h))
            g=E.transform(h,M)
            assert all(c<p for c in g)
            assert len(g) in (6,7)
            lead=g[-1];desired=1 if E.chi[lead]==1 else F.nu
            scale=desired*pow(lead,-1,p)%p
            assert E.chi[scale]==1
            g=[c*scale%p for c in g]
            assert E.cm(g)==[0,0,0,0]
            n=p+1+(E.chi[g[-1]] if len(g)==7 else 0)
            n+=sum(E.chi[sum(c*pow(x,i,p) for i,c in enumerate(g))%p] for x in range(p))
            assert n==p+1
            rec['coefficients']=g+[0]*(7-len(g))
            result['kernels'].append(rec);row['Jacobian_kernels']+=1
        assert row['Jacobian_kernels']+row['decomposed_graph_kernels']==4
        result['models'].append(row)
    result['seconds']=round(time.time()-start,3)
    return result


if __name__=='__main__':
    for p in map(int,sys.argv[1:] or [7,31,431]):
        d=generate_restriction_seeds(p)
        out=OUT/f'p{p}_restriction_gluing_seeds.json';out.write_text(json.dumps(d,indent=2)+'\n')
        print(p,len(d['kernels']),'smooth',len(d['decomposed_graph_kernels']),'decomposed',d['seconds'],'seconds',flush=True)

"""Independent verification of all split/restriction Jacobian boundaries.

Compare explicit Howe--Leprevost--Poonen gluing, including exchanged
descent for restrictions, against reversed degenerate Richelot quotients.
Also recover the actual two elliptic Fp quotient forms from every split
seed. Thus a scalar error that changes a rational twist is detected.

Run after split_gluing_seeds.py, restriction_gluing_seeds.py and
rational_census.py. Outputs results/p{p}_gluing_verification.json.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import json
from pathlib import Path
import time

from rational_census import Census

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'


def verify(p):
    started=time.time()
    E=Census(p);F=E.F;P=E.P
    splits=json.loads((OUT/f'p{p}_split_gluing_seeds.json').read_text())
    restrictions=json.loads((OUT/f'p{p}_restriction_gluing_seeds.json').read_text())
    census=json.loads((OUT/f'p{p}_rational_census.json').read_text())
    elliptics=splits['elliptic_classes']
    elliptic_ids={E.canonical_elliptic(v['a'],v['b']):i for i,v in enumerate(elliptics)}
    assert len(elliptic_ids)==len(elliptics)
    roots=[P.roots([v['b'],v['a'],0,1]) for v in elliptics]
    canons=[];buckets=defaultdict(list)
    for v in census['vertices']:
        c=E.canonical(v['roots']);canons.append(c)
        buckets[c[0]].append(v['id'])

    def identify(coefficients):
        f=list(E.normalize(coefficients));can=E.canonical(E.roots(f))
        hits=[i for i in buckets[can[0]]
              if E.rational_isomorphisms(f,census['vertices'][i]['coefficients'],can,canons[i],True)]
        assert len(hits)==1,(p,f,hits)
        return hits[0]

    # A separate description of the two degree-two elliptic quotients
    # verifies the Fp twist in every HLP split gluing.
    split_actual=Counter();restriction_actual=Counter()
    for k in splits['kernels']:
        f=k['coefficients'];assert len(f)==7
        assert E.cm(f)==[0,0,0,0]
        a=E.short([f[0],f[2],f[4],f[6]])
        b=E.short([f[6],f[4],f[2],f[0]])
        recovered=sorted([elliptic_ids[E.canonical_elliptic(*a)],
                          elliptic_ids[E.canonical_elliptic(*b)]])
        assert recovered==k['split_ids'],(p,k,recovered)
        split_actual[tuple(recovered),identify(f)]+=1
    for k in restrictions['kernels']:
        f=list(E.normalize(k['coefficients']))
        assert E.cm(f)==[0,0,0,0]
        restriction_actual[tuple(k['j_orbit']),identify(f)]+=1

    # Duality converts the independently computed J->S/R Richelot rows
    # into expected S/R->J entries. Compare every ordered entry.
    split_expected=Counter();restriction_expected=Counter()
    for edge in census['decomposed_edges']:
        j=edge['source'];aut_j=census['vertices'][j]['aut']
        reverse=Fraction(edge['multiplicity']*edge['aut'],aut_j)
        assert reverse.denominator==1
        if edge['type']=='split':
            pair=tuple(sorted(elliptic_ids[tuple(v)] for v in edge['elliptic_models']))
            split_expected[pair,j]+=reverse
        else:
            restriction_expected[tuple(edge['j_orbit']),j]+=reverse
    assert split_actual==split_expected,(p,'split boundary differs')
    assert restriction_actual==restriction_expected,(p,'restriction boundary differs')

    # Independently classify each affine graph kernel on a split source.
    # x-scalar square => rational elliptic isomorphism => split quotient.
    # x-scalar nonsquare => Frobenius negates it => restriction quotient.
    affine_records=[];reverse_to_split=defaultdict(Fraction)
    for rec in splits['decomposed_graph_kernels']:
        i,j=rec['split_ids'];alpha=roots[i]
        beta=[roots[j][r] for r in rec['root_permutation']]
        scale=F.div(F.sub(beta[1],beta[0]),F.sub(alpha[1],alpha[0]))
        shift=F.sub(beta[0],F.mul(scale,alpha[0]))
        assert 0<scale<p and shift==0
        assert all(beta[r]==F.mul(scale,alpha[r]) for r in range(3))
        assert elliptics[j]['a']==scale*scale*elliptics[i]['a']%p
        assert elliptics[j]['b']==scale*scale*scale*elliptics[i]['b']%p
        row=dict(source_split_ids=[i,j],x_scalar=scale,
                 root_permutation=rec['root_permutation'])
        if E.chi[scale]==1:
            assert i==j
            row.update(type='split',target_split_ids=[j,j])
        else:
            jj=elliptics[j]['j'];aut_r=12 if jj==0 else 8 if jj==1728%p else 4
            reverse_to_split[jj]+=Fraction(aut_r,8 if i==j else 4)
            row.update(type='restriction',j_orbit=[jj])
        affine_records.append(row)
    for jj in sorted({e['j'] for e in elliptics}):
        assert reverse_to_split[jj]==(3 if jj==0 else 2 if jj==1728%p else 1)

    # Recheck actual Fp and Fp2 point counts for samples throughout each
    # independent gluing list. All other maximality claims follow from
    # the explicit rational polarized quotient construction.
    point_counts=[]
    for kind,data in [('split',splits),('restriction',restrictions)]:
        n=len(data['kernels'])
        indices=sorted({0,n//4,n//2,3*n//4,n-1}) if n else []
        for index in indices:
            f=list(E.normalize(data['kernels'][index]['coefficients']))
            check=E.verify(f)
            point_counts.append(dict(source_type=kind,kernel_index=index,**check))
    report=dict(p=p,field={'p':p,'nu':F.nu},
                split_Jacobian_kernels=len(splits['kernels']),
                restriction_Jacobian_kernels=len(restrictions['kernels']),
                split_boundary_ordered_pairs=len(split_actual),
                restriction_boundary_ordered_pairs=len(restriction_actual),
                recovered_elliptic_quotient_pairs=len(splits['kernels']),
                vanishing_Cartier_Manin_matrices=len(splits['kernels'])+len(restrictions['kernels']),
                split_decomposed_graph_kernels=affine_records,
                restriction_reverse_split_totals={str(j):int(n) for j,n in reverse_to_split.items()},
                independent_point_counts=point_counts,
                checks=['exact Fp elliptic quotients of every split gluing',
                        'every split and restriction HLP kernel matched by exact Fp curve isomorphism',
                        'every explicit boundary entry agrees with reversed Richelot quotients',
                        'all affine split graph kernels classified by the rational square test',
                        'all seed Cartier-Manin matrices vanish',
                        'sample Fp and Fp2 point counts are maximal'],
                seconds=round(time.time()-started,3))
    fn=OUT/f'p{p}_gluing_verification.json';fn.write_text(json.dumps(report,indent=2)+'\n')
    print(p,'all',len(splits['kernels']),'split and',len(restrictions['kernels']),
          'restriction HLP kernels independently verified;',report['seconds'],'seconds',flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('primes',type=int,nargs='*',default=[431])
    args=parser.parse_args()
    for p in args.primes:verify(p)

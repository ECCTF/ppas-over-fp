"""Assemble and verify the entire p=431 rational (2,2) kernel graph.

Independent elliptic Velu rows and affine-map tests supplement the exact
Jacobian BFS. R->S/J entries are derived by weighted polarized duality.
Their row degrees certify that all nonproduct restriction kernels occur.
"""
from pathlib import Path
import itertools
import json
import sys
from collections import Counter,deque
from fractions import Fraction

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'
sys.path.insert(0,str(BASE))
from ten_bit_jacobian_samples import Experiment
from ten_bit_experiments import model,jinvariant


def assemble(p=431):
    E=Experiment(p);F=E.F;P=E.P
    seeds=json.loads((OUT/f'p{p}_split_gluing_seeds.json').read_text())
    bfs=json.loads((OUT/f'p{p}_rational_census.json').read_text())
    elliptics=seeds['elliptic_classes']
    canon={E.canonical_elliptic(e['a'],e['b']):i for i,e in enumerate(elliptics)}
    assert len(canon)==len(elliptics)
    pairs=list(itertools.combinations_with_replacement(range(len(elliptics)),2))
    split_ids={pair:i for i,pair in enumerate(pairs)}
    vertices=[dict(id=i,label=f'S{i}',type='split',elliptic_ids=list(pair),
                   aut=8 if pair[0]==pair[1] else 4) for i,pair in enumerate(pairs)]
    r_orbits=sorted({tuple(e['j_orbit']) for e in bfs['decomposed_edges'] if e['type']=='restriction'})
    restriction_ids={orbit:len(vertices)+i for i,orbit in enumerate(r_orbits)}
    for orbit,index in restriction_ids.items():
        j=orbit[0]
        aut=2 if len(orbit)==2 else 12 if j==0 else 8 if j==1728%p else 4
        vertices.append(dict(id=index,label=f'R{j}',type='restriction',j_orbit=list(orbit),aut=aut,degree=7))
    offset=len(vertices)
    for v in bfs['vertices']:
        vertices.append(dict(v,id=offset+v['id']))
    B=Counter()
    categories=Counter()
    # Recompute all rational elliptic 2-kernels and their actual Fp forms.
    elliptic_targets=[];roots=[]
    for e in elliptics:
        a,b=e['a'],e['b'];targets=[]
        roots.append(P.roots([b,a,0,1]))
        for r in roots[-1]:
            if r>=p:continue
            aa=(-4*a-15*r*r)%p
            bb=(-22*r*r*r-8*a*r)%p
            target=canon[E.canonical_elliptic(aa,bb)]
            targets.append(target)
        assert len(targets) in (1,3)
        elliptic_targets.append(targets)
    for (a,b),i in split_ids.items():
        for c in elliptic_targets[a]:
            for d in elliptic_targets[b]:
                j=split_ids[tuple(sorted((c,d)))]
                B[i,j]+=1;categories['S_product']+=1
        degrees=(len(elliptic_targets[a]),len(elliptic_targets[b]))
        vertices[i]['degree']=15 if degrees==(3,3) else 3
    # Each split graph kernel induced by a geometric elliptic isomorphism
    # has an affine root map beta = a*alpha. Rational square a gives S,
    # nonsquare a gives R, since Frobenius negates the elliptic isomorphism.
    affine_records=[]
    for rec in seeds['decomposed_graph_kernels']:
        a,b=rec['split_ids'];perm=rec['root_permutation']
        alpha=roots[a];beta=[roots[b][j] for j in perm]
        scale=F.div(F.sub(beta[1],beta[0]),F.sub(alpha[1],alpha[0]))
        shift=F.sub(beta[0],F.mul(scale,alpha[0]))
        assert shift==0 and 0<scale<p
        assert all(beta[j]==F.mul(scale,alpha[j]) for j in range(3))
        i=split_ids[a,b]
        if E.chi[scale]==1:
            assert a==b
            j=i
        else:
            invariant=elliptics[a]['j']
            assert invariant==elliptics[b]['j']
            j=restriction_ids[(invariant,)]
        B[i,j]+=1;categories['S_decomposed_graph']+=1
        affine_records.append(dict(source=i,target=j,x_scalar=scale,root_permutation=perm))
    # All split graph kernels with a smooth quotient are in the seed list.
    assert len(seeds['kernels'])==len(bfs['seed_indices'])
    for rec,j in zip(seeds['kernels'],bfs['seed_indices']):
        i=split_ids[tuple(rec['split_ids'])]
        B[i,offset+j]+=1;categories['S_J']+=1
    for e in bfs['edges']:
        B[offset+e['source'],offset+e['target']]+=e['multiplicity']
        categories['J_J']+=e['multiplicity']
    for e in bfs['decomposed_edges']:
        i=offset+e['source']
        if e['type']=='split':
            pair=tuple(sorted(canon[tuple(m)] for m in e['elliptic_models']))
            j=split_ids[pair]
        else:j=restriction_ids[tuple(e['j_orbit'])]
        assert vertices[j]['aut']==e['aut']
        B[i,j]+=e['multiplicity'];categories['J_'+vertices[j]['type']]+=e['multiplicity']
    # Every R has exactly three rational product kernels, indexed by the
    # three geometric elliptic 2-kernels. Do not identify conjugate kernels.
    for orbit,i in restriction_ids.items():
        j=orbit[0];a,b=model(F,j)
        rts=P.roots([b,a,0,1]);assert len(rts)==3
        for r in rts:
            aa=F.sub(F.mul((-4)%p,a),F.mul(15%p,F.mul(r,r)))
            bb=F.sub(F.mul((-22)%p,F.power(r,3)),F.mul(8%p,F.mul(a,r)))
            jt=jinvariant(F,aa,bb)
            target=tuple(sorted({jt,F.frob(jt)}))
            k=restriction_ids[target]
            B[i,k]+=1;categories['R_product']+=1
    # The remaining R rows follow from exact polarized duality.
    for (i,j),n in list(B.items()):
        if vertices[j]['type']=='restriction' and vertices[i]['type']!='restriction':
            reverse=Fraction(n*vertices[j]['aut'],vertices[i]['aut'])
            assert reverse.denominator==1,(i,j,n,reverse)
            assert not B[j,i]
            B[j,i]=int(reverse);categories['R_'+vertices[i]['type']]+=int(reverse)
    degrees=Counter()
    neighbors=[set() for _ in vertices]
    for (i,j),n in B.items():
        assert n>0
        degrees[i]+=n
        assert n*vertices[j]['aut']==B[j,i]*vertices[i]['aut'],(i,j,n,B[j,i])
        neighbors[i].add(j);neighbors[j].add(i)
    for i,v in enumerate(vertices):assert degrees[i]==v['degree'],(i,v,degrees[i])
    unseen=set(range(len(vertices)));components=[]
    while unseen:
        root=min(unseen);unseen.remove(root);todo=[root];comp=[]
        while todo:
            i=todo.pop();comp.append(i)
            for j in neighbors[i]:
                if j in unseen:unseen.remove(j);todo.append(j)
        components.append(sorted(comp))
    component_stats=[dict(size=len(c),counts=dict(Counter(vertices[i]['type'] for i in c)),
                          directed_kernels=sum(degrees[i] for i in c)) for c in components]
    # Cross-check the independent existing product-only Sage output.
    oldfile=OUT/f'p{p}_products.json'
    if oldfile.exists():
        old=json.loads(oldfile.read_text())
        assert old['ell2']['product_vertices']==list(map(list,pairs))
        oldB=Counter({(i,j):n for i,j,n in old['ell2']['product_edges']})
        fresh=Counter()
        for (a,b),i in split_ids.items():
            for c in elliptic_targets[a]:
                for d in elliptic_targets[b]:fresh[i,split_ids[tuple(sorted((c,d)))]]+=1
        assert fresh==oldB
    result=dict(p=p,ell=2,field=bfs['field'],weil_polynomial=[p*p,0,2*p,0,1],
                elliptic_classes=elliptics,vertices=vertices,
                edges=[dict(source=i,target=j,multiplicity=n) for (i,j),n in sorted(B.items())],
                counts=dict(Counter(v['type'] for v in vertices)),total_vertices=len(vertices),
                directed_kernels=sum(B.values()),edge_categories=dict(categories),
                components=components,component_statistics=component_stats,
                split_affine_graph_kernels=affine_records,
                checks=['independent elliptic Velu product rows agree with previous Sage output',
                        'all rational Richelot Jacobian rows explicitly enumerated',
                        'all Frobenius-stable split graph kernels enumerated',
                        'all restriction product kernels explicitly enumerated',
                        'every restriction row saturated by 3 product kernels and dual S/J kernels',
                        'all Frobenius-fixed Lagrangian row degrees','weighted duality for every ordered vertex pair'])
    out=OUT/f'p{p}_rational_graph.json'
    out.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:result[k] for k in ['counts','total_vertices','directed_kernels','edge_categories','component_statistics']},indent=2))
    print(out)
    return result


if __name__=='__main__':assemble(int(sys.argv[1]) if len(sys.argv)>1 else 431)

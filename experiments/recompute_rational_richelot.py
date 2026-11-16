"""Recompute small complete rational Richelot rows without SageMath.

Existing exhaustive curve representatives and split rows are retained.
Every Jacobian row is freshly evaluated with the rational correspondence
normalization -[B,C][C,A][A,B]/det(q2,q1,q0). The old kernel records are
updated once, preserving their original branch ordering.
"""
from collections import Counter
import json
from pathlib import Path
from rational_census import Census

OUT=Path(__file__).resolve().parent/'results'
NORMALIZATION='-[B,C][C,A][A,B]/det(rows(q2,q1,q0))'


def run(p):
    path=OUT/f'p{p}_maximal_graph.json'
    data=json.loads(path.read_text())
    E=Census(p);vertices=data['vertices'];oldB=data['adjacency']
    mapping={}
    for old,v in enumerate(vertices):
        if v['type']=='Jacobian':
            new=E.insert(v['coefficients'])
            assert new not in mapping
            mapping[new]=old
    count=len(E.vertices)
    for i in range(count):E.expand(i)
    assert len(E.vertices)==count
    B=[list(row) for row in oldB]
    for i in mapping.values():B[i]=[0]*len(vertices)
    for (i,j),n in E.edges.items():B[mapping[i]][mapping[j]]+=n
    elliptics=data['elliptic_models']
    canonical={E.canonical_elliptic(v['a'],v['b']):i for i,v in enumerate(elliptics)}
    for (i,key),n in E.decomposed.items():
        e=E.decomposed_records[key]
        if e['type']=='restriction':
            assert len(e['j_orbit'])==1
            candidates=[j for j,v in enumerate(vertices) if v['type']=='restriction' and v['j']==e['j_orbit'][0]]
        else:
            pair=sorted(canonical[tuple(model)] for model in e['elliptic_models'])
            candidates=[j for j,v in enumerate(vertices) if v['type']=='split' and v['elliptic_ids']==pair]
        assert len(candidates)==1
        B[mapping[i]][candidates[0]]+=n
    needs_update=data.get('richelot_codomain_normalization')!=NORMALIZATION
    for rec in data['richelot_kernels']:
        if rec['degenerate']:continue
        f=rec['quotient_coefficients']
        if needs_update:f=[-c%p for c in f]
        rec['quotient_coefficients']=f
        new=E.insert(f);assert new<count
        rec['target']=mapping[new]
    records=Counter((r['source'],r['target']) for r in data['richelot_kernels'])
    for i in mapping.values():
        for j in range(len(vertices)):assert records[i,j]==B[i][j]
    for i,row in enumerate(B):
        assert sum(row)==data['outdegrees'][i]
        for j,n in enumerate(row):assert n*vertices[j]['aut']==B[j][i]*vertices[i]['aut']
    unseen=set(range(len(vertices)));components=[]
    while unseen:
        root=min(unseen);unseen.remove(root);pending=[root];component=[]
        while pending:
            i=pending.pop();component.append(i)
            for j in list(unseen):
                if B[i][j] or B[j][i]:unseen.remove(j);pending.append(j)
        components.append(sorted(component))
    changed=sum(B[i][j]!=oldB[i][j] for i in range(len(B)) for j in range(len(B)))
    data['adjacency']=B;data['components']=components
    data['richelot_codomain_normalization']=NORMALIZATION
    check='rational correspondence sign and freshly recomputed Jacobian rows'
    if check not in data['checks']:data['checks'].append(check)
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    print(p,'changed matrix entries',changed,'component sizes',list(map(len,components)),flush=True)
    return dict(p=p,changed_matrix_entries=changed,components=list(map(len,components)))


if __name__=='__main__':
    results=[run(p) for p in (7,11,13,31)]
    (OUT/'rational_richelot_recomputation.json').write_text(json.dumps(results,indent=2)+'\n')

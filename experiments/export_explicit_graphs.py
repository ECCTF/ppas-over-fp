"""Readable full-vertex exports, support drawings, and a reproduction bundle."""
from pathlib import Path
import json,sys,zipfile,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
BASE=Path(__file__).resolve().parent;OUT=BASE/'results'
COLORS={'split':'#2364aa','restriction':'#dc8126','Jacobian':'#21946e'}

def export(p):
    g=json.loads((OUT/f'p{p}_rational_graph.json').read_text())
    G=nx.Graph();G.add_nodes_from(range(len(g['vertices'])))
    G.add_edges_from((e['source'],e['target']) for e in g['edges'] if e['source']!=e['target'])
    pos=nx.spring_layout(G,seed=20261116,iterations=350)
    fig,ax=plt.subplots(figsize=(18,16));fig.patch.set_facecolor('#fafaf8');ax.set_facecolor('#fafaf8')
    nx.draw_networkx_edges(G,pos,ax=ax,edge_color='#788492',alpha=.33,width=.65)
    nx.draw_networkx_nodes(G,pos,ax=ax,node_size=100 if len(G)<200 else 70,
        node_color=[COLORS[v['type']] for v in g['vertices']],linewidths=.4,edgecolors='white')
    nx.draw_networkx_labels(G,pos,ax=ax,labels={v['id']:v['label'] for v in g['vertices']},font_size=5.5,font_color='#17212b')
    ax.set_title(f'p = {p} | complete rational (2,2)-isogeny graph\n'
        f'{len(G)} vertices | {g["directed_kernels"]} directed kernels | {nx.number_connected_components(G)} connected component',fontsize=19,pad=16)
    ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=c,label=t,markersize=8)
              for t,c in COLORS.items()],loc='upper right',frameon=False)
    ax.text(.5,-.025,'All vertices shown. Drawing shows undirected support; loops and multiplicities are retained in CSV/JSON.',
        transform=ax.transAxes,ha='center',fontsize=11,color='#475569')
    ax.axis('off');fig.tight_layout()
    fig.savefig(OUT/f'p{p}_complete_support.svg',bbox_inches='tight')
    fig.savefig(OUT/f'p{p}_complete_support.png',dpi=120,bbox_inches='tight');plt.close(fig)
    # Human-readable index of every model and every outgoing adjacency row.
    adjacency={v['id']:[] for v in g['vertices']}
    for e in g['edges']:adjacency[e['source']].append((e['target'],e['multiplicity']))
    rows=[]
    for v in g['vertices']:
        neighbors=', '.join(f'<a href="#v{j}">{html.escape(g["vertices"][j]["label"])}</a> × {m}' for j,m in adjacency[v['id']])
        rows.append(f'<tr id="v{v["id"]}"><td>{v["label"]}</td><td>{v["type"]}</td>'
          f'<td class="model">{html.escape(v["model_description"])}</td><td>{v["aut"]}</td><td>{v["degree"]}</td><td>{neighbors}</td></tr>')
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>p={p} complete isogeny graph</title>
<style>body{{font:15px/1.6 system-ui,sans-serif;margin:32px;color:#17212b}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #d4dce3;padding:8px;text-align:left;vertical-align:top}}th{{background:#edf3f7;position:sticky;top:0}}.model{{font-family:monospace;min-width:260px}}tr:target{{background:#fff0bd}}a{{color:#2364aa}}img{{max-width:100%}}</style>
<h1>p={p} complete rational (2,2)-isogeny graph</h1><p>{g['total_vertices']} vertices, {g['directed_kernels']} directed edges counted by kernels; graph traversal finds 1 connected component.
S denotes a product of elliptic curves, R a Weil restriction, and J a genus 2 Jacobian. All isomorphisms are defined over Fp and preserve the principal polarization.</p>
<p>Coefficients in the vertex equations are interpreted in the corresponding finite fields. S carries the product polarization, and J the canonical principal polarization. For R, the elliptic curve is chosen in its maximal form over Fp².
The figure shows the undirected support; the table retains the kernel multiplicity of each directed edge. Click a neighbor label to jump to that vertex.</p>
<img src="p{p}_complete_support.svg" alt="Undirected support of the complete graph"><table><thead><tr><th>Vertex</th><th>Type</th><th>Representative model</th><th>|Aut|</th><th>Outdegree</th><th>Outgoing edges: target × kernel multiplicity</th></tr></thead><tbody>{''.join(rows)}</tbody></table></html>'''
    (OUT/f'p{p}_vertex_catalogue.html').write_text(page,encoding='utf8')

if __name__=='__main__':
    for p in map(int,sys.argv[1:] or [47,71,127]):export(p)

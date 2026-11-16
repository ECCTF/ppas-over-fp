"""Publication figures from computed adjacency data; no schematic edges."""
from pathlib import Path
import json, math
from collections import Counter
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';FIG=ROOT.parent/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
blue='#276b98';orange='#d98836';gray='#7f8d98'
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',dpi=180)
    plt.close(fig)

# All 18 Jacobian vertices at p=7; omitted decomposed kernels are annotations.
d=json.loads((OUT/'p7_richelot.json').read_text());vs={v['class_id']:v for v in d['vertices']}
G=nx.Graph();G.add_nodes_from(vs);G.add_edges_from((i,j) for i,j,w in d['edges'])
components=sorted(nx.connected_components(G),key=lambda c:(vs[min(c)]['a_number'],min(c)))
fig,ax=plt.subplots(figsize=(10,5.3));pos={}
for k,c in enumerate(components):
    cx=(k%4)*2.6;cy=-(k//4)*2.3
    for t,v in enumerate(sorted(c)):pos[v]=(cx+(t-(len(c)-1)/2)*1.25,cy)
nx.draw_networkx_nodes(G,pos,ax=ax,node_size=1100,node_color=[blue if vs[v]['a_number']==1 else orange for v in G])
nx.draw_networkx_labels(G,pos,ax=ax,labels={v:str(v) for v in G},font_color='white',font_size=15)
nx.draw_networkx_edges(G,pos,ax=ax,edgelist=[(i,j) for i,j in G.edges if i!=j],width=1.5,edge_color='#314452',node_size=1100)
for v in nx.nodes_with_selfloops(G):
    x,y=pos[v]
    ax.add_patch(FancyArrowPatch((x-.18,y+.15),(x+.18,y+.15),connectionstyle='arc3,rad=-1.8',
                                arrowstyle='-',lw=1.5,color='#314452'))
front=Counter(r['source'] for r in d['geometrically_decomposed_kernel_frontier'])
for v,n in front.items():
    x,y=pos[v];ax.text(x,y-.57,f'+ {n}\ndecomposed',ha='center',va='top',fontsize=13,color='#5d401f')
ax.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=blue,label='a = 1',markersize=9),
                   Line2D([0],[0],marker='o',color='w',markerfacecolor=orange,label='a = 2',markersize=9)],
          loc='upper center',bbox_to_anchor=(.5,1.13),ncol=2,frameon=False,fontsize=15)
ax.margins(.13,.2);ax.axis('off');save(fig,'p7_richelot')

# The actual product-kernel subgraphs on all 903 split vertices.
d=json.loads((OUT/'p431_products.json').read_text());stats={}
fig,axes=plt.subplots(1,2,figsize=(11.4,5.6))
for ell,ax in zip([2,3],axes):
    entry=d[f'ell{ell}'];G=nx.Graph();G.add_nodes_from(range(len(entry['product_vertices'])))
    G.add_edges_from((i,j) for i,j,w in entry['product_edges'] if i!=j)
    comps=sorted(nx.connected_components(G),key=lambda c:(-len(c),min(c)))
    stats[str(ell)]={'vertices':len(G),'directed_kernel_edges':sum(w for i,j,w in entry['product_edges']),
                     'support_components':len(comps),'component_size_histogram':dict(sorted(Counter(map(len,comps)).items())),
                     'loop_kernels':sum(w for i,j,w in entry['product_edges'] if i==j),
                     'degree_histogram':dict(sorted(Counter(sum(w for i,j,w in entry['product_edges'] if i==v) for v in G).items()))}
    # Place each connected component in its own grid cell. Coordinates carry no arithmetic meaning.
    pos={};cols=math.ceil(math.sqrt(len(comps)));rows=math.ceil(len(comps)/cols)
    for k,c in enumerate(comps):
        H=G.subgraph(c);local=nx.spring_layout(H,seed=431,iterations=140) if len(c)>1 else {next(iter(c)):(0,0)}
        for v,(x,y) in local.items():pos[v]=((k%cols)*2.6+x,-(k//cols)*2.6+y)
    nx.draw_networkx_edges(G,pos,ax=ax,width=.48,alpha=.65,edge_color=blue)
    nx.draw_networkx_nodes(G,pos,ax=ax,node_size=4.5,node_color=blue,linewidths=0)
    ax.set_title(f'$({ell},{ell})$ product kernels\n{len(comps)} support components',fontsize=18)
    ax.margins(.08);ax.set_aspect('equal');ax.axis('off')
fig.tight_layout();save(fig,'p431_product_graphs')
(OUT/'graph_statistics.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats,indent=2))

# Full outgoing neighborhood of one Jacobian, with its three computed kernels.
local=json.loads((OUT/'p431_local_richelot.json').read_text())
assert [r['target'] for r in local['kernels']]==['same_jacobian','split_product','split_product']
fig,ax=plt.subplots(figsize=(8,2.8));G=nx.DiGraph()
G.add_edges_from([('J','J'),('J','P26'),('J','P27')]);pos={'J':(0,0),'P26':(-1.8,-.65),'P27':(1.8,-.65)}
nx.draw_networkx_nodes(G,pos,ax=ax,node_color=[orange,blue,blue],node_size=[1400,1300,1300])
nx.draw_networkx_labels(G,pos,ax=ax,labels={'J':'J(C)','P26':r'$E_{26}^{\,2}$','P27':r'$E_{27}^{\,2}$'},font_color='white',font_size=13)
nx.draw_networkx_edges(G,pos,ax=ax,edgelist=[('J','P26'),('J','P27')],node_size=1400,arrows=True,arrowsize=20,width=1.6,edge_color='#314452')
ax.add_patch(FancyArrowPatch((-.22,.13),(.22,.13),connectionstyle='arc3,rad=-1.8',
                            arrowstyle='-|>',mutation_scale=17,lw=1.6,color='#314452'))
ax.text(0,.81,r'$C: y^2=x^5-x$',ha='center',fontsize=13)
for xpos,idx in [(-1.8,26),(1.8,27)]:
    model=d['elliptic_classes'][idx]
    ax.text(xpos,-1.2,rf'$y^2=x^3+{model["a"]}x+{model["b"]}$',ha='center',fontsize=12)
ax.set_xlim(-3.1,3.1);ax.set_ylim(-1.45,1.15);ax.axis('off');save(fig,'p431_local_richelot')

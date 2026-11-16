"""Figures and TeX tables from the complete maximal-block adjacency matrices."""
from pathlib import Path
import json,math
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
BASE=Path(__file__).resolve().parent;OUT=BASE/'results';FIG=BASE.parent/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'font.size':11})
colors={'split':'#246b96','restriction':'#c57a20','Jacobian':'#518a68'}
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');fig.savefig(FIG/(name+'.png'),bbox_inches='tight',dpi=190);plt.close(fig)
def support(d):
    B=d['adjacency'];G=nx.Graph();G.add_nodes_from(range(len(B)))
    G.add_edges_from((i,j) for i in G for j in G if i<j and B[i][j]);return G
def small(d,ax,pos,loop_angles=None):
    loop_angles=loop_angles or {}
    G=support(d);B=d['adjacency'];vs=d['vertices']
    nx.draw_networkx_nodes(G,pos,ax=ax,node_size=950,node_color=[colors[v['type']] for v in vs],edgecolors='white',linewidths=1)
    nx.draw_networkx_labels(G,pos,ax=ax,labels={i:v['label'] for i,v in enumerate(vs)},font_size=12,font_color='white')
    for i,j in G.edges:
        for a,b in [(i,j),(j,i)]:
            nx.draw_networkx_edges(G,pos,ax=ax,edgelist=[(a,b)],arrows=True,arrowstyle='-|>',arrowsize=14,
                                   connectionstyle='arc3,rad=0.15',node_size=1050,width=1.2,edge_color='#566773')
            xa,ya=pos[a];xb,yb=pos[b];dx=xb-xa;dy=yb-ya
            ax.text((xa+xb)/2+dy*.10,(ya+yb)/2-dx*.10,str(B[a][b]),fontsize=11,ha='center',va='center',
                    bbox={'facecolor':'white','edgecolor':'none','pad':.3})
    for i in G:
        if not B[i][i]:continue
        x,y=pos[i];angle=loop_angles.get(i,0)
        def rotate(dx,dy):
            return (x+dx*math.cos(angle)-dy*math.sin(angle),y+dx*math.sin(angle)+dy*math.cos(angle))
        ax.add_patch(FancyArrowPatch(rotate(-.20,.25),rotate(.20,.25),connectionstyle='arc3,rad=-1.8',arrowstyle='-|>',mutation_scale=12,color='#566773'))
        ax.text(*rotate(0,.77),str(B[i][i]),ha='center',va='center',fontsize=11)
    ax.set_aspect('equal');ax.margins(.24,.3);ax.axis('off')
d2=json.loads((OUT/'p7_maximal_graph.json').read_text());d3=json.loads((OUT/'p7_maximal_ell3_graph.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11.5,5.6))
pos2={0:(-1.8,1),1:(1.8,1),2:(0,0),3:(-1.8,-1.1),4:(0,-2.2),5:(1.8,-1.1)}
pos3={0:(-1.3,1.4),4:(1.3,1.4),1:(-1.3,-.5),5:(1.3,-.5),2:(-1.3,-2.4),3:(1.3,-2.4)}
for d,ax,pos in [(d2,axes[0],pos2),(d3,axes[1],pos3)]:
    small(d,ax,pos,{3:math.pi/2,4:math.pi} if d['ell']==2 else None);ax.set_title(f"p = 7, ell = {d['ell']} | {d['directed_kernels']} kernels",fontsize=14,pad=20)
fig.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markersize=9,label=t) for t,c in colors.items()],loc='lower center',ncol=3,frameon=False)
fig.subplots_adjust(bottom=.13,wspace=.22);save(fig,'maximal_p7_full')
for p in [11,13,31]:
    d=json.loads((OUT/f'p{p}_maximal_graph.json').read_text());G=support(d);B=d['adjacency'];vs=d['vertices'];pos={}
    for ci,comp in enumerate(d['components']):
        sub=G.subgraph(comp);local=nx.spring_layout(sub,seed=173,k=1.3/math.sqrt(max(1,len(comp))),iterations=700)
        scale=math.sqrt(len(comp));offset=ci*(math.sqrt(len(G))+2)
        for i,(a,b) in local.items():pos[i]=(a*scale+offset,b*scale)
    fig,ax=plt.subplots(figsize=(11,7 if p==31 else 5))
    nx.draw_networkx_edges(G,pos,ax=ax,width=.9,edge_color='#899aa5',alpha=.8)
    size=420 if p==31 else 760
    for typ,col in colors.items():
        nodes=[i for i,v in enumerate(vs) if v['type']==typ]
        nx.draw_networkx_nodes(G,pos,ax=ax,nodelist=nodes,node_color=col,node_size=size,edgecolors='white',linewidths=.7)
    nx.draw_networkx_labels(G,pos,ax=ax,labels={i:v['label'] for i,v in enumerate(vs)},font_color='white',font_size=10 if p==31 else 13)
    component_word='component' if len(d['components'])==1 else 'components'
    ax.set_title(f"Complete support: p = {p}, ell = 2 | {len(G)} vertices | {len(d['components'])} {component_word}",fontsize=14,pad=15)
    ax.text(.5,-.025,'All vertices and non-loop adjacencies shown; the supplied matrix retains every direction, loop and multiplicity.',transform=ax.transAxes,ha='center',fontsize=9)
    ax.axis('off');ax.margins(.08,.12)
    fig.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=c,markersize=9,label=t) for t,c in colors.items()],loc='lower center',ncol=3,frameon=False)
    fig.subplots_adjust(bottom=.14);save(fig,f'maximal_p{p}_full_support')
    if p==31:
        fig,ax=plt.subplots(figsize=(8.5,7.5));im=ax.imshow(B,cmap='Blues',vmin=0,vmax=max(map(max,B)),interpolation='nearest')
        for b in [20.5,23.5]:ax.axhline(b,color='#c57a20',lw=.8);ax.axvline(b,color='#c57a20',lw=.8)
        ax.set_xticks([10,22,42],['split (21)','Res (3)','Jacobian (36)']);ax.set_yticks([10,22,42],['split (21)','Res (3)','Jacobian (36)'])
        ax.set_xlabel('Target');ax.set_ylabel('Source');ax.set_title('p = 31, ell = 2: complete directed kernel adjacency')
        fig.colorbar(im,ax=ax,label='Number of kernels',shrink=.78);save(fig,'maximal_p31_adjacency')
print('Maximal-block figures written',flush=True)

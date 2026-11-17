"""Plot exact outgoing Richelot neighborhoods from the recorded 10-bit sources."""
from pathlib import Path
import json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

BASE=Path(__file__).resolve().parent
FIG=BASE.parent/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig,grid=plt.subplots(2,2,figsize=(8,7.6))
axes=grid.ravel()
for letter,ax,p in zip('abc',axes,(719,751,823)):
    data=json.loads((BASE/'results'/f'p{p}_jacobian_neighborhood.json').read_text())
    vertices=data['vertices'];n=len(vertices)-1
    xy={0:(0,0)}
    for k in range(1,len(vertices)):
        angle=2*math.pi*(k-1)/n+math.pi/2
        xy[k]=(math.cos(angle),math.sin(angle))
    for edge in data['edges']:
        i,j=edge['source'],edge['target'];x,y=xy[i];u,v=xy[j]
        arrow=FancyArrowPatch((x,y),(u,v),arrowstyle='-|>',mutation_scale=14,
                              shrinkA=19,shrinkB=19,lw=1.7,color='#42566a',zorder=1)
        ax.add_patch(arrow)
        # All labels count distinct source kernels, including multiplicity.
        ax.text(.53*u,.53*v,str(edge['forward_kernels']),ha='center',va='center',
                fontsize=14,bbox={'boxstyle':'round,pad=.16','fc':'white','ec':'none'},zorder=4)
    for k,vertex in enumerate(vertices):
        x,y=xy[k];color='#f4bc73' if vertex['type']=='split' else '#b7d9ef'
        ax.scatter([x],[y],s=760,color=color,edgecolors='#253e53',linewidths=1.4,zorder=2)
        ax.text(x,y,vertex['label'],ha='center',va='center',fontsize=14,fontweight='bold',zorder=3)
    ax.text(.5,1.02,r'(%s) $p=%d,\ a=%d$'%(letter,p,data['source_parameter_a']),ha='center',transform=ax.transAxes,fontsize=14)
    ax.set_xlim(-1.4,1.4);ax.set_ylim(-1.4,1.4);ax.set_aspect('equal');ax.axis('off')
axes[3].axis('off')
axes[3].text(.5,.64,'J: Jacobian\n\nS: split product',ha='center',va='center',fontsize=14,transform=axes[3].transAxes)
fig.subplots_adjust(left=.04,right=.96,top=.91,bottom=.04,wspace=.18,hspace=.30)
for ext in ('pdf','png'):fig.savefig(FIG/f'ten_bit_jacobian_neighborhoods.{ext}',dpi=200,bbox_inches='tight')
print('Saved ten_bit_jacobian_neighborhoods.pdf/png')

"""Reproducible vector drawings of the complete restriction sectors.

Uses ReportLab only. Edges in the support plots are undirected and unweighted.
The accompanying matrix (and all CSV files) retain every directed kernel.
"""
from pathlib import Path
import json, math, random
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, Color, black, white

BASE=Path(__file__).resolve().parent
FIG=BASE.parent/'figures'
BLUE=HexColor('#b7d7ec');GOLD=HexColor('#f2d48a');GREY=HexColor('#89949b')

def load(p):return json.loads((BASE/'results'/f'p{p}_ten_bit_restriction_ell3.json').read_text())

def layout(B):
    n=len(B);rng=random.Random(20261105)
    pos=[[rng.uniform(-1,1),rng.uniform(-1,1)] for _ in range(n)]
    edges=[(i,j) for i in range(n) for j in range(i+1,n) if B[i][j] or B[j][i]]
    k=math.sqrt(4/n)
    for it in range(900):
        disp=[[0.,0.] for _ in pos]
        for i in range(n):
            for j in range(i+1,n):
                dx=pos[i][0]-pos[j][0];dy=pos[i][1]-pos[j][1];d=max(math.hypot(dx,dy),.0001)
                z=k*k/d/d
                disp[i][0]+=dx*z;disp[i][1]+=dy*z;disp[j][0]-=dx*z;disp[j][1]-=dy*z
        for i,j in edges:
            dx=pos[i][0]-pos[j][0];dy=pos[i][1]-pos[j][1];d=max(math.hypot(dx,dy),.0001)
            z=d/k
            disp[i][0]-=dx*z;disp[i][1]-=dy*z;disp[j][0]+=dx*z;disp[j][1]+=dy*z
        temp=.10*(1-it/900)+.001
        for i in range(n):
            d=max(math.hypot(*disp[i]),.0001)
            pos[i][0]+=disp[i][0]*min(d,temp)/d
            pos[i][1]+=disp[i][1]*min(d,temp)/d
    lo=[min(z[c] for z in pos) for c in (0,1)];hi=[max(z[c] for z in pos) for c in (0,1)]
    return [[(z[c]-lo[c])/(hi[c]-lo[c]) for c in (0,1)] for z in pos]

def support(c,d,x,y,w,h,r=8):
    B=d['adjacency'];pos=layout(B)
    pos=[(x+18+a*(w-36),y+18+b*(h-36)) for a,b in pos]
    c.setStrokeColor(GREY);c.setLineWidth(.7)
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            if B[i][j] or B[j][i]:c.line(*pos[i],*pos[j])
    for i,(xx,yy) in enumerate(pos):
        if B[i][i]:
            c.setStrokeColor(HexColor('#343e44'));c.setLineWidth(.8)
            c.circle(xx,yy+r+4,5,stroke=1,fill=0)
    for v,(xx,yy) in zip(d['vertices'],pos):
        c.setFillColor(BLUE if len(v['j_orbit'])==1 else GOLD)
        c.setStrokeColor(HexColor('#33434d'));c.setLineWidth(.65)
        c.circle(xx,yy,r,stroke=1,fill=1)
        c.setFillColor(black);c.setFont('Helvetica',max(6,r-.2));c.drawCentredString(xx,yy-2.5,str(v['id']))
    c.setFillColor(black);c.setFont('Helvetica-Bold',12)
    c.drawCentredString(x+w/2,y+h+16,f"p = {d['p']}: {len(B)} vertices, {d['directed_kernels']} kernels")

def legend(c,x,y):
    for col,label,dx in [(BLUE,'j in F_p',0),(GOLD,'conjugate nonrational pair',100)]:
        c.setFillColor(col);c.circle(x+dx,y,5,stroke=1,fill=1)
        c.setFillColor(black);c.setFont('Helvetica',9);c.drawString(x+dx+9,y-3,label)

def main():
    FIG.mkdir(exist_ok=True)
    d=load(991);out=FIG/'ten_bit_restriction_p991.pdf'
    c=canvas.Canvas(str(out),pagesize=(950,540));c.setTitle('Complete maximal Weil-restriction (3,3) sector at p=991')
    support(c,d,15,80,505,410,8)
    B=d['adjacency'];cell=7.4;xx=550;yy=95;n=len(B)
    for i in range(n):
        for j in range(n):
            val=B[i][j]
            c.setFillColor(white if not val else Color(.08,.34,.58,alpha=1))
            c.rect(xx+j*cell,yy+(n-1-i)*cell,cell,cell,stroke=0,fill=1)
            if val:c.setFillColor(white);c.setFont('Helvetica',5.5);c.drawCentredString(xx+(j+.5)*cell,yy+(n-1-i)*cell+1.8,str(val))
    c.setStrokeColor(black);c.setLineWidth(.4);c.rect(xx,yy,n*cell,n*cell,stroke=1,fill=0)
    c.setFillColor(black);c.setFont('Helvetica-Bold',12);c.drawCentredString(xx+n*cell/2,506,'Directed kernel adjacency matrix')
    c.setFont('Helvetica',8)
    for i in range(0,n,5):
        c.drawRightString(xx-4,yy+(n-1-i)*cell+1,str(i));c.drawCentredString(xx+(i+.5)*cell,yy-11,str(i))
    c.drawCentredString(xx+n*cell/2,yy-27,'codomain vertex');c.saveState();c.translate(xx-26,yy+n*cell/2);c.rotate(90);c.drawCentredString(0,0,'source vertex');c.restoreState()
    legend(c,95,43);c.setFont('Helvetica',9);c.drawCentredString(475,20,'Left: undirected support, with loops. Right: all directed multiplicities. One connected component.')
    c.showPage();c.save()
    standalone=FIG/'ten_bit_restriction_p991_support.pdf'
    c=canvas.Canvas(str(standalone),pagesize=(560,550))
    c.setTitle('Complete restriction (3,3) graph support at p=991')
    support(c,d,10,75,540,420,9)
    legend(c,125,40);c.setFont('Helvetica',9)
    c.drawCentredString(280,18,'Complete undirected support, including loops. Edge multiplicities are omitted.')
    c.showPage();c.save()
    out2=FIG/'ten_bit_restriction_comparison.pdf';c=canvas.Canvas(str(out2),pagesize=(1120,440))
    c.setTitle('Complete maximal restriction sectors for three 10-bit primes')
    for k,p in enumerate([523,673,991]):support(c,load(p),12+k*370,80,356,310,7)
    legend(c,370,42);c.setFont('Helvetica',9);c.drawCentredString(560,18,'Complete undirected supports of the closed (3,3) restriction sectors. Each vertex has four outgoing kernels.')
    c.showPage();c.save()
    print(out);print(standalone);print(out2)

if __name__=='__main__':main()

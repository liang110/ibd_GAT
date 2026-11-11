"""Figure 1 — auditable directional-transfer framework."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manuscript_rebuild" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"svg.fonttype":"none","pdf.fonttype":42,"ps.fonttype":42,"font.size":7.0})
INK="#1F2933"; MUTED="#5B6770"; LINE="#6B7280"; NAVY="#28577D"; TEAL="#3FA7A3"; ORANGE="#D77A32"; PURPLE="#7A5A9E"
PALE_BLUE="#EAF2F8"; PALE_TEAL="#EAF7F5"; PALE_ORANGE="#FBEFE5"; PALE_PURPLE="#F2EDF7"; PALE_GRAY="#F5F7F8"

def card(ax,x,y,w,h,title,lines,edge,fill,title_size=7.2,body_size=6.4,radius=.018):
    p=FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.006,rounding_size={radius}",facecolor=fill,edgecolor=edge,linewidth=1.0,zorder=2); ax.add_patch(p)
    ax.text(x+.035*w,y+h-.23*h,title,ha="left",va="center",color=edge,fontweight="bold",fontsize=title_size,zorder=3)
    ax.text(x+.035*w,y+h-.49*h,"\n".join(lines),ha="left",va="top",color=INK,fontsize=body_size,linespacing=1.25,zorder=3)
    return p

def plain_box(ax,x,y,w,h,text,edge,fill,fs=6.5,weight="normal"):
    p=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.006,rounding_size=.014",facecolor=fill,edgecolor=edge,linewidth=1.0,zorder=2); ax.add_patch(p)
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",color=INK,fontsize=fs,fontweight=weight,linespacing=1.25,zorder=3); return p

def arrow(ax,a,b,color=LINE,lw=1.0,ls="-",rad=0.0):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=8,linewidth=lw,linestyle=ls,color=color,connectionstyle=f"arc3,rad={rad}",zorder=1))

def graph_glyph(ax,cx,cy,scale=.8):
    pts=[(cx-.07*scale,cy+.045*scale),(cx-.005*scale,cy+.085*scale),(cx+.07*scale,cy+.045*scale),(cx+.045*scale,cy-.065*scale),(cx-.05*scale,cy-.065*scale)]
    for i,j in [(0,1),(1,2),(2,3),(3,4),(4,0),(1,4),(1,3)]: ax.plot([pts[i][0],pts[j][0]],[pts[i][1],pts[j][1]],color=PURPLE,lw=.55,alpha=.75,zorder=4)
    for x,y in pts: ax.add_patch(Circle((x,y),.011*scale,facecolor="white",edgecolor=PURPLE,lw=.7,zorder=5))
    ax.plot([pts[0][0],pts[1][0]],[pts[0][1],pts[1][1]],color=ORANGE,lw=1.3,zorder=5); ax.plot([pts[1][0],pts[3][0]],[pts[1][1],pts[3][1]],color=TEAL,lw=1.3,zorder=5)

def main():
    fig,ax=plt.subplots(figsize=(5.8,5.15),facecolor="white"); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    ax.text(.035,.970,"Source-only, confound-audited directional transfer",fontsize=10.4,fontweight="bold",color=INK,ha="left",va="top")
    ax.text(.035,.935,"A fixed microbiome contract separates representation effects from cohort, participant and topology artifacts.",fontsize=6.5,color=MUTED,ha="left",va="top")
    ax.text(.035,.882,"a",fontsize=8.2,fontweight="bold",color=INK); ax.text(.058,.882,"Input contract",fontsize=7.8,fontweight="bold",color=INK)
    card(ax,.035,.735,.205,.105,"HMP2 / IBDMDB",["1,638 samples · 578 features"],NAVY,PALE_BLUE)
    card(ax,.035,.595,.205,.105,"PRISM",["220 samples · 209 features"],TEAL,PALE_TEAL)
    plain_box(ax,.065,.450,.145,.085,"126 shared microbiome\nfeatures\nmetabolite overlap = 0",ORANGE,PALE_ORANGE,6.1,"bold")
    # Route HMP2 around the PRISM card; the two inputs must converge on the
    # shared contract without an arrow visually passing through the cohort box.
    ax.plot([.205,.252,.252,.180],[.735,.705,.545,.535],color=NAVY,lw=.9,zorder=1)
    ax.add_patch(FancyArrowPatch((.252,.545),(.180,.535),arrowstyle="-|>",mutation_scale=8,color=NAVY,lw=.9,zorder=1))
    arrow(ax,(.137,.595),(.137,.536),TEAL,lw=.9)
    ax.text(.137,.425,"microbiome-only primary transfer",fontsize=5.9,color=ORANGE,ha="center",va="top")
    ax.text(.270,.882,"b",fontsize=8.2,fontweight="bold",color=INK); ax.text(.293,.882,"Source-only fitting",fontsize=7.8,fontweight="bold",color=INK)
    card(ax,.270,.735,.205,.105,"Harmonize",["canonical IDs · closure → CLR","denylist clinical / ID fields"],NAVY,PALE_BLUE,body_size=6.1)
    card(ax,.270,.595,.205,.105,"Split and fit",["participant-grouped source splits","prevalence · scale · MI fit in source"],TEAL,PALE_TEAL,body_size=6.1)
    card(ax,.270,.455,.205,.105,"Freeze",["feature order + hash recorded","target remains locked"],ORANGE,PALE_ORANGE,body_size=6.1)
    arrow(ax,(.240,.492),(.270,.492),ORANGE,lw=1.0); arrow(ax,(.3725,.735),(.3725,.700),NAVY,lw=.85); arrow(ax,(.3725,.595),(.3725,.560),TEAL,lw=.85)
    ax.text(.505,.882,"c",fontsize=8.2,fontweight="bold",color=INK); ax.text(.528,.882,"Representation branches",fontsize=7.8,fontweight="bold",color=INK)
    card(ax,.505,.700,.215,.125,"Tabular references",["L2 logistic · elastic net","random forest"],ORANGE,PALE_GRAY,body_size=6.2)
    for i,h in enumerate([.018,.030,.023,.037]): ax.add_patch(Rectangle((.666+i*.011,.724),.007,h,facecolor=ORANGE,alpha=.75,edgecolor="none",zorder=4))
    card(ax,.505,.535,.215,.125,"No-Prior graph branch",["Spearman top-k graph","graph attention → mean pool"],PURPLE,PALE_PURPLE,body_size=6.0); graph_glyph(ax,.695,.575,.58)
    arrow(ax,(.475,.482),(.505,.762),ORANGE,lw=.9,rad=.12); arrow(ax,(.475,.482),(.505,.597),PURPLE,lw=.9,rad=-.12)
    plain_box(ax,.755,.585,.205,.105,"Controls\nno-edge · identity\nrewired topology · MLP",PURPLE,PALE_GRAY,6.0); arrow(ax,(.720,.597),(.755,.637),PURPLE,lw=.8,ls="--",rad=.10)
    ax.text(.035,.378,"d",fontsize=8.2,fontweight="bold",color=INK); ax.text(.058,.378,"Directional transfer and evaluation",fontsize=7.8,fontweight="bold",color=INK)
    plain_box(ax,.035,.235,.185,.095,"source\nHMP2 / IBDMDB",NAVY,PALE_BLUE,6.3,"bold"); plain_box(ax,.270,.235,.205,.095,"target locked after fitting\nPRISM",TEAL,PALE_TEAL,6.3,"bold")
    arrow(ax,(.220,.300),(.270,.300),NAVY,lw=1.0); ax.text(.245,.360,"HMP2 → PRISM",fontsize=5.8,color=NAVY,fontweight="bold",ha="center")
    arrow(ax,(.270,.270),(.220,.270),TEAL,lw=.9); ax.text(.245,.215,"PRISM → HMP2",fontsize=5.6,color=TEAL,fontweight="bold",ha="center",va="top")
    plain_box(ax,.035,.105,.185,.075,"PRISM discovery\n155 samples",TEAL,PALE_TEAL,5.9); plain_box(ax,.270,.105,.205,.075,"PRISM validation\n65 samples",NAVY,PALE_BLUE,5.9)
    arrow(ax,(.220,.142),(.270,.142),TEAL,lw=.9); ax.text(.245,.168,"origin-held-out",fontsize=5.7,color=TEAL,ha="center")
    ax.text(.505,.378,"e",fontsize=8.2,fontweight="bold",color=INK); ax.text(.528,.378,"Evidence gates",fontsize=7.8,fontweight="bold",color=INK)
    card(ax,.505,.258,.205,.082,"Validity",["leakage + known-truth controls"],NAVY,PALE_BLUE,body_size=5.9); card(ax,.735,.258,.225,.082,"Uncertainty",["paired participant-cluster bootstrap"],ORANGE,PALE_ORANGE,body_size=5.9)
    card(ax,.505,.135,.205,.082,"Topology",["edge Jaccard + degree diagnostic"],PURPLE,PALE_PURPLE,body_size=5.9); card(ax,.735,.135,.225,.082,"Scope",["candidate features; no causal claim"],TEAL,PALE_TEAL,body_size=5.9)
    arrow(ax,(.475,.282),(.505,.299),ORANGE,lw=.8,ls="--"); arrow(ax,(.475,.282),(.505,.176),PURPLE,lw=.8,ls="--")
    ax.text(.035,.042,"Solid arrows: data/model flow",fontsize=5.6,color=MUTED,ha="left"); ax.text(.340,.042,"Dashed arrows: audit controls",fontsize=5.6,color=MUTED,ha="left")
    fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
    for ext,dpi in [("png",600),("tiff",600),("pdf",None),("svg",None)]: fig.savefig(OUT/f"Figure1_framework.{ext}",bbox_inches="tight",facecolor="white",**({"dpi":dpi} if dpi else {}))
    plt.close(fig); print(OUT/"Figure1_framework.svg")

if __name__ == "__main__": main()

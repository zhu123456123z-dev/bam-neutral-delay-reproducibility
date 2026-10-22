from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parents[1]/'outputs'/'python'
out.mkdir(parents=True,exist_ok=True)
x=np.linspace(.05,5,700);y=np.linspace(.05,12,700);X,Y=np.meshgrid(x,y)
region=np.zeros_like(X);region[.105125*X*Y<1]=1;region[.92*X*Y<1]=2;region[(2.3*X<1)&(.4*Y<1)]=3
from matplotlib.colors import ListedColormap
fig,ax=plt.subplots(figsize=(8.8,4.8));ax.contourf(X,Y,region,levels=[-.5,.5,1.5,2.5,3.5],cmap=ListedColormap(['white','#fde5de','#dceaf8','#e1e1e1']))
ax.plot(x,1/(.105125*x),color='#b94227',label='Separated spectral boundary')
ax.plot(x,1/(.92*x),color='#246491',label='Merged spectral boundary')
ax.plot([0,1/2.3,1/2.3],[2.5,2.5,0],color='#444444',label='Unweighted row boundary')
ax.scatter([1,1],[1,2],c='black',s=28,zorder=5);ax.annotate('case 1',(1,1),xytext=(1.1,.6));ax.annotate('case 2',(1,2),xytext=(1.1,2.2))
ax.set(xlim=(0,5),ylim=(0,12),xlabel=r'$\xi$',ylabel=r'$\eta$');ax.legend(loc='upper right');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'fig4_regions.png',dpi=240)

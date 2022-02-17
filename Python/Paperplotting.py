# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.colorbar
path = "C:/Users/KingofPi/Desktop/Masterthemen/Kont.-limes/Kuramoto-Sakaguchi-master/Python/Data"
def Figure2():
    annotations = ["a","b","c","d","e","f","g","h"]
    data = [(18,0.9),(18,0.95),(18,1.05),(18,1.75),(29,0.55),(29,0.9),(29,1.75),(29,2.75)]
    fig = plt.figure()
    grid = ImageGrid(fig,111,
    				 nrows_ncols=(2,4),
    				 axes_pad=0.,
    				 share_all=True,
    				 cbar_location="right",
    				 cbar_mode="single",
    				 cbar_size="7%",
    				 cbar_pad=0.15
    				 )
    for ax,data,anno in zip(grid,data,annotations):
        im = ax.imshow(np.loadtxt(path+f"/{data[0]}_{data[1]}_coupling.txt"),vmin=-1,vmax=1,interpolation="nearest",cmap=plt.get_cmap("bwr"),origin="lower",extent=(1,50,1,50))
        ax.set_xlabel("Index j")
        ax.set_ylabel("Index i")
        ax.annotate(anno,xy=(2,49),va="top",ha="left",backgroundcolor="white",bbox=dict(boxstyle="square,pad=0.0",fc="white",ec="white"),fontsize=8)
        #ax.annotate(1,50,anno,va="top",ha="left",backgroundcolor="white")
    
    cbar=matplotlib.colorbar.Colorbar(ax.cax,im)
    ax.cax.toggle_label(True)
    cbar.ax.set_ylabel("$\kappa_{ij}$")
    grid[0].annotate("$\sigma$", xy=(4.5, 1.1), xycoords='axes fraction', xytext=(0, 1.1), 
                arrowprops=dict(arrowstyle="->", color="black"),fontsize=10)
    plt.savefig("Figure2.png")
    plt.savefig("Figure2.svg")
    plt.savefig("Figure2.pdf")
    
    
Figure2()
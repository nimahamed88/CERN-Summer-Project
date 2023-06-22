import matplotlib.pyplot as plt 
import mplhep as hep 
hep.style.use("CMS") 
import numpy as np
import pandas as pd

def inspectDistribution(dflist,
                        var,xlabel,bins,ylabel='Events/bin',xscale='linear',yscale='linear'):
    """
    make a simple histogram out of a variable
    dflist is a list of (dataframe, mask, label) to overlay in the plot
    if mask is None, a trivial mask (all events passing) is used
    if documul=True, the cumulative distribution is added
    """
    
    fig,ax=plt.subplots(figsize=(8,8))
    for df,mask,label in dflist:
        if mask is None: mask=np.ones(df.shape[0],dtype=bool)
        ax.hist(df[mask][var],label=label,histtype='step',lw=2,bins=bins)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    ax.grid()
    if len(dflist)>3:
        ax.legend(bbox_to_anchor=(1.,1.),fontsize=16)
        fig.set_figwidth(12)
    else:
        ax.legend()
    hep.cms.label(loc=1,data=True,rlabel=r'$\sqrt{s_{NN}}=13.6 TeV$')
    fig.tight_layout()
    plt.show()

# import libraries
import pandas as pd
import os
import sys
import math
import numpy as np    
import scanpy as sc
import matplotlib.pyplot as plt
#from test_functions import *
import matplotlib.backends.backend_pdf


def PlotActivity(adata=None,prop=None, outpath="./", fraction=0.05, ncol=5, minnx=None, maxxx=None, minny=None, maxxy=None, spot_size=100, random_state=0, sample_name='sample', label='HDP', ntopics=0):# n_neighbors=0, n_components=0, alpha_W=0):
    if outpath == './':
        outpath = str(outpath) + "/Plot_activity_sub"+str(fraction)+".pdf"

    prop_W=np.array(prop)
    # for plotting each HDP
    for ii in range(prop_W.shape[1]):
        adata.obs[prop.columns[ii]]=prop_W[:,ii]
    adata.obsm['spatial']=np.array(adata.obsm['spatial'])
    adata_sub_sub=sc.pp.subsample(adata, fraction=fraction, copy=True, random_state=random_state)

    if (minnx is None or maxxx is None or minny is None or maxxy is None):
        minnx=min(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
        maxxx=max(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
        minny=min(np.ndarray.flatten(adata.obsm['spatial'][:,1]))
        maxxy=max(np.ndarray.flatten(adata.obsm['spatial'][:,1]))

    #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
    pdf = matplotlib.backends.backend_pdf.PdfPages(outpath)
    nrow = math.ceil(ntopics/ncol)
    fig, ax = plt.subplots(nrow,ncol, figsize=(6*ncol+2,6*nrow))
    for ii in range(nrow):
        for jj in range(ncol):
            if (ii*ncol+jj)<ntopics:
                ax[ii,jj].set_xlim(minnx,maxxx)
                ax[ii,jj].set_ylim(minny,maxxy)
    for ii in range(nrow):
        for jj in range(ncol):
            if (ii*ncol+jj)<ntopics:
                sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[ii*ncol+jj]], ax=ax[ii,jj], show=False)

    plt.tight_layout(pad=3.0)
    pdf.savefig( fig )
    pdf.close()



def PlotMajCluster(adata=None,majcluster='0', outpath="./", fraction=1, minnx=None, maxxx=None, minny=None, maxxy=None, spot_size=100, sample_name='sample', label='HDP', ntopics=0):# n_neighbors=0, n_components=0, alpha_W=0):
    if outpath == './':
        outpath = str(outpath) + "/Plot_major_cluster_sub"+str(fraction)+".pdf"

    # for plotting
    adata.obs['HDP_cluster']=list(map(str, majcluster))
    adata.obsm['spatial']=np.array(adata.obsm['spatial'])
    if (minnx is None or maxxx is None or minny is None or maxxy is None):
        minnx=min(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
        maxxx=max(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
        minny=min(np.ndarray.flatten(adata.obsm['spatial'][:,1]))
        maxxy=max(np.ndarray.flatten(adata.obsm['spatial'][:,1]))

    #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
    pdf = matplotlib.backends.backend_pdf.PdfPages(outpath)
    fig, (ax1) = plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'wspace':0.9})
    ax1.set_xlim(minnx,maxxx)
    ax1.set_ylim(minny,maxxy)
    sc.pl.spatial(adata, spot_size=spot_size, color=['HDP_cluster'],title="Major Cell programs" , show=False, ax=ax1)
    pdf.savefig( fig )
    pdf.close()


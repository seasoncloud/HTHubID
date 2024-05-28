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


def PlotActivity(adata=None,prop=None, outpath="./", fraction=0.05, ncol=5, minnx=None, maxxx=None, minny=None, maxxy=None, spot_size=100, random_state=0, sample_name='sample', label='HDP', ntopics=0, multi_samples=False, plot_set=None):# n_neighbors=0, n_components=0, alpha_W=0):

    prop_W=np.array(prop)
    # for plotting each HDP
    for ii in range(prop_W.shape[1]):
        adata.obs[prop.columns[ii]]=prop_W[:,ii]
    adata.obsm['spatial']=np.array(adata.obsm['spatial'])

    # # # select prop
    # prop2=np.array(prop)
    # sel_idx=[]
    # for ii in range((prop2.shape[1])):
    #     rr=np.ptp(prop2[:,ii])
    #     if rr>0.000000001:
    #         sel_idx.append(ii)

    # prop=prop        

    if multi_samples== False:
        if prop.shape[1]<=30:
            if fraction is not None:
                adata_sub_sub=sc.pp.subsample(adata, fraction=fraction, copy=True, random_state=random_state)
            if (minnx is None or maxxx is None or minny is None or maxxy is None):
                minnx=min(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
                maxxx=max(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
                minny=min(np.ndarray.flatten(adata.obsm['spatial'][:,1]))
                maxxy=max(np.ndarray.flatten(adata.obsm['spatial'][:,1]))

            #if outpath == './':
            outpath = str(outpath) + "Plot_activity_sub"+str(fraction)+".pdf"

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
                        if fraction is None:
                            adata_sub_sub= adata[prop.iloc[:,(30*pp+(ii*ncol+jj))]>0,:]#sc.pp.subsample(adata, fraction=fraction, copy=True, random_state=random_state)
                        sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[ii*ncol+jj]], ax=ax[ii,jj], show=False)

            plt.tight_layout(pad=3.0)
            pdf.savefig( fig )
            pdf.close()
        else:
            if fraction is not None:
                adata_sub_sub=sc.pp.subsample(adata, fraction=fraction, copy=True, random_state=random_state)
            
            if (minnx is None or maxxx is None or minny is None or maxxy is None):
                minnx=min(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
                maxxx=max(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
                minny=min(np.ndarray.flatten(adata.obsm['spatial'][:,1]))
                maxxy=max(np.ndarray.flatten(adata.obsm['spatial'][:,1]))

            #if outpath == './':
            outpath = str(outpath) + "Plot_activity_sub"+str(fraction)+".pdf"

            #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
            pdf = matplotlib.backends.backend_pdf.PdfPages(outpath)
            nrow_all = math.ceil(ntopics/ncol)
            nrow=6

            for pp in range(math.ceil(nrow_all/nrow)):
                fig, ax = plt.subplots(nrow,ncol, figsize=(6*ncol+2,6*nrow))
                for ii in range(nrow):
                    for jj in range(ncol):
                        if (30*pp+(ii*ncol+jj))<ntopics:
                            ax[ii,jj].set_xlim(minnx,maxxx)
                            ax[ii,jj].set_ylim(minny,maxxy)
                for ii in range(nrow):
                    for jj in range(ncol):
                        if (30*pp+(ii*ncol+jj))<ntopics:
                            if fraction is None:
                                adata_sub_sub= adata[prop.iloc[:,(30*pp+(ii*ncol+jj))]>0,:]#sc.pp.subsample(adata, fraction=fraction, copy=True, random_state=random_state)
                            sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[(30*pp+(ii*ncol+jj))]], ax=ax[ii,jj], show=False)
                plt.tight_layout(pad=3.0)
                pdf.savefig( fig )
            pdf.close()

    else:
        if prop.shape[1]<=30:
            multi_sample_names=[item.split("_")[0] for item in adata.obs['index']]
            adata.obs['multi_sample_names']=multi_sample_names

            if plot_set is None:
                plot_set=list(set(multi_sample_names))
            else:
                print("plot_set should be a list")
            for ss in plot_set:
                adata_sub_sub=adata[adata.obs['multi_sample_names']==ss]
                if fraction is not None:
                    adata_sub_sub=sc.pp.subsample(adata_sub_sub, fraction=fraction, copy=True, random_state=random_state)

                if (minnx is None or maxxx is None or minny is None or maxxy is None):
                    minnx=min(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,0]))
                    maxxx=max(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,0]))
                    minny=min(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,1]))
                    maxxy=max(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,1]))

                #if outpath == './':
                outpath2 = str(outpath) + "Plot_activity_sample_"+str(ss)+"_sub"+str(fraction)+".pdf"

                #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
                pdf = matplotlib.backends.backend_pdf.PdfPages(outpath2)
                nrow = math.ceil(ntopics/ncol)
                fig, ax = plt.subplots(nrow,ncol, figsize=(6*ncol+2,6*nrow*2))
                if nrow>1:
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (ii*ncol+jj)<ntopics:
                                ax[ii,jj].set_xlim(minnx,maxxx)
                                ax[ii,jj].set_ylim(minny,maxxy)
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (ii*ncol+jj)<ntopics:
                                sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[ii*ncol+jj]], ax=ax[ii,jj], show=False)
                else:
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (ii*ncol+jj)<ntopics:
                                ax[jj].set_xlim(minnx,maxxx)
                                ax[jj].set_ylim(minny,maxxy)
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (ii*ncol+jj)<ntopics:
                                sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[ii*ncol+jj]], ax=ax[jj], show=False)
                plt.tight_layout(pad=3.0)
                pdf.savefig( fig )
                pdf.close()

        else:
            multi_sample_names=[item.split("_")[0] for item in adata.obs['index']]
            adata.obs['multi_sample_names']=multi_sample_names

            if plot_set is None:
                plot_set=list(set(multi_sample_names))
            else:
                print("plot_set should be a list")
            for ss in plot_set:
                adata_sub_sub=adata[adata.obs['multi_sample_names']==ss]
                if fraction is not None:
                    adata_sub_sub=sc.pp.subsample(adata_sub_sub, fraction=fraction, copy=True, random_state=random_state)

                if (minnx is None or maxxx is None or minny is None or maxxy is None):
                    minnx=min(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,0]))
                    maxxx=max(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,0]))
                    minny=min(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,1]))
                    maxxy=max(np.ndarray.flatten(adata_sub_sub.obsm['spatial'][:,1]))

                #if outpath == './':
                outpath2 = str(outpath) + "Plot_activity_sample_"+str(ss)+"_sub"+str(fraction)+".pdf"

                #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
                pdf = matplotlib.backends.backend_pdf.PdfPages(outpath2)
                nrow_all = math.ceil(ntopics/ncol)
                nrow=6
                for pp in range(math.ceil(nrow_all/nrow)):
                    fig, ax = plt.subplots(nrow,ncol, figsize=(6*ncol+2,6*nrow*2))
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (30*pp+(ii*ncol+jj))<ntopics:
                                ax[jj].set_xlim(minnx,maxxx)
                                ax[jj].set_ylim(minny,maxxy)
                    for ii in range(nrow):
                        for jj in range(ncol):
                            if (30*pp+(ii*ncol+jj))<ntopics:
                                sc.pl.spatial(adata_sub_sub, spot_size=spot_size, color=[prop.columns[30*pp+(ii*ncol+jj)]], ax=ax[jj], show=False)
                    plt.tight_layout(pad=3.0)
                    pdf.savefig( fig )
                pdf.close()





def PlotMajCluster(adata=None, majcluster='0', outpath="./", fraction=1, minnx=None, maxxx=None, minny=None, maxxy=None, spot_size=100, sample_name='sample', label='HDP', ntopics=0, multi_samples=False, multi_sample_names=None, plot_set=None, palette=None):# n_neighbors=0, n_components=0, alpha_W=0):

    adata.obs['HDP_cluster']=list(map(str, majcluster))
    adata.obsm['spatial']=np.array(adata.obsm['spatial'])

    if multi_samples== False:
        # for plotting
        if (minnx is None or maxxx is None or minny is None or maxxy is None):
            minnx=min(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
            maxxx=max(np.ndarray.flatten(adata.obsm['spatial'][:,0]))
            minny=min(np.ndarray.flatten(adata.obsm['spatial'][:,1]))
            maxxy=max(np.ndarray.flatten(adata.obsm['spatial'][:,1]))

        #if outpath == './':
        outpath = str(outpath) + "Plot_major_cluster_sub"+str(fraction)+".pdf"
        #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
        pdf = matplotlib.backends.backend_pdf.PdfPages(outpath)
        fig, (ax1) = plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'wspace':0.9})
        ax1.set_xlim(minnx,maxxx)
        ax1.set_ylim(minny,maxxy)
        sc.pl.spatial(adata, spot_size=spot_size, color=['HDP_cluster'],palette=palette,title="Major Cell programs" , show=False, ax=ax1)
        pdf.savefig( fig )
        pdf.close()

    else:
        #if multi_sample_names is None:
        multi_sample_names=[item.split("_")[0] for item in adata.obs['index']]
        # else:
        #     print("multi_sample_names sould be a list.")
        adata.obs['multi_sample_names']=multi_sample_names
        print(adata.obs)
        print(multi_sample_names[0:10])

        if plot_set is None:
            plot_set=list(set(multi_sample_names))
        else:
            print("plot_set should be a list")
        for ss in plot_set:
            adata_sub=adata[adata.obs['multi_sample_names']==ss]
            # for plotting
            if (minnx is None or maxxx is None or minny is None or maxxy is None):
                minnx=min(np.ndarray.flatten(adata_sub.obsm['spatial'][:,0]))
                maxxx=max(np.ndarray.flatten(adata_sub.obsm['spatial'][:,0]))
                minny=min(np.ndarray.flatten(adata_sub.obsm['spatial'][:,1]))
                maxxy=max(np.ndarray.flatten(adata_sub.obsm['spatial'][:,1]))

            outpath2 = str(outpath) + "Plot_major_cluster_sample_"+str(ss)+"_sub"+str(fraction)+".pdf"

            #pdf = matplotlib.backends.backend_pdf.PdfPages(outpath+"/Plots/"+str(sample_name)+"/merFISH_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub_all.pdf")
            pdf = matplotlib.backends.backend_pdf.PdfPages(outpath2)
            fig, (ax1) = plt.subplots(1, 1, figsize=(20,10), gridspec_kw={'wspace':0.9})
            ax1.set_xlim(minnx,maxxx)
            ax1.set_ylim(minny,maxxy)
            sc.pl.spatial(adata_sub, spot_size=spot_size, color=['HDP_cluster'],palette=palette,title="Major Cell programs" , show=False, ax=ax1)
            pdf.savefig( fig )
            pdf.close()
            print(outpath2)


# import packageW
from gensim.test.utils import common_corpus, common_dictionary
from gensim.models import HdpModel
from gensim import corpora
import pandas as pd
import argparse
import os
import sys
import numpy as np
#import anndata as ad
import h5py      
import scanpy as sc
from scanpy import read_h5ad
import scipy as scipy
from .Cell_program_functions import *
from .Modules import *
from .Visualizations import *
#os.chdir("/gladstone/engelhardt/home/chwu/Research/HubID/Gensim/Scripts")



class CellProgramEstimator:
    def __init__(self, adata_sub=None, sample_name='Sample', assay='merFISH',
                 outdir='./', is_filter=True, K=100,T=100, gamma=1,alpha=1,kappa=1, random_state=0, fraction=0.05, ncol=5, spot_size_cluster=100, spot_size_activity=100):
        self.adata_sub = adata_sub
        self.sample_name = sample_name
        self.assay = assay
        self.outdir = outdir
        self.is_filter= is_filter
        self.K = K
        self.T = T
        self.gamma = gamma
        self.alpha = alpha
        self.kappa = kappa
        self.random_state = random_state
        self.fraction = fraction
        self.ncol = ncol
        self.spot_size_cluster = spot_size_cluster
        self.spot_size_activity = spot_size_activity


    def EstCellPrograms(self):
        adata_sub = self.adata_sub
        sample_name = self.sample_name
        assay = self.assay
        outdir = self.outdir
        is_filter= self.is_filter
        K = self.K
        T = self.T
        gamma = self.gamma
        alpha = self.alpha
        kappa = self.kappa
        random_state = self.random_state
        fraction = self.fraction
        ncol = self.ncol
        spot_size_cluster = self.spot_size_cluster
        spot_size_activity = self.spot_size_activity

        os.makedirs(outdir+"/"+str(sample_name)+"/Tables/", exist_ok=True)
        os.makedirs(outdir+"/"+str(sample_name)+"/Plots/", exist_ok=True)
        os.makedirs(outdir+"/"+str(sample_name)+"/Objects/", exist_ok=True)
    
        # extract info/ process matrices
        if adata_sub is not None:
            metadata=pd.DataFrame(adata_sub.obsm['spatial'])
            metadata.index=adata_sub.obs.index
            metadata.columns=["center_x", "center_y"]
            cellbygene=pd.DataFrame(adata_sub.X)
            cellbygene.columns=adata_sub.var.index
            cellbygene.index=adata_sub.obs.iloc[:,0]
            print(cellbygene)
            gene_names=adata_sub.var.iloc[:,0].tolist()
        else:
            sys.exit('adata_sub is None.')

        # identify cell programs
        print("Identifying Cell Programs...")
        topic_path=outdir+"/"+str(sample_name)+"/Tables/"
        cellprogram = IdentCellProgram(cellbygene=cellbygene,sample_name=sample_name, assay=assay,
                                        doc_tokenized_path=outdir+"/"+str(sample_name)+"/Objects/",
                                        BoW_corpus_path=outdir+"/"+str(sample_name)+"/Objects/",
                                        model_path=outdir+"/"+str(sample_name)+"/Objects/",
                                        prop_path=outdir+"/"+str(sample_name)+"/Tables/",
                                        topic_path=topic_path,
                                        K=K,T=T,
                                        gamma=gamma,alpha=alpha,kappa=kappa,random_state=random_state)
        Cell_Program = cellprogram['Cell_Program']
        Gene_Program = cellprogram['Gene_Program']
        print("Cell program estimation completed!")


        # get the main Cell program of each cell for plotting
        clusters_W = GetMain(Cell_Program)
        
        # Weight Gene_Program
        gene_prop=Gene_Program.sum(axis=1)
        sel_idx=gene_prop!=0
        Gene_Program_weighted=pd.DataFrame(index=Gene_Program.index[sel_idx], columns=Gene_Program.columns)
        for ii in range(Gene_Program.loc[sel_idx].shape[0]):
            rr = Gene_Program.loc[sel_idx].iloc[ii, :]
            m1 = rr.max()
            m2 = rr[rr != m1].max()
            mm = pd.Series(m1, index=rr.index)
            mm[rr == m1] = m2
            ns = rr * np.log((rr + 0.000000001) / (mm + 0.000000001))
            Gene_Program_weighted.iloc[ii,:]=ns

        Gene_Program_weighted.to_csv(topic_path+'/'+"/topic_info_weighted_"+str(sample_name)+"_"+str(assay)+".csv")
        
        # file printing out top genes
        rank20 = pd.DataFrame(index=Gene_Program_weighted.columns, columns=['0'])
        for ii in range(Gene_Program_weighted.shape[1]):
            tmp=Gene_Program_weighted.iloc[:, ii].sort_values(ascending=False).head(20)
            gene20=tmp.index
            ww20=pd.DataFrame(tmp)
            ww20=ww20.astype(float).round(4)
            ww20=ww20.replace({-0.0: 0.0})
            ww20=ww20.iloc[:,0]
            cbn = [f"{ww20.iloc[x]}*{gene20[x]}" for x in range(20)]
            rank20.iloc[ii,:]=pd.Series(" + ".join(cbn))
            
        rank20.to_csv(topic_path+"/top20_genes_weighted_"+str(sample_name)+"_"+str(assay)+".csv", index=True, header=False, quoting=1)

        
        # filter out non informative topics
        if is_filter==True:
            ndec=4
            idx_topic=Sel_topics(program=Gene_Program_weighted, ndec=ndec)
            Gene_Program_weighted_filtered=Gene_Program_weighted.iloc[:,idx_topic]
            print(str(len(idx_topic))+" topics are selected after filtering.")
            n_components=len(idx_topic)
            Cell_Program_filtered=Cell_Program.iloc[:, idx_topic]
            
            # save the filtered matrix
            Gene_Program_weighted_filtered.to_csv(topic_path+'/'+"/topic_info_weighted_filtered_"+str(ndec)+"_"+str(sample_name)+"_"+str(assay)+".csv")
            rank20.iloc[idx_topic,:].to_csv(topic_path+"/top20_genes_weighted_filtered"+str(ndec)+"_"+str(sample_name)+"_"+str(assay)+".csv", index=True, header=False, quoting=1)
        else:
            n_components=K
            Cell_Program_filtered=Cell_Program

        # plot each cell program
        print("Plotting...")


        tmp1=PlotMajCluster(adata=adata_sub,majcluster=clusters_W, outpath=outdir+"/"+str(sample_name)+"/Plots/"+str(assay)+"_"+str(sample_name)+"_major_topics_all.pdf",  spot_size=spot_size_cluster, sample_name=sample_name, label='CellProgram', ntopics=K)
        tmp2=PlotActivity(adata=adata_sub,prop=Cell_Program_filtered,outpath=outdir+"/"+str(sample_name)+"/Plots/"+str(assay)+"_"+str(sample_name)+"_sub"+str(fraction)+"_showall.pdf",
                         fraction=fraction, ncol=ncol, spot_size=spot_size_activity, random_state=random_state, sample_name=sample_name, label='CellProgram', ntopics=n_components)

        print("Cell programs estimation completed!")


# K=100,T=100, gamma=1,alpha=1,kappa=1

def main():
    # add argument
    parser = argparse.ArgumentParser(description='test script for the parser')
    parser.add_argument('-a', '--anndata', type=str,default=None, help='Path for the anndata')
    parser.add_argument('-s', '--sample_name', type=str, default="Sample" ,help='sample name')
    parser.add_argument('--assay', type=str, default="merFISH" ,help='Spatial transcriptomic platform')
    parser.add_argument('-o', '--outdir', type=str, default='./', help='Output path for tables and plots')
    parser.add_argument('--is_filter', type=bool, default=True ,help='Output the filtered topics for tables and plots')
    parser.add_argument('-k', '--K', type=int, default=100 ,help='Max number of topics in HDP')
    parser.add_argument('-t', '--T', type=int, default=100 ,help='Max number of tables in HDP')
    parser.add_argument('-g', '--gamma', type=int, default=1 ,help='parameter in HDP')
    parser.add_argument('--alpha', type=int, default=1 ,help='parameter in HDP')
    parser.add_argument('--kappa', type=int, default=1 ,help='parameter in HDP')
    parser.add_argument('-r', '--random_state', type=int, default=0 ,help='set the random state for NMF')
    parser.add_argument('-f', '--fraction', type=float, default=0.03 ,help='subset cells for plotting')
    parser.add_argument('-c', '--ncol', type=int, default=5 ,help='ncol for plotting')
    parser.add_argument('--spot_size_cluster', type=float, default=100 ,help='spot size for plotting')
    parser.add_argument('--spot_size_activity', type=float, default=100 ,help='spot size for plotting')

    # Parsing the arguments
    args = parser.parse_args()

    # Assign the arguments
    anndata_path = args.anndata
    sample_name =  args.sample_name
    assay = args.assay
    outdir = args.outdir
    is_filter= args.is_filter
    K = args.K
    T = args.T 
    gamma = args.gamma
    alpha = args.alpha
    kappa = args.kappa
    random_state = args.random_state
    fraction = args.fraction
    ncol = args.ncol
    spot_size_cluster = args.spot_size_cluster
    spot_size_activity = args.spot_size_activity

    # read input data
    adata_sub = sc.read_h5ad(anndata_path)
    # mtx = pd.read_csv(cellbygene_path)
    # mtx.index=mtx.iloc[:,0]# rownames(mtx)=mtx[,1]
    # mtx = mtx.iloc[:,1:]
 
    # est cell program
    cellprogram_estimator = CellProgramEstimator(adata_sub=adata_sub, sample_name=sample_name, assay=assay, outdir=outdir, is_filter=is_filter, K=K, T=T, gamma=gamma, alpha=alpha, kappa=kappa, random_state=random_state, fraction=fraction, ncol=ncol, spot_size_cluster=spot_size_cluster, spot_size_activity=spot_size_activity)
    # Call the est cell program method
    cellprogram_estimator.EstCellPrograms()

if __name__ == "__main__":
    main()

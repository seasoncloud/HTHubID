# import pacakges
import argparse
import numpy as np
import pandas as pd
#from sklearn.decomposition import NMF
import os
import scanpy as sc
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from .Landscape_functions import *
from .Modules import *
from .Visualizations import *
import matplotlib.backends.backend_pdf 

class LandscapeEstimator:
    def __init__(self, adata_sub=None, prop=None, program=None, sample_name='Sample', assay='merFISH', cellid=None, topicid=None,
                 neighbor_mode='NNeighbors', n_neighbors=500, eps=500, col_type='program', outdir='./', n_components=10, init='random',
                 random_state=0, alpha_W=0, fraction=0.05, ncol=5, spot_size=100):
        self.adata_sub = adata_sub
        self.prop = prop
        self.program = program
        self.sample_name = sample_name
        self.assay = assay
        self.cellid = cellid
        self.topicid = topicid
        self.neighbor_mode = neighbor_mode
        self.n_neighbors = n_neighbors
        self.eps = eps
        self.col_type = col_type
        self.outdir = outdir
        self.n_components = n_components
        self.init = init
        self.random_state = random_state
        self.alpha_W = alpha_W
        self.fraction = fraction
        self.ncol = ncol
        self.spot_size = spot_size

    def EstLandscape(self):
        adata_sub = self.adata_sub
        prop = self.prop
        program = self.program
        sample_name = self.sample_name
        assay = self.assay
        cellid = self.cellid
        topicid = self.topicid
        neighbor_mode = self.neighbor_mode
        n_neighbors = self.n_neighbors
        col_type = self.col_type
        eps = self.eps
        outdir = self.outdir
        n_components = self.n_components
        init = self.init
        random_state = self.random_state
        alpha_W = self.alpha_W
        fraction = self.fraction
        ncol = self.ncol
        spot_size = self.spot_size

        os.makedirs(outdir+"/"+str(sample_name)+"/Tables/", exist_ok=True)
        os.makedirs(outdir+"/"+str(sample_name)+"/Plots/", exist_ok=True)

        # select informative topics
        if program is not None:
            idx_topic=Sel_topics(program=program, ndec=4)
        else:
            idx_topic=np.array(range(prop.shape[1]))

        prop_filtered=prop[:,idx_topic]
        topicid=topicid[idx_topic]
        print(str(len(idx_topic))+" topics are selected after filtering.")
        #print(prop_filtered)
    
        # extract info/ process matrices
        metadata=pd.DataFrame(adata_sub.obsm['spatial'])
        metadata.index=adata_sub.obs.index
        metadata.columns=["center_x", "center_y"]
        cellbygene=pd.DataFrame(prop_filtered)
        
        if topicid is not None:
            cellbygene.columns=topicid
        else:
            cellbygene.columns=["Topic"+str(idx_topic.tolist()[x]) for x in range(cellbygene.shape[1])]
        
        if cellid is not None:
            cellbygene.index=cellid
        else:
            cellbygene.index=np.array(range(cellbygene.shape[0])).astype(str).tolist()
        print(cellbygene)

        # get neighborhood info
        print("Running GenNeighborhood...")
        cellbygene_smooth=GenNeighborhood(cellbygene=cellbygene, metadata=metadata, mode=neighbor_mode, n_neighbors=n_neighbors, sample_name=sample_name,col_type=col_type, eps=eps, outdir=outdir, assay=assay) 
        print("GenNeighborhood succeeded.")

        # identify landscapes
        print("Identifying Landscapes...")
        Landscapes = IdentLandscape(cellbygene_smooth=cellbygene_smooth, n_components=n_components, init=init, random_state=random_state, alpha_W=alpha_W)
        W=Landscapes['Cell_Landscape']
        H=Landscapes['Landscape_Program']
        print("Landscape estimation completed!")

        # save the landscapes
        if neighbor_mode=='NNeighbors':
            H.transpose().to_csv(outdir+"/"+str(sample_name)+"/Tables/"+"/Sectopic_info_"+str(sample_name)+"_"+str(assay)+"_n_neighbors_"+str(n_neighbors)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
            W.to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_prop_"+str(sample_name)+"_"+str(assay)+"_n_neighbors_"+str(n_neighbors)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
        elif neighbor_mode=='radius':
            H.transpose().to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_info_"+str(sample_name)+"_"+str(assay)+"_eps_"+str(eps)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
            W.to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_prop_"+str(sample_name)+"_"+str(assay)+"_eps_"+str(eps)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")


        # get the main landscape of each cell for plotting
        clusters_W = GetMain(W)

        # plot each landscapes
        print("Plotting...")
        if neighbor_mode=='NNeighbors':
            outpath = outdir+"/"+str(sample_name)+"/Plots/"+"/"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub"+str(fraction)+"_all.pdf"
        elif neighbor_mode=='radius':
            outpath = outdir+"/"+str(sample_name)+"/Plots/"+"/"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub"+str(fraction)+"_all.pdf"

        tmp=PlotActivity(adata=adata_sub,prop=W, outpath=outpath, fraction=fraction, ncol=ncol, spot_size=spot_size, random_state=random_state, sample_name=sample_name, label='SecTopic', ntopics=n_components)

        print("Landscape estimation completed!")




# def EstLandscape(adata_sub=None, prop=None, program=None, sample_name='Sample', assay='merFISH', cellid=None, neighbor_mode='NNeighbors', n_neighbors=500, eps=500, outdir='./', n_components=10, init='random',random_state=0, alpha_W=0, fraction=0.05, ncol=5, spot_size=100):
#     os.makedirs(outdir+"/"+str(sample_name)+"/Tables/", exist_ok=True)
#     os.makedirs(outdir+"/"+str(sample_name)+"/Plots/", exist_ok=True)

#     # select informative topics
#     idx_topic=Sel_topics(program=program, ndec=4)
#     prop_filtered=prop[:,idx_topic]
#     print(str(len(idx_topic))+" topics are selected after filtering.")
#     #print(prop_filtered)
    
#     # extract info/ process matrices
#     metadata=pd.DataFrame(adata_sub.obsm['spatial'])
#     metadata.index=adata_sub.obs.index
#     metadata.columns=["center_x", "center_y"]
#     cellbygene=pd.DataFrame(prop_filtered)
#     cellbygene.columns=["Topic"+str(idx_topic.tolist()[x]) for x in range(cellbygene.shape[1])]
#     cellbygene.index=cellid
#     print(cellbygene)

#     # get neighborhood info
#     print("Running GenNeighborhood...")
#     cellbygene_smooth=GenNeighborhood(cellbygene=cellbygene, metadata=metadata, mode=neighbor_mode, n_neighbors=n_neighbors, sample_name=sample_name,eps=eps, outdir=outdir, assay=assay) 
#     print("GenNeighborhood succeeded.")

#     # identify landscapes
#     print("Identifying Landscapes...")
#     Landscapes = IdentLandscape(cellbygene_smooth=cellbygene_smooth, n_components=n_components, init=init, random_state=random_state, alpha_W=alpha_W)
#     W=Landscapes['Cell_Landscape']
#     H=Landscapes['Landscape_Program']
#     print("Landscape estimation completed!")

#     # save the landscapes
#     if neighbor_mode=='NNeighbors':
#         H.transpose().to_csv(outdir+"/"+str(sample_name)+"/Tables/"+"/Sectopic_info_"+str(sample_name)+"_"+str(assay)+"_n_neighbors_"+str(n_neighbors)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
#         W.to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_prop_"+str(sample_name)+"_"+str(assay)+"_n_neighbors_"+str(n_neighbors)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
#     elif neighbor_mode=='radius':
#         H.transpose().to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_info_"+str(sample_name)+"_"+str(assay)+"_eps_"+str(eps)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")
#         W.to_csv(outdir+"/"+str(sample_name)+"/Tables/Sectopic_prop_"+str(sample_name)+"_"+str(assay)+"_eps_"+str(eps)+"_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+".csv")


#     # get the main landscape of each cell for plotting
#     clusters_W = GetMain(W)

#     # plot each landscapes
#     print("Plotting...")
#     if neighbor_mode=='NNeighbors':
#         outpath = outdir+"/"+str(sample_name)+"/Plots/"+"/"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub"+str(fraction)+"_all.pdf"
#     elif neighbor_mode=='radius':
#         outpath = outdir+"/"+str(sample_name)+"/Plots/"+"/"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+"_Sectopics_ncomp_"+str(n_components)+"_alpha_W_"+str(alpha_W)+"_sub"+str(fraction)+"_all.pdf"

#     tmp=PlotActivity(adata=adata_sub,prop=W, outpath=outpath, fraction=fraction, ncol=ncol, spot_size=spot_size, random_state=random_state, sample_name=sample_name, label='SecTopic', ntopics=n_components)

#     print("Landscape estimation completed!")


def main():
    # add argument
    parser = argparse.ArgumentParser(description='test script for the parser')
    parser.add_argument('-a', '--anndata', type=str, help='Path for the anndata')
    parser.add_argument('-p', '--prop', type=str, help='Path for the prop-- cell by program matrix')
    parser.add_argument('-d', '--program', type=str, help='Path for the program-- program by gene matrix')
    parser.add_argument('-s', '--sample_name', type=str, default="Sample" ,help='sample name')
    parser.add_argument('--assay', type=str, default="merFISH" ,help='Spatial transcriptomic platform')
    parser.add_argument('-m', '--neighbor_mode', type=str, default='NNeighbors', help='define neighborhood based on NNeighbors or radius')
    parser.add_argument('-n', '--n_neighbors', type=int, default=500 ,help='Set number of nearest neighbors when neighbor_mode is NNeighbors')
    parser.add_argument('-e', '--eps', type=int, default=500 ,help='Set radius when neighbor_mode is radius')
    parser.add_argument('--col_type', type=str, default="program" ,help='gene (1st level) or program (2nd level)')
    parser.add_argument('-o', '--outdir', type=str, default='./', help='Output path for tables and plots')
    parser.add_argument('-k', '--n_components', type=int, default=10 ,help='number of landscapes')
    parser.add_argument('-i', '--init', type=str, default="random" ,help='NMF initiation mode')
    parser.add_argument('-r', '--random_state', type=int, default=0 ,help='set the random state for NMF')
    parser.add_argument('-W', '--alpha_W', type=float, default=0.01 ,help='parameter used in NMF')
    parser.add_argument('-f', '--fraction', type=float, default=0.05 ,help='subset cells for plotting')
    parser.add_argument('-c', '--ncol', type=int, default=5 ,help='ncol for plotting')
    parser.add_argument('--spot_size', type=float, default=100 ,help='spot size for plotting')

    # Parsing the arguments
    args = parser.parse_args()

    # Assign the arguments
    anndata_path = args.anndata
    prop_path = args.prop
    program_path = args.program
    sample_name =  args.sample_name
    assay = args.assay
    neighbor_mode = args.neighbor_mode
    n_neighbors = args.n_neighbors
    eps = args.eps
    col_type = args.col_type
    outdir = args.outdir
    n_components = args.n_components
    init = args.init
    random_state = args.random_state
    alpha_W = args.alpha_W
    fraction = args.fraction
    ncol = args.ncol
    spot_size = args.spot_size


    # read input data
    adata_sub = sc.read_h5ad(anndata_path)
    prop=pd.read_csv(prop_path)
    cellid=prop.iloc[:,0]
    prop=prop.iloc[:,1:]
    topicid=prop.columns
    prop=np.array(prop)
    
    # if not showing the weights for programs, use all programs provided
    if program_path is not None:
        program=pd.read_csv(program_path, index_col=0)
    else:
        program=None
    
    # remove NA
    rm_id=np.union1d(np.where((np.isnan(prop).any(axis=1))), np.where(pd.DataFrame(adata_sub.obsm['spatial']).isna().any(axis=1)))
    sel_id=np.delete(np.array(range(prop.shape[0])),rm_id)
    adata_sub=adata_sub[sel_id,:]
    prop=prop[sel_id,:]
    cellid=cellid.loc[sel_id]
 
    # run EstLandscape
    landscape_estimator = LandscapeEstimator(adata_sub=adata_sub, prop=prop, program=program, sample_name=sample_name, assay=assay, cellid=cellid, topicid=topicid, neighbor_mode=neighbor_mode, n_neighbors=n_neighbors, eps=eps, col_type=col_type, outdir=outdir, n_components=n_components, init=init, random_state=random_state,alpha_W=alpha_W, fraction=fraction, ncol=ncol, spot_size=spot_size)
    # Call the EstLandscape method
    landscape_estimator.EstLandscape()

if __name__ == "__main__":
    main()

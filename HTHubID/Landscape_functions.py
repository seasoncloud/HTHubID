# import libraries
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors




def GenNeighborhood(cellbygene=None, metadata=None, mode='NNeighbors', n_neighbors=500, eps=500,sample_name="sample",col_type='program',outdir='./', assay="merFISH"):
    cell_idx=cellbygene.index
    topic_names=cellbygene.columns
    cellbygene=np.array(cellbygene)
    if mode=='NNeighbors': ## revised from Jack Demaray's codes
        if os.path.exists(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+".csv"):
            cellbygene_smooth= pd.read_csv(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+".csv", index_col=0)
            print("Read cellbygene_smooth from "+ outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+".csv")
        else:
            fit = NearestNeighbors(n_neighbors=n_neighbors).fit(metadata.loc[:, ["center_x", "center_y"]])
            m = fit.kneighbors(metadata.loc[:, ["center_x", "center_y"]])
            m = m[0], m[1]
            args = m[0].argsort(axis=1)
            add = np.arange(m[1].shape[0])*n_neighbors
            sorted_indices = m[1].flatten()[args+add[:, None]]
            cellbygene_smooth=np.zeros((cellbygene.shape[0],cellbygene.shape[1]))
            for ii in range(sorted_indices.shape[0]):
                rr=cellbygene[sorted_indices[ii,:],:]
                cellbygene_smooth[ii,:]=rr.mean(axis=0)
            cellbygene_smooth=pd.DataFrame(cellbygene_smooth)
            cellbygene_smooth.index=cell_idx
            cellbygene_smooth.columns=topic_names
            cellbygene_smooth.to_csv(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_n_neighbors_"+str(n_neighbors)+".csv")
    elif mode=='radius':
        if os.path.exists(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+".csv"):
            cellbygene_smooth= pd.read_csv(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+".csv",  index_col=0)
            print("Read cellbygene_smooth from "+outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+".csv")
        else:
            nn = NearestNeighbors(radius=eps)
            nn.fit(metadata.loc[:, ["center_x", "center_y"]])
            distances,neighbors = nn.radius_neighbors(metadata.loc[:, ["center_x", "center_y"]], return_distance=True)
            sorted_data = [(np.array(neigh)[np.argsort(dist)], np.array(dist)[np.argsort(dist)])
                for neigh, dist in zip(neighbors, distances)]
            # Unzip the sorted results
            sorted_indices, sorted_distances = zip(*sorted_data)
            cellbygene_smooth=np.zeros((cellbygene.shape[0],cellbygene.shape[1]))
            for ii in range(len(sorted_indices)):
                rr=cellbygene[sorted_indices[ii],:]
                cellbygene_smooth[ii,:]=rr.mean(axis=0)
            cellbygene_smooth=pd.DataFrame(cellbygene_smooth)
            cellbygene_smooth.index=cell_idx
            cellbygene_smooth.columns=topic_names
            cellbygene_smooth.to_csv(outdir+"/"+str(sample_name)+"/Tables/neighborhoodby"+str(col_type)+"_"+str(assay)+"_"+str(sample_name)+"_eps_"+str(eps)+".csv")
            
            # plot number of neighbors within eps radius
            nneighbors=[]
            for ii in sorted_indices:
                nneighbors.append(len(ii))
            plt.hist(nneighbors, bins=50)
            plt.savefig(outdir+"/"+str(sample_name)+"/Tables/"+str(sample_name)+"_"+str(assay)+"_nneighbors_distn_eps_"+str(eps)+".pdf", format="pdf")
            plt.show() 
            plt.close()


    else:
        print("Please specify a valid neighborhood mode")
        raise SystemExit()
 
    return cellbygene_smooth

def IdentLandscape(cellbygene_smooth=None, n_components=10, init="random", random_state=0, alpha_W=0.01):
    model = NMF(n_components=n_components, init=init, random_state=random_state, alpha_W=alpha_W)
    W = model.fit_transform(cellbygene_smooth) # cell * component
    H = model.components_   # component * programs
    W = pd.DataFrame(W)
    W.index=cellbygene_smooth.index
    W.columns=["SecTopic"+str(x) for x in range(W.shape[1])]
    H=pd.DataFrame(H)
    H.index=["SecTopic"+str(x) for x in range(W.shape[1])]
    #H.columns=["Topic"+str(idx_topic.tolist()[x]) for x in range(cellbygene_smooth.shape[1])]
    H.columns=cellbygene_smooth.columns
    print(W)
    print(H)
    results = {'Cell_Landscape':W, 'Landscape_Program':H}
    return results



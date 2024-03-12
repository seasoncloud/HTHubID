# import packageW
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import geopandas as gpd
import random
import seaborn as sns
import os
import sys
#os.chdir("/gladstone/engelhardt/home/chwu/Research/HubID/Gensim/Scripts")



def main():
    # add argument
    parser = argparse.ArgumentParser(description='test script for the parser')
    parser.add_argument('-n', '--sample_name', type=str,default="Sample", help='Sample name for the object')
    parser.add_argument('-a', '--assay', type=str,default="Sample", help='Assay for the object')
    parser.add_argument('-m', '--input_mtx', type=str,default= None , help='Path for the matrix')
    parser.add_argument('-d', '--input_meta', type=str,default= None , help='Path for the metadata')
    parser.add_argument('-i', '--index_col_mtx', type=int, default= None , help='Column(s) to use as row label(s) or None for the matrix')
    parser.add_argument('-j', '--index_col_meta', type=int, default= None , help='Column(s) to use as row label(s) or None for the matrix')
    parser.add_argument('-s1', '--spatial_col1', type=str, default='center_x', help='column names for spatial info')
    parser.add_argument('-s2', '--spatial_col2', type=str, default='center_y', help='column names for spatial info')
    parser.add_argument('-e', '--seed', type=int, default=0 ,help='set the random state for shuffling')
    parser.add_argument('-f1', '--min_genes', type=int, default=0 ,help='filter cells based on the number of genes')
    parser.add_argument('-f2', '--min_counts', type=int, default=0 ,help='filter the cells based counts')
    parser.add_argument('-p', '--outpath', type=str, default='./' ,help='output path for the object')


    # Parsing the arguments
    args = parser.parse_args()

    # Assign the arguments
    sample_name = args.sample_name
    assay = args.assay
    input_mtx = args.input_mtx
    input_meta = args.input_meta
    index_col_mtx = args.index_col_mtx
    index_col_meta = args.index_col_meta
    spatial_col = [args.spatial_col1, args.spatial_col2]
    seed = args.seed
    min_genes = args.min_genes
    min_counts = args.min_counts
    outpath = args.outpath


    # if index_col_mtx != "None":
    #     index_col_mtx = int(index_col_mtx)

    # if index_col_meta != "None":
    #     index_col_meta = int(index_col_meta) 

    # set seed
    #random.seed(seed)

    # read input data
    metadata=pd.read_csv(input_meta, index_col=index_col_meta) # 11 for baysor
    #metadata.set_index("cell_id", inplace=True)
    # There is ONE cell with NaN as its x, y coordinates; I just drop that and move on
    metadata.dropna(axis=0, how="any", subset= spatial_col, inplace=True)

    # shuffle
    np.random.seed(seed)
    idx=np.random.permutation(metadata.index)
    metadata=metadata.reindex(idx)
    
    cellbygene = pd.read_csv(input_mtx, index_col=index_col_mtx)
    cellbygene = cellbygene.reindex(idx)

    # Select the 'center_x' and 'center_y' columns into the 'spatial' DataFrame
    coordinates = metadata[spatial_col]
    coordinates.index.name = 'cell'
    coordinates.index = coordinates.index.map(str)

    # create andata object
    adata = sc.AnnData(cellbygene, 
        cellbygene.index.to_frame(), 
        cellbygene.columns.to_frame())  ## index to string
    adata.obsm['spatial']=coordinates

    # filtering
    adata_sub=adata#sc.pp.subsample(adata, fraction=0.1, copy=True, random_state=0)
    if min_genes>0:
        sc.pp.filter_cells(adata_sub, min_genes=min_genes)
    elif min_counts>0:
        sc.pp.filter_cells(adata_sub,  min_counts= min_counts)
    
    # save the object
    if outpath == './':
        outpath = outpath+str(sample_name)+"_"+str(assay)+"_adata_mingenes"+str(min_genes)+"_mincounts"+str(min_counts)+".h5ad"
    adata_sub.obs.rename(columns={0: "index"}, inplace=True)
    adata_sub.var.rename(columns={0: "gene"}, inplace=True)
    adata_sub.write_h5ad(outpath)

    print("Object successfully saved!")

    
if __name__ == "__main__":
    main()

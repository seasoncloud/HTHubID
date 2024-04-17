# import libraries
from gensim.test.utils import common_corpus, common_dictionary
from gensim.models import HdpModel
from gensim import corpora
import pandas as pd
import os
import sys
import numpy as np
import pandas as pd
#import anndata as ad
import h5py      
import scanpy as sc
from scanpy import read_h5ad
import scipy as scipy
import os
import pickle
import time
import matplotlib.pyplot as plt
#from test_functions import *
import matplotlib.backends.backend_pdf
from sklearn.neighbors import NearestNeighbors



def IdentCellProgram(cellbygene=None, sample_name='Sample', assay='assay' ,doc_tokenized_path='./', BoW_corpus_path='./', model_path='./', prop_path='./',topic_path='./', K=100, T=100, gamma=1, alpha=1, kappa=1, random_state=100):
    gene_names = cellbygene.columns
    cellid = cellbygene.index
    dat=np.array(cellbygene)
    # doc_tokentize
    #print("Generating doc_tokenized!")
    if os.path.isfile(doc_tokenized_path+'/'+str(sample_name)+"_"+str(assay)+'_doc_tokenized.pkl'):
        print("Importing doc_tokenized!")
        with open(doc_tokenized_path+'/'+str(sample_name)+"_"+str(assay)+'_doc_tokenized.pkl', 'rb') as f2:
            doc_tokenized = pickle.load(f2)
        print("import from "+doc_tokenized_path+'/'+str(sample_name)+"_"+str(assay)+'_doc_tokenized.pkl' )
    else:
        print("Generating doc_tokenized!")
        doc_tokenized=[]
        for ii in range(dat.shape[0]):
            rr=dat[ii,:]
            rr=rr.astype(int)
            ind=0
            doc=[]
            for jj in rr:
                if jj>0:
                    for tt in range(jj):
                        doc.append(str(gene_names[ind]))
                ind+=1
            doc_tokenized.append(doc)
        with open(doc_tokenized_path+'/'+str(sample_name)+"_"+str(assay)+'_doc_tokenized.pkl', 'wb') as ff:
            pickle.dump(doc_tokenized, ff)
        print("doc_tokenized has been generated and stored in "+ doc_tokenized_path+'/'+str(sample_name)+"_"+str(assay)+'_doc_tokenized.pkl')
    
    # get BoW_corpus
    dictionary = corpora.Dictionary()
    
    if os.path.isfile(BoW_corpus_path+'/'+str(sample_name)+'_'+str(assay)+'_BoW_corpus_1.pkl'):
        with open(BoW_corpus_path+'/'+str(sample_name)+'_'+str(assay)+'_BoW_corpus_1.pkl', 'rb') as f2:
            BoW_corpus = pickle.load(f2)
        with open(BoW_corpus_path+'/'+str(sample_name)+'_'+str(assay)+'_BoW_corpus_2.pkl', 'rb') as f3:
            BoW_corpus2 = pickle.load(f3)
        BoW_corpus= BoW_corpus + BoW_corpus2
    else:
        BoW_corpus = [dictionary.doc2bow(doc, allow_update=True) for doc in doc_tokenized]
        with open(BoW_corpus_path+'/'+str(sample_name)+'_'+str(assay)+'_BoW_corpus_1.pkl', 'wb') as f2:
            pickle.dump(BoW_corpus[0:int(len(BoW_corpus)/2)], f2)
        with open(BoW_corpus_path+'/'+str(sample_name)+'_'+str(assay)+'_BoW_corpus_2.pkl', 'wb') as f3:
            pickle.dump(BoW_corpus[int(len(BoW_corpus)/2):len(BoW_corpus)], f3)

    # train HDP model
    print("Training the HDP model.")
    if os.path.isfile(model_path+'/'+'hdp_'+str(sample_name)+'_'+str(assay)+'.model'):
        hdp = HdpModel.load(model_path+'/'+'hdp_'+str(sample_name)+'_'+str(assay)+'.model')
        print("import from "+model_path+'/'+'hdp_'+str(sample_name)+'_'+str(assay)+'.model')
        total=0
    else:
        t0 = time.time()
        hdp = HdpModel(BoW_corpus, dictionary, T=T,K=K, gamma=gamma, alpha=alpha, kappa=kappa, random_state=random_state)
        t1 = time.time()
        total=t1-t0
        hdp.save(model_path+'/'+'hdp_'+str(sample_name)+'_'+str(assay)+'.model')
        print("HDP model has been trained and stored in "+model_path+'/'+'hdp_'+str(sample_name)+'_'+str(assay)+'.model')
    
    # get prop prop
    result_dist=hdp[BoW_corpus]
    # get prop
    if os.path.isfile(prop_path+"/"+"/prop_"+str(sample_name)+"_"+str(assay)+".csv"):
        prop=pd.read_csv(prop_path+"/"+"/prop_"+str(sample_name)+"_"+str(assay)+".csv")
        prop=prop.iloc[:,1:]
        topicid=prop.columns
        #prop=np.array(prop)
        if len(topicid) != Ｋ:
            sys.exit("nrows of imported prop is not equal to the number of topics!")
    else: 
        prop = np.empty((dat.shape[0],T), float)#pd.DataFrame([])
        for ii in range(dat.shape[0]):
            rr=np.zeros(T)
            result=np.array(result_dist[ii])
            #print(result)
            if len(result)>0:
                rr[result[:,0].astype(int)]=result[:,1]
            rr=np.reshape(rr,(1,T)) 
            prop[ii,:]=rr
        # save prop
        prop= pd.DataFrame(prop)
        prop.index=cellid
        prop.columns=["Program"+str(x) for x in range(prop.shape[1])]
        prop.to_csv(prop_path+"/"+"/prop_"+str(sample_name)+"_"+str(assay)+".csv")


    # get topic info
    topic_info=hdp.get_topics()
    df_topic_info=pd.DataFrame(topic_info)
    df_topic_info.columns=hdp.id2word.id2token.values()
    df_topic_info=df_topic_info.transpose()
    df_topic_info.to_csv(topic_path+'/'+"/topic_info_"+str(sample_name)+"_"+str(assay)+".csv")
    
    # retun results
    results = {'Cell_Program':prop, 'Gene_Program':df_topic_info, 'HDP_time':total}
    return results


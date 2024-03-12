import pandas as pd
import os
import sys
import numpy as np   




def GetMain(W=None):
    prop_W=np.array(W)
    clusters_W=[]
    for ii in range(prop_W.shape[0]):
        clusters_W.append(str(prop_W[ii,:].argmax()))
    return clusters_W


def Sel_topics(program=None, ndec=4):
    wmax=[]
    for ii in range(program.shape[1]):
        tmp=program.iloc[:, ii].sort_values(ascending=False).head(2)
        gene20=tmp.index
        ww20=pd.DataFrame(tmp)
        ww20=ww20.astype(float).round(ndec)
        ww20=ww20.replace({-0.0: 0.0})
        ww20=ww20.iloc[0,0]
        wmax.append(ww20)

    wmax = np.array(wmax)  #
    idx_topic = np.where((wmax != 0) & (wmax>0))[0]

    return idx_topic


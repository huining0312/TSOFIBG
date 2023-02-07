from cmath import log
from math import log2
import pandas as pd
import numpy as np
# --- NDCG function 找大的比較好的 --- #
def get_dcg(y_true,y_pred,k):
    df = pd.DataFrame({"y_true":y_true,"y_pred":y_pred}) #降序排列
    y1 = df.sort_values(by="y_pred", ascending=False) 
    dcg = 0
    for i in range(1,k+1):
        a = y1.iloc[i-1]["y_true"] / log2(i+1)
        dcg = dcg + a
    return dcg

def get_ndcg(y_true,y_pred,k):
    dcg=get_dcg(y_true,y_pred,k)
    idcg=get_dcg(y_true,y_true,k)
    ndcg=dcg/idcg
    return ndcg

def get_ndcg_mean(y_true,y_pred,k):
    a=0
    for i in range(1,k+1):
        b=get_ndcg(y_true,y_pred,k=i)
        a=a+b
    return (a/k)


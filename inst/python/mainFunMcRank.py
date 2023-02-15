# function for trait
def main_fun_Mcrank(snp_matrix,TRsize,pheno,trainsetCandi,RFC):
  xrange = range
  gblup_Mcrank_NDCG = list()
  
  for trait in pheno.columns.values[1:]:
    
    #Mcrank_NDCG = pandas.DataFrame(index=TRsize, columns=['RE_NDCG@k=1','RE_NDCG@k=5','RE_NDCG@k=10','RE_meanNDCG@k=10']) 
    Mcrank_res = pandas.DataFrame(index = snp_matrix.index, columns=[str(i) for i in TRsize])
    
    for size in range(0,len(TRsize)):
      X = snp_matrix.loc[snp_matrix.index.isin(trainsetCandi[size])]
      Y = pheno[trait]
      Y2 = Y.loc[pheno.index.isin(trainsetCandi[size])]
      Y2=Y2.values.ravel()  ###轉換變一維
      mcrank = RFC(n_estimators=300,max_features=0.6,max_depth=10,criterion="entropy")
      model2 = OrdinalMcRank(mcrank) ##random forest model
      model2.fit(X,Y2)
      Y_pred2=model2.predict(snp_matrix)
      Mcrank_res[str(TRsize[size])] = Y_pred2
      #g_true = pheno[trait]

      #g_true = g_true.values.ravel()
      #result1=round(get_ndcg(g_true,Y_pred2,k=1),4)
      #result2=round(get_ndcg(g_true,Y_pred2,k=5),4)
      #result3=round(get_ndcg(g_true,Y_pred2,k=10),4)
      #result4=round(get_ndcg_mean(g_true,Y_pred2,k=10),4)
      #AA=[result1,result2,result3,result4]
      # print(AA)
      
      #Mcrank_NDCG.iloc[size,range(4)] = AA
      
    gblup_Mcrank_NDCG.append(Mcrank_res)
      
  return gblup_Mcrank_NDCG
  

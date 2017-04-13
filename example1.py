#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 00:10:20 2017

@author: jmfuseli
"""

import numpy as np
from pandas import DataFrame 
from sklearn import linear_model

#YORF stands for yeast open reading frame

df = DataFrame.from_csv("file.txt", sep="\t",header=1)
#see the dimensions of the df

#extract row names
row_names = list(df.index)

#these are the column indices to separate the df by the 3 timecourses
#ex. [0:9], [9:18], [18:27]
a = [0,9,18]
b = [9,18,27]

def get_halflife(dataframe):
    #drop zeroes since log(0)= -inf
    dataframe = dataframe.loc[(dataframe!=0)]
    #get the natural log of data
    dataframe_log = np.log(dataframe)
    #drop NaNs
    dataframe_noNaN = dataframe_log.dropna(axis=0, how='all')
    #check length of df_noNaN since at least one YORF has no data
    dataframe_noNaN_len = len(dataframe_noNaN.index)
    if dataframe_noNaN_len != 0:
       reg = linear_model.LinearRegression()
       x = np.matrix(dataframe_noNaN.index).reshape(dataframe_noNaN_len,1)
       reg.fit(x,dataframe_noNaN)
       slope = reg.coef_[0]
       if slope != 0:
          half_life=(np.log(2)/slope)*-1
          return half_life
    else:
       pass


def write(tuple_tr, filename):
   sub_list = []
   for tr in range(len(tuple_tr)):
      sub_list.append(tuple_tr[tr][0] + "\n")
   #write to a file 
   f = open(filename,"w")
   f.writelines(sub_list)
   f.close()   
   
d = {} #empty dict to store output YORF and its halflife

for i in range(len(row_names)):
   halflives = []
   for j in range(len(a)):
      df1 = df.loc[row_names[i]] 
      #check num of dimensions since the same transcript appears more than once   
      df1_len = len(df1.shape)
      if df1_len > 1:
          sub_hl = []
          #takes avg of the halflives if the YORF appears more than once
          for k in range(df1_len):
              df2 = df1.iloc[k][a[j]:b[j]]
              hl2 = get_halflife(df2)
              if hl2 is not None:
                 sub_hl.append(hl2)
          one_hl = sum(sub_hl)/len(sub_hl)
          halflives.append(one_hl)
      else: 
          df2 = df1[a[j]:b[j]]
          hl = get_halflife(df1)
          if hl is not None: #check if halflives is empty
             halflives.append(hl)  
   #calculate avg halflife
   #check len of halflives since some YORFs have no data
   if len(halflives) != 0:
      avg_hl = sum(halflives)/len(halflives)
      if row_names[i] not in d:
         d[row_names[i]] = avg_hl
   else:
       pass

#sort the dictionary based on the halflife or value   
sorted_YORFs = sorted(d.items(), key = lambda tup: tup[1])  

#get length of sorted_YORFs
result_len = len(sorted_YORFs)
print "The number of genes after halflife analysis is %i." %result_len

#take the bottom 10% of genes which is about 616 genes
subset_low = sorted_YORFs[:616]

#top 10% of 6163 genes, so 616 from the top
subset_high = sorted_YORFs[5546:6162]

#write top 10% and bottom 10% YORFs to a file
write(subset_low, "bottom10.txt")
write(subset_high, "top10.txt")   



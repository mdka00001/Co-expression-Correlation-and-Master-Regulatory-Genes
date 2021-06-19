import pandas as pd
import csv
import operator
import math
import numpy as np
from scipy import stats
import sys
import statistics
from matplotlib import pyplot
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import r2_score

class Correlation:
    def __init__(self, path):
        with open(path, "r") as handle:
            table=csv.reader(handle, delimiter="\t")
            df1=pd.DataFrame(table)
            df=df1.T
            df.drop(df.columns[[0]], axis=1, inplace=True)
            new_header = df.iloc[0]
            df_head = df[1:]
            df_head.columns=new_header

            x_replace=df_head.replace('',0, regex=True)
            
            df_3=x_replace.astype(float)
            array_4=(df_3.T).to_numpy()
            

            x=[]

            diff=[]

            for i in df_3.columns:
                x.append(df_3[i].mean())
            for i in range(0,len(x)):
                d=array_4[i]-x[i]
                diff.append(d)
            
            
            diff_df=pd.DataFrame(diff)
            diff_df2=diff_df.T
            diff_df2.columns=new_header
            
            final_list={}

            list_itgb={}
            for i in diff_df2.columns:
                for j in diff_df2.columns:
                    if i != j:
                        x=diff_df2[i]
                        y=diff_df2[j]

                        xx=x**2
                        yy=y**2

                        z=x*y

                        val=sum(z)/(math.sqrt(sum(xx))*math.sqrt(sum(yy)))
                        
                        name=i+" "+j
                        
                        final_list[name]=val
                        
                        if i == "ITGB2":
                            list_itgb[name]=val
            print(list_itgb)
                            
            
            temp1=[]
            res1 = dict()
            for key, val in final_list.items():
                if val not in temp1:
                    temp1.append(val)
                    res1[key] = val
            
            corr=list(final_list.values())

            """ITGB2"""
            itgb_val=list(list_itgb.values())
            

            sorted_itgb ={}

            for key, value in sorted(list_itgb.items(), key=lambda item: item[1]):
                data = {key: value}
                sorted_itgb.update(data)

            high=max(itgb_val)
            low=min(itgb_val)
            near_zero=min(i for i in sorted_itgb if  sorted_itgb[i] != 0.000000)

            def get_key(val):
                for key, value in list_itgb.items():
                    if val == value:
                        return key
 
                return "key doesn't exist"

            stdoutOrigin=sys.stdout
            sys.stdout=open("log2.txt","w")
            
            print("Highest correlation: %s Pair" % get_key(high))
            print("Lowest correlation: %s Pair" % get_key(low))
            print("Near to zero correlation: %s Pair" % near_zero)
            

            sys.stdout.close()
            sys.stdout=stdoutOrigin

            sns.displot(corr, bins=50, kde=True, edgecolor="k")
            plt.savefig("Histogram.png")
            plt.close()

            df_reg=x_replace[["ITGB2", "TRIM25", "CKMT1", "ABHD15"]]
            df_reg2=df_reg.reset_index()
            df_reg2.drop(df_reg2.columns[[0]], axis=1, inplace=True)
            df_reg3=df_reg2.astype(float)
            
            
            for i in df_reg3.columns:
                for j in df_reg3.iloc[:,1:]:
                    if i == "ITGB2":

                        
                        sns.regplot(x=i, y=j, data=df_reg3)

                        if j=="TRIM25":
                            plt.savefig("TRIM25.png")
                            plt.close()
                        if j=="CKMT1":
                            plt.savefig("CKMT1.png")
                            plt.close()
                        if j=="ABHD15":
                            plt.savefig("ABHD15.png")
                            plt.close()
            
           
            w=df_reg3["ITGB2"].to_numpy()
            x=df_reg3["TRIM25"].to_numpy()
            y=df_reg3["CKMT1"].to_numpy()
            z=df_reg3["ABHD15"].to_numpy()

            stdoutOrigin=sys.stdout
            sys.stdout=open("log1.txt","w")
            print("ITGB2 TRIM25 r-score: %f" % r2_score(w,x))
            print("ITGB2 ABDH15 r-score: %f" % r2_score(w,z))
            print("ITGB2 CKMT1 r-score: %f" % r2_score(w,y))
            
            sys.stdout.close()
            sys.stdout=stdoutOrigin
           

                        


a=Correlation(r"D:\BIII\Assignment_7\BI3_A7_Karim_Elamaldeniya\gene_expression.csv")
print(a)

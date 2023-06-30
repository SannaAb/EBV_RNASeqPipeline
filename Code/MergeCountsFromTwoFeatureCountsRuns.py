import pandas as pd
import sys

firstcounts=sys.argv[1]
unnassignedcounts=sys.argv[2]
data1=pd.read_csv(firstcounts,sep="\t",comment='#',skiprows=1)
data2=pd.read_csv(unnassignedcounts,sep="\t",comment='#',skiprows=1)
res=pd.merge(data1,data2,how="left",on="Geneid")
columnstosave=[0,1,2,3,4,5,6,12]
res=res.iloc[:,columnstosave]
c1=res.columns[-2]
c2=res.columns[-1]
res[c2] = res[c2].fillna(0)
res[c2] = res[c2].astype(int)
c2new=res.columns[-1].replace("_unassigned_sorted.bam","")
res[c2new]=res[c1]+res[c2]
columnstosave=[0,1,2,3,4,5,8]
res=res.iloc[:,columnstosave]
outfile=c2new+"_Merge_unassigned_FeatureCounts.txt"
res.to_csv(outfile,index=False,sep="\t",header=True)
print("Fin")

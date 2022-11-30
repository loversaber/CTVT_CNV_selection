import pandas as pd

def generate_matrix(df):
 df=df.copy()
 #df[""]=df["CHROM"]+"_"+df["dataset_segment_id"]
 df1=df[df["anc_index"]==91]
 df2=df[df["desc_index"]==91]

 df1.loc[:,"changes"]=df1["anc_cnv"]-2
 df1["abs_changes"]=df1["changes"].abs()
 df11=df1[df1["changes"]>0]
 df11.loc[:,"type"]="GAIN"
 df12=df1[df1["changes"]<0]
 df12.loc[:,"type"]="LOSS"

 df2.loc[:,"changes"]=df2["desc_cnv"]-2
 df2.loc[:,"abs_changes"]=df2["changes"].abs()
 df21=df2[df2["changes"]>0]
 df21["type"]="GAIN"
 df22=df2[df2["changes"]<0]
 df22["type"]="LOSS"
 df3=pd.concat([df11,df12,df21,df22])[["CHROM","dataset_segment_id","abs_changes","type"]]
 df4=pd.merge(df3,df_coordinate,on=["CHROM","dataset_segment_id"],how="left")
 df4=df4.copy()
 df4["WIDTH"]=df4["END"]-df4["START"]
 df4=df4.drop_duplicates().reset_index(drop=True)
 df4.to_csv("Anc_test.csv",sep=",",header=True,index=False)
 df41=df4[df4["abs_changes"]==1]
 df42=df4[df4["abs_changes"]>1]
 df_end=pd.DataFrame(columns=['CHROM', 'START', 'END', 'WIDTH', 'EVENT TYPE'])
 ll=[]#'CNV EVENT ID'
 for index,row in df4.iterrows():
  abs_diff=row["abs_changes"]
  row_values=[row["CHROM"],row["START"],row["END"],row["WIDTH"],row["type"]]
  for i in range(abs_diff):
   df_end.loc[len(df_end)]=pd.Series(row_values,index=df_end.columns)
 df_end=df_end.copy().reset_index(drop=True)
 df_end=df_end.rename(columns={"type":"EVENT TYPE"})
 df_end['CNV EVENT ID']=df_end.index+1
 df_end=df_end[['CNV EVENT ID','CHROM', 'START', 'END', 'WIDTH', 'EVENT TYPE']]
 df_end.to_csv("Anc_chrom_map.csv",sep="\t",header=True,index=False)

if __name__=="__main__":
 df=pd.read_csv("All_changes_new_type.csv",sep=",",header=0)#wrong events have been removed
 df_coordinate=pd.read_csv("All_Gain_Loss_counts.csv",sep=",",header=0)[["CHROM","dataset_segment_id","START","END"]]
 generate_matrix(df)
 

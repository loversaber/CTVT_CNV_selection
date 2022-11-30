import pandas as pd

df=pd.read_excel("Events_Matrix_v1_new_copy_correct_copy2.xlsx",sheet_name="Events_Matrix_v1_new_copy_corre",header=0)

df=df.copy()
#df["CNV EVENT ID"]=df["ID"].str.split(r"L|G",expand=True)[0]
df["WIDTH"]=df["END"]-df["START"]+1
df=df.rename(columns={"DIRECTION":"EVENT TYPE"})
df=df.fillna("N")
df1=df.copy()
df1["CNV EVENT ID"]=df1.index+1
df_end1=df1[["CNV EVENT ID","CHROM","START","END","WIDTH","EVENT TYPE"]]
df_end1.to_csv("Events_correct_with_wrong_v1.txt",sep="\t",header=True,index=False)

df=df[~df["NOTE"].str.contains(r'wrong')].reset_index(drop=True)
print(df.index)
df["CNV EVENT ID"]=df.index+1

df_end=df[["CNV EVENT ID","CHROM","START","END","WIDTH","EVENT TYPE"]]
df_end.to_csv("Events_correct_v1.txt",sep="\t",header=True,index=False)


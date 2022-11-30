import pandas as pd
import numpy as np

###check the Back mutation
def get_every_path():
 d={}
 for i in range(1,len(all_ordered_samples)+2):#plus REF
  '''
  if i==1:
   print(f"leaf index {i},'REF'")
  else:
   print(f"leaf index,label:{i} {all_ordered_samples[i-2]}")
  '''
  s_path=""
  anc_i=df_edge[df_edge["desc"]==i]["anc"].values[0]
  s_path+=str(i)+"_"+str(anc_i)
  while anc_i!=47:###cannot use integer value to compare
   anc_i2=df_edge[df_edge["desc"]==anc_i]["anc"].values[0]
   s_path+="_"+str(anc_i2)
   anc_i=anc_i2
  print(s_path)
  d[i]=s_path
 return d

def modify_type(df_site_pattern_new,opt_csv):###The best way should be tracking based on the tree node index; not other way can be used to fulfill all situations
 df_segments_counts2=df_site_pattern_new[["CHROM","dataset_segment_id"]].value_counts().to_frame(name="count").reset_index()
 df_segs=df_segments_counts2[df_segments_counts2["count"]>=2].reset_index(drop=True)
 df_site_pattern_new=df_site_pattern_new.copy().reset_index(drop=True)
 ###change direction of REF event first
 for index, row in df_site_pattern_new.iterrows():
  if row["desc_label"]=="REF":
   anc_cnv=row["anc_cnv"]
   desc_cnv=row["desc_cnv"]
   df_site_pattern_new.loc[index,"anc_cnv"]=desc_cnv
   df_site_pattern_new.loc[index,"desc_cnv"]=anc_cnv
 ###
 df_site_pattern_sorted=df_site_pattern_new.copy()
 df_site_pattern_sorted["direction"]=df_site_pattern_sorted["anc_cnv"].astype(str)+"->"+df_site_pattern_sorted["desc_cnv"].astype(str)
 ###also contain one change
 df_segs_directions=df_site_pattern_sorted[["CHROM","dataset_segment_id","direction"]].value_counts().reset_index().rename(columns={0:"count"})
 ###
 df_segment_record=pd.DataFrame(columns=["CHROM","segments","nGain","nLoss","direction"])
 d_direction={}#restore 3->2
 l_new_types=[]#restore 24G1BM
 for index, row in df_site_pattern_sorted.iterrows():
  chrom=row["CHROM"]
  desc_label=row["desc_label"]
  seg=row["dataset_segment_id"]
  seg_mark=str(chrom)+"_"+str(seg)
  type_i=row["type"]
  anc_cnv=row["anc_cnv"]
  desc_cnv=row["desc_cnv"]
  direction_i=f"{anc_cnv}->{desc_cnv}"
  ###how many duplications of this direction
  direction_count=df_segs_directions[(df_segs_directions["CHROM"]==chrom)&(df_segs_directions["dataset_segment_id"]==seg)&(df_segs_directions["direction"]==direction_i)]["count"].values[0]
  if seg_mark not in d_direction:
   d_direction[seg_mark]=[]
  d_direction[seg_mark].append(direction_i)##record the cnv change like 2->3 (gain)
  if type_i=="GAIN":#
   mark="G"
   nGain=1
   nLoss=0
  else:
   mark="L"
   nGain=0
   nLoss=1
  #print('df_segment_record')
  #print(df_segment_record)
  ###Run changing type
  if seg_mark not in df_segment_record["segments"].to_list():#the first Gain or Loss
   df_segment_record.loc[len(df_segment_record)]=pd.Series([chrom,seg_mark,nGain,nLoss,direction_i],index=df_segment_record.columns)
   if direction_count==1:
    new_type=str(seg)+mark+"1"
   else:
    new_type=str(seg)+mark+"1R1"
  else:#
   '''
   gain0=df_segment_record[df_segment_record["segments"]==seg]["nGain"].values[0]
   loss0=df_segment_record[df_segment_record["segments"]==seg]["nLoss"].values[0]
   df_segment_record.loc[df_segment_record["segments"]==seg,"nGain"]=gain0+nGain
   df_segment_record.loc[df_segment_record["segments"]==seg,"nLoss"]=loss0+nLoss
   '''
   if nGain==1:#
    if direction_count==1:#first Gain
     #if desc_label!="REF":
     step=desc_cnv-2
     new_type=str(seg)+"G"+str(step)
    else:#R1...
     L_direct=d_direction[seg_mark]
     repeat_times=L_direct.count(direction_i)
     step=desc_cnv-2
     new_type=str(seg)+"G"+str(step)+"R"+str(repeat_times)
   if nLoss==1:
    if direction_count==1:
     if anc_cnv<=2:
      step=2-desc_cnv
      new_type=str(seg)+"L"+str(step)
     else:###back mutation 3->2 or 4->3
      step=anc_cnv-2
      new_type=str(seg)+"G"+str(step)+"BM"
    else:#repeat R1 R2
     L_direct=d_direction[seg_mark]
     repeat_times=L_direct.count(direction_i)
     if anc_cnv<=2:#2->1 1->0
      step=2-desc_cnv
      new_type=str(seg)+"L"+str(step)+"R"+str(repeat_times)
     else:###BM BM1 BM2 3->2 two times
      step=anc_cnv-2
      new_type=str(seg)+"G"+str(step)+"BM"+str(repeat_times)
  #print(seg_mark,direction_i,new_type,direction_count,nGain,nLoss)
  l_new_types.append(new_type) 
 df_site_pattern_sorted=df_site_pattern_sorted.copy()
 df_site_pattern_sorted["new_type"]=l_new_types 
 df_site_pattern_sorted.to_csv(opt_csv,sep=",",header=True,index=False)
 return df_site_pattern_sorted

def compare_nodes_anc_desc_status(anc_index_list):#store every [anc_i,desc_i] pair
 df_edge_tmp=df_edge[(df_edge["anc"].isin(anc_index_list))|(df_edge["desc"].isin(anc_index_list))]
  

###########################Filter out (Or correct) wrong events 
def check_single_sample_events(df_site_pattern_new):#check and find in df_segment_samples_residue
 df_segs=df_segments_counts[df_segments_counts["count"]>=2].reset_index(drop=True)
 l_df=[]
 for index,row in df_segs.iterrows():
  chrom=str(row["CHROM"])
  #seg=str(row["dataset_segment_id"])
  seg=row["dataset_segment_id"]
  df_tmp=df_site_pattern_new[(df_site_pattern_new["CHROM"]==chrom)&(df_site_pattern_new["dataset_segment_id"]==seg)&(df_site_pattern_new["desc_label"]!="REF")].reset_index(drop=True)
  l_event_rec=[]
  l_residur_i=[]
  l_modify_desc_cnv=[]
  l_modify_desc_cnv_mark=[]
  for index1,row1 in df_tmp.iterrows():
   anc_index=row1["anc_index"]
   desc_index=row1["desc_index"]
   anc_cnv=row1["anc_cnv"]
   desc_cnv=row1["desc_cnv"]
   desc_sample=row1["desc_label"]
   if desc_sample!="N":
    #print(chrom,seg,desc_sample)
    residue_item=df_segment_samples_residue[(df_segment_samples_residue["CHROM"]==chrom)&(df_segment_samples_residue["dataset_segment_id"]==seg)][desc_sample].values[0]#median|round(median)|positive|negative
    l_residur_i.append(residue_item)
    residue_item_list=residue_item.split("|")
    residue_item_list1=list(map(float,residue_item_list))
    median=residue_item_list1[0]
    round_median=residue_item_list1[1]
    bottom_integer=np.floor(median)
    t1=round_median-0.25
    t2=round_median+0.25
    positive=residue_item_list1[2]
    negative=residue_item_list1[3]
    high_prop=positive/(positive+negative)
    low_prop=negative/(positive+negative)
    if median>=t1 and median<=t2:#in range
     event_rec="correct"
     final_mark="YES"
     new_desc_cnv=desc_cnv
    else:
     if bottom_integer==round_median:#correctly chosed
      event_rec="correct"
      final_mark="YES"
      new_desc_cnv=desc_cnv
     else:#modified
      new_desc_cnv=bottom_integer
      if bottom_integer==anc_cnv:#no change
       final_mark="NO"
       event_rec="wrong"
      else:#modified
       final_mark="Modified"###maintain this event!!!!!
       event_rec="corrected"
    l_event_rec.append(event_rec)
    l_modify_desc_cnv.append(new_desc_cnv)
    l_modify_desc_cnv_mark.append(final_mark)
    #print(seg,desc_sample,t1,t2,median,event_rec)
   else:##not only one sample one event
    l_residur_i.append("N")
    event_rec="correct"
    l_event_rec.append(event_rec)
    l_modify_desc_cnv.append(desc_cnv)
    l_modify_desc_cnv_mark.append("YES")
  df_tmp=df_tmp.copy()
  #print("Seg:",seg)
  #print(df_tmp)
  #print(l_event_rec)
  #print(len(df_tmp),len(l_event_rec))
  df_tmp["event_record"]=l_event_rec
  df_tmp["residue"]=l_residur_i
  df_tmp["modified_desc"]=l_modify_desc_cnv
  df_tmp["modified_desc_syb"]=l_modify_desc_cnv_mark
  l_df.append(df_tmp)
 dfs=pd.concat(l_df)
 dfs.to_csv("All_events_at_least_2.csv",sep=",",header=True,index=False)
 return dfs

def remove_wrong_events(df_site_pattern_new,df_at_least_2_correct_residue):
 df_merge=pd.merge(df_site_pattern_new,df_at_least_2_correct_residue,on=list(df_site_pattern_new.columns),how="left")
 df_merge=df_merge.copy()
 df_merge["modified_desc"]=df_merge["modified_desc"].fillna(df_merge["desc_cnv"])
 df_merge["desc_cnv"]=df_merge["modified_desc"].astype(int)
 df_merge=df_merge.fillna("N")
 df_merge["diff"]=df_merge["desc_cnv"]-df_merge["anc_cnv"]
 df_merge["abs_diff"]= df_merge["diff"].abs()
 df_merge.to_csv("All39chrs_changes_info_decomposed_Corrected_mark.csv",sep=",",header=True,index=False)
 df_final=df_merge[df_merge["event_record"]!="wrong"]
 df_final.to_csv("All39chrs_changes_info_decomposed_Corrected.csv",sep=",",header=True,index=False)
 return df_final

def get_residue(df_All_changes_info_decomposed):###residue matrix for all samples in all segments
 l_segment=[]
 df_store_residue=pd.DataFrame(columns=["CHROM","dataset_segment_id","START","END"]+all_ordered_samples)
 for index,row in df_original[["CHROM","dataset_segment_id","START","END"]].drop_duplicates().iterrows():
  chrom=row["CHROM"]
  start=row["START"]
  end=row["END"]
  segment_id=row["dataset_segment_id"]
  df_cnv=df_CNV[(df_CNV["CHROM"]==chrom)&(df_CNV["dataset_segment_id"]==segment_id)&(df_CNV["START"]>=start)&(df_CNV["END"]<=end)&(df_CNV["excluded_from_segmentation"]==False)]
  l_values=[chrom,segment_id,start,end]
  for sample_i in all_ordered_samples:
   median_value=df_cnv[sample_i].median()
   predicted_value=round(median_value)
   n_negative=len(df_cnv[df_cnv[sample_i]<predicted_value])
   n_positive=len(df_cnv[df_cnv[sample_i]>predicted_value])
   str_median="%.4f"%median_value
   item_i=f"{str_median}|{predicted_value}|{n_positive}|{n_negative}"
   l_values.append(item_i)
  df_store_residue.loc[len(df_store_residue)]=pd.Series(l_values,index=df_store_residue.columns)
 df_store_residue.to_csv("All_residue_matrix_new.csv",sep=",",header=True,index=False) 
 return df_store_residue

def pro_anc_desc(df_site_pattern_new,opt_csv):###generate events matrix based on tree structure
 df2=pd.DataFrame(columns=["CHROM","dataset_segment_id","START","END","type","new_type"]+all_ordered_samples)
 df2=df2.copy()
 for index,row in df_site_pattern_new.iterrows():
  anc_index=row["anc_index"]
  desc_index=row["desc_index"]
  desc_label=row["desc_label"]
  abs_diff=row["abs_diff"]
  type_i0=row["type"]
  type_i1=row["new_type"]
  segment=row["dataset_segment_id"]
  chrom=row["CHROM"]
  df_chr_segment=df_original[(df_original["CHROM"]==str(chrom))&(df_original["dataset_segment_id"]==int(segment))]
  seg_start=df_chr_segment["START"].values[0]
  seg_end = df_chr_segment["END"].values[0]
  row_values=[chrom,segment,seg_start,seg_end,type_i0,type_i1]###
  if desc_label=="REF":
   all_samples=df_anc_desc_nodes[df_anc_desc_nodes["anc_node"]==anc_index]["tip_labels"].values[0]
   all_samples_list=[i for i in all_samples.split("|") if i!="REF"]
   event_0_1_new1=event_0_1.copy()
   for i in all_samples_list:
    i_index=d_all_ordered_samples[i]
    event_0_1_new1[i_index]=1
   row_values1=row_values+event_0_1_new1
  else:##anc not the root
   event_0_1_new2=event_0_1.copy()
   if desc_label !="N" :##leaf node
    print("leaf node:",desc_label)
    i_index=d_all_ordered_samples[desc_label]
    event_0_1_new2[i_index]=1
   else:##parent node
    all_samples=df_anc_desc_nodes[df_anc_desc_nodes["anc_node"]==desc_index]["tip_labels"].values[0]
    all_samples_list=[i for i in all_samples.split("|")]
    for i in all_samples_list:
     i_index=d_all_ordered_samples[i]
     event_0_1_new2[i_index]=1
    n_1=len([i for i in event_0_1_new2 if i==1])
    print("parent node:",desc_index,len(all_samples_list),n_1)
    print(all_samples)
   row_values1=row_values+event_0_1_new2
  for step_i in range(abs_diff):##step higher than 1 should be duplicated 2 times
   #print("###Step higher than 2:",chrom,segment,anc_index,desc_index)
   df2.loc[len(df2)]=pd.Series(row_values1,index=df2.columns)
 #reorder the events according the segment_id
 df2["index"]=df2.index
 df2=df2.sort_values(by=["dataset_segment_id","index"],ascending=[True,True])
 df2.drop("index",axis=1,inplace=True)
 df2.to_csv(opt_csv,sep=",",header=True,index=False)
 return df2

def calculate_G_L_depth(df):
 df=df[["CHROM","dataset_segment_id","START","END","type"]]
 df1=df[["CHROM","dataset_segment_id","START","END","type"]].value_counts().reset_index(name="count")
 df2=pd.pivot_table(df1,values="count",index=["CHROM","dataset_segment_id","START","END"],columns=["type"]).fillna(0).reset_index()
 df2.to_csv("All_Gain_Loss_counts.csv",sep=",",header=True,index=False)

def generate_chrom_map_plot_input(df,opt_csv):
 df=df.copy()
 df=df.reset_index(drop=True)
 df["WIDTH"]=df["END"]-df["START"]+1
 df=df.rename(columns={"type":"EVENT TYPE"})
 df=df.fillna("N")
 df1=df.copy()
 chrom_order=list(map(str,list(range(1,39))))+["X"]
 df1["CHROM"]=pd.Categorical(df1["CHROM"],categories=chrom_order)
 df1=df1.sort_values(["CHROM","START","END"]).reset_index(drop=True)
 df1["CNV EVENT ID"]=df1.index+1
 df_end1=df1[["CNV EVENT ID","CHROM","START","END","WIDTH","EVENT TYPE"]]
 df_end1.to_csv(opt_csv,sep="\t",header=True,index=False)

def pro_chrN(chrN):
 chrN_mark="chr"+chrN#chrN should be string
 opt_csv3=chrN_mark+"_all_changes_info_decomposed_with_new_type.csv"
 df_site_pattern_new00=modify_type(df_site_pattern_new,opt_csv3)
 df_chr=df_original[df_original["CHROM"]==chrom]
 opt_csv4="events_matrix.csv"
 df_event_matrix=pro_anc_desc(df_site_pattern_new00,opt_events_matrix_csv)
 generate_chrom_map_plot_input(df_event_matrix,"events_with_wrong_chrom_map.csv")
 df_event_matrix=df_event_matrix.copy()
 #remove wrong events from change
 df_site_pattern_new1=remove_wrong_events(df_site_pattern_new,df_at_least_2_correct_residue)
 #add BM G1R1 G1R2
 df_event_matrix2=pro_anc_desc(df_site_pattern_new2,"events_matrix_correct.csv")
 #generate chromosome map plot

if __name__=="__main__":
 #edge info with node index from the tree
 df_edge=pd.read_csv("edge_table.csv",sep=",",skiprows=1,names=["anc","desc"])
 all_ordered_samples=['401T', '1382_399T', '982T', '1322_399T', '809T', '1381T', '1380T', '24T', '1532T', '2T', '1208T', '560T', '464T', '349T2', '645T', '459T', '550T', '1281T', '556T', '559T', '1315T', '1247T', '2070T1', '1940T1', '423T', '1210T', '666T', '468T', '609T', '608T', '131T', '126T', '79Ta', '773T1', '439T', '410T', '365T', '355T', '627T1', '683T', '652T', '855T', '851T', '341Ta', '335Ta']
 d_all_ordered_samples={i:n for n,i in enumerate(all_ordered_samples)}
 event_0_1=[0]*len(all_ordered_samples)
 ###anc with all desc nodes and labels if desc is tip node
 df_anc_desc_nodes=pd.read_csv("anc_desc_nodes_v1.csv",sep=",",header=0)
 df_original=pd.read_csv("All_totalCN.csv",sep=",",header=0,dtype={"CHROM":str})
 df_CNV=pd.read_csv("All_cns_raw.csv",sep=",",header=0,dtype={"CHROM":str})
 ###residue matrix##can be done for ALL chroms
 df_All_changes_info=pd.read_csv("All39chrs_changes_info.csv",sep=",",header=0)
 df_All_changes_info.drop("new_type",axis=1,inplace=True)
 df_All_changes_info_decomposed=df_All_changes_info.assign(segments=df_All_changes_info["segments"].str.split("|")).explode("segments").reset_index(drop=True)
 df_All_changes_info_decomposed=df_All_changes_info_decomposed.rename(columns={"segments":"dataset_segment_id"})
 df_All_changes_info_decomposed["dataset_segment_id"]=df_All_changes_info_decomposed["dataset_segment_id"].astype(int)
 df_segments_counts=df_All_changes_info_decomposed[["CHROM","dataset_segment_id"]].value_counts().to_frame(name="count").reset_index()
 df_All_changes_info_decomposed.to_csv("All39chrs_changes_info_decomposed.csv",sep=",",header=True,index=False)
 ##########Get segements;Should use new function to save time
 #df_segment_samples_residue=get_residue(df_All_changes_info_decomposed)
 df_segment_samples_residue=pd.read_csv("All_residue_matrix_new.csv",sep=",",header=0)
 df_segment_samples_residue=df_segment_samples_residue.copy()
 #df_segment_samples_residue["dataset_segment_id"]=df_segment_samples_residue["dataset_segment_id"].astype(str)
 #df_All_changes_info_decomposed_new_type=modify_type(df_All_changes_info_decomposed,"All39chrs_changes_info_decomposed_NewType.csv")
 print("###Correct Events:")
 ########
 df_at_least_2_correct_residue=check_single_sample_events(df_All_changes_info_decomposed)
 #df_at_least_2_correct_residue=pd.read_csv("All_events_at_least_2.csv",sep=",",header=0)
 #df_at_least_2_correct_residue["dataset_segment_id"]=df_at_least_2_correct_residue["dataset_segment_id"].astype(str)
 print("###Remove wrong events:")
 df_site_pattern_new1=remove_wrong_events(df_All_changes_info_decomposed,df_at_least_2_correct_residue)
 print("#nEvents:",len(df_site_pattern_new1))
 ###############Modify type and generate event matrix
 print("###Modify Type:")
 df_site_pattern_new2=modify_type(df_site_pattern_new1,"All_changes_new_type.csv")
 print("###Generate Events Matrix:")
 df_event_matrix2=pro_anc_desc(df_site_pattern_new2,"All_events_matrix_correct.csv")
 calculate_G_L_depth(df_event_matrix2)
 print("###Generate chromosome map plot")
 generate_chrom_map_plot_input(df_event_matrix2,"All_events_without_wrong_chrom_map.csv")

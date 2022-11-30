library(data.table)

#args=commandArgs(trailingOnly=TRUE)

filter_2T<-function(df_cnv0,df_raw0,chrN){
print(chrN)
chrN_mark<-paste("chr",chrN,sep="")

df_cnv<-df_cnv0[CHROM==chrN]
df_cnv1<-copy(df_cnv)
setkey(df_cnv,dataset_segment_id)
diff_adjacent<-df_cnv[abs(`2T`-`1208T`)==1,c("dataset_segment_id","2T","1208T")]

df_raw1<-df_raw0[CHROM==chrN]
df_raw<-df_raw0[excluded_from_segmentation=="FALSE" & CHROM=='35']
raw_segment_median<- df_raw[,lapply(.SD,median),.SDcols=c("2T","1208T"),by=dataset_segment_id]

df_compare<-merge(diff_adjacent,raw_segment_median,by="dataset_segment_id")
df_compare[,big_round:=ifelse(`2T.y`<round(`2T.y`),1,0)]
df_compare[,in_range:=ifelse(((big_round==1)&(`2T.x`-`2T.y`>0.25))|((big_round==0)&(`2T.y`-`2T.x`>0.25)),1,0)]

df_modify<-df_compare[which(in_range==1),]

for (i in unique(df_modify[,dataset_segment_id])){
raw_2T<-df_cnv[which(dataset_segment_id==i),"2T"]
raw_1208T<-df_cnv[which(dataset_segment_id==i),"1208T"]
df_cnv[which(dataset_segment_id==i),"2T"]<-df_cnv[which(dataset_segment_id==i),"1208T"]
new_2T<-df_cnv[which(dataset_segment_id==i),"2T"]
cat(chrN_mark,i,raw_2T$`2T`,raw_2T$`1208`,new_2T$`2T`,"\n")
}

new_order=c("CHROM","START","END","dataset_segment_id",'401T', '1382_399T', '982T', '1322_399T', '809T', '1381T', '1380T', '24T', '1532T', '2T', '1208T', '560T', '464T', '349T2', '645T', '459T', '550T', '1281T', '556T', '559T', '1315T', '1247T', '2070T1', '1940T1', '423T', '1210T', '666T', '468T', '609T', '608T', '131T', '126T', '79Ta', '773T1', '439T', '410T', '365T', '355T', '627T1', '683T', '652T', '855T', '851T', '341Ta', '335Ta')
setcolorder(df_cnv,new_order)
setcolorder(df_cnv1,new_order)

new_CNV_csv<-paste(chrN_mark,"totalCN.csv",sep="_")
old_CNV_csv<-paste(chrN_mark,"totalCN_old.csv",sep="_")
print(c("Nsegments:",nrow(df_cnv)))
#fwrite(df_cnv,new_CNV_csv)
#fwrite(df_cnv1,old_CNV_csv)

new_order1=c("CHROM","START","END","excluded_from_segmentation","dataset_segment_id",'401T', '1382_399T', '982T', '1322_399T', '809T', '1381T', '1380T', '24T', '1532T', '2T', '1208T', '560T', '464T', '349T2', '645T', '459T', '550T', '1281T', '556T', '559T', '1315T', '1247T', '2070T1', '1940T1', '423T', '1210T', '666T', '468T', '609T', '608T', '131T', '126T', '79Ta', '773T1', '439T', '410T', '365T', '355T', '627T1', '683T', '652T', '855T', '851T', '341Ta', '335Ta')
setcolorder(df_raw1,new_order1)
raw_CNV_csv<-paste(chrN_mark,"cns_raw.csv",sep="_")
#fwrite(df_raw1,raw_CNV_csv)
my_list<-list("raw_matrix"=df_raw1,"new_matrix"=df_cnv)
return(my_list)
}


df_cnv0<-fread("copynumber_matrix_2022-07-01.csv")
df_raw0<-fread("copynumber_rawdata_2022-07-01.csv.gz")
print("Nbins of tumor summary CNV and raw CN:")
print(c(nrow(df_cnv0),nrow(df_raw0)))

empty_list1<-list()
empty_list2<-list()
all_chroms<-unique(df_cnv0[,.(CHROM)])
n_chroms<-nrow(all_chroms)
for (i in 1:n_chroms){
chrN<-all_chroms[i]
my_list<-filter_2T(df_cnv0,df_raw0,chrN)
raw_matrix<-my_list$raw_matrix
new_matrix<-my_list$new_matrix
empty_list1<-append(empty_list1,list(raw_matrix))
empty_list2<-append(empty_list2,list(new_matrix))
}

df_summary_CNV<-rbindlist(empty_list2)
df_raw_CN<-rbindlist(empty_list1)
fwrite(df_summary_CNV,"All_totalCN.csv")
fwrite(df_raw_CN,"All_cns_raw.csv")


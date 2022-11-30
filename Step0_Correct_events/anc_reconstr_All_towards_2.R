library(phangorn)
library(ape)
library(data.table)
library(RColorBrewer)
library(phytools)

setwd("/Users/ql4/Documents/work2022/work05/work_0705_chr1/")

tree <- read.tree("/Users/ql4/Documents/work2022/work05/work_0526_CTVT_CNV/ctvt_c2t_cpg.fa.treefile")
fwrite(data.table(tree$edge),"edge_table.csv")
tree <- midpoint(tree, "REF")
# tree <- drop.tip(tree, "REF")
tree <- ladderize(tree)
plot(tree,font=1,adj=1,show.node.label = T,cex=1,align.tip.label = T)
nodelabels(frame = "none")
tiplabels(frame = "none")

new_order=c("CHROM","START","END","dataset_segment_id",
            '401T', '1382_399T', '982T', '1322_399T', '809T', 
            '1381T', '1380T', '24T', '1532T', '2T', '1208T', 
            '560T', '464T', '349T2', '645T', '459T', '550T', 
            '1281T', '556T', '559T', '1315T', '1247T', '2070T1', 
            '1940T1', '423T', '1210T', '666T', '468T', '609T', 
            '608T', '131T', '126T', '79Ta', '773T1', '439T', 
            '410T', '365T', '355T', '627T1', '683T', '652T', 
            '855T', '851T', '341Ta', '335Ta')

args<-commandArgs(trailingOnly = TRUE)

data <- fread("All_totalCN.csv")
data[, REF := 2]
n <- colnames(data)[5:50]
data[, (n) := lapply(.SD, function(col) pmax(0, pmin(30, col))), .SDcols = n]


###########Function to get all changes info
get_chrN_changes_info<-function(data,chrN){
data1<-data[CHROM==chrN]
#as.character(args[1])
print("CHROM:")
print(chrN)
chrN_mark<-paste("chr",chrN,sep="")

rownum_to_sitepattern <- function(n, anc_data) {
    attr(anc_data, "index")[n]
}

sitepattern_to_rownums <- function(p, anc_data) {
    which(attr(anc_data, "index") == p)
}

fix_vector2<-function(v){
  v_cnv<-sum(v*(seq_along(v)-1))
  print(v_cnv)
  max_05<-max(which(v<1 & v>0))
  min_05<-min(which(v<1 & v>0))
  cnv<-sum(v*(seq_along(v)-1))
  v[v<1]<-0
  v[floor(cnv)+1]<-1
  new_cnv<-sum(v*(seq_along(v)-1))
  print(new_cnv)
  return(v)
}

fix_vector<-function(v){
  v_cnv<-sum(v*(seq_along(v)-1))
  print(v_cnv)
  max_05<-max(which(v<1 & v>0))
  min_05<-min(which(v<1 & v>0))
  cnv<-sum(v*(seq_along(v)-1))
  if (cnv<2){
    v[v<1]<-0
    v[ceiling(cnv)+1]<-1
  }else{
    v[v<1]<-0
    v[floor(cnv)+1]<-1
  }
  new_cnv<-sum(v*(seq_along(v)-1))
  print(new_cnv)
  return(v)
}

count_changes <- function(tree, anc_data, site_pattern) {
    changes <- 0
    site_pattern_i_segments<-paste(which(attr(anc_data,"index")==site_pattern),collapse = "|")
    record_df<-setNames(data.table(matrix(nrow=0,ncol=11)),c("CHROM","site_pattern","anc_index","desc_index","desc_label","anc_cnv","desc_cnv","diff","type","new_type","segments"))
    for (index in rev(ape::postorder(tree))) {
        edge <- tree$edge[index, ]
        ancestor <- edge[1]
        descendent <- edge[2]
        
        reconstr_anc_name <- attr(anc_data, "names")[ancestor]
        reconstr_desc_name <- attr(anc_data, "names")[descendent]
        
        ancestral_data0 <- anc_data[[ancestor]][site_pattern, ]
        descendent_data0 <- anc_data[[descendent]][site_pattern, ]
        ###Fix the ambiguity to max copy number;0,0.5,0.5 -->0,0,1
        if (length(ancestral_data0[ancestral_data0>0 & ancestral_data0<1])>0){
          anc_data[[ancestor]][site_pattern, ]<-fix_vector(anc_data[[ancestor]][site_pattern, ])
          print(list(site_pattern,ancestor,descendent))
          print(ancestral_data0)
          print(anc_data[[ancestor]][site_pattern, ])
        }
        if (length(descendent_data0[descendent_data0>0 & descendent_data0<1])>0){
          anc_data[[descendent]][site_pattern, ]<-fix_vector(anc_data[[descendent]][site_pattern, ])
        }
        
        ancestral_data <- anc_data[[ancestor]][site_pattern, ]
        descendent_data <- anc_data[[descendent]][site_pattern, ]
        
        if (descendent<=46){
          desc_label<-attr(anc_data,"names")[descendent]
        }else{
          desc_label<-"N"
        }
        anc_CNV<-sum(ancestral_data * (seq_along(ancestral_data)-1))
        desc_CNV<-sum(descendent_data * (seq_along(descendent_data)-1))
        diff<-desc_CNV-anc_CNV
        if (diff>0){
          if(desc_label=="REF"){
            type<-"LOSS"}
          else{
            type<-"GAIN"
          }
        }else{
          if(desc_label=="REF"){
            type<-"GAIN"}
          else{
            type<-"LOSS"
          }
        }
        abs_diff <- abs( anc_CNV-desc_CNV )
        if (abs_diff>0){
          all_desc<-getDescendants(tree,descendent)
          tip_desc<-all_desc[all_desc<47]
          tip_nodes_labels<-tree$tip.label[tip_desc]
          desc_nodes_i<-paste(all_desc,collapse = "|")
          tip_nodes_i<-paste(tip_desc,collapse = "|")
          tip_nodes_labels_i<-paste(tip_nodes_labels,collapse = "|")
          new_type<-""
          #print(c(ancestor,desc_nodes_i))
          #print(list(site_pattern,ancestor,descendent,anc_CNV,desc_CNV,diff))
          record_df<-rbind(record_df,list(chrN,site_pattern,ancestor,descendent,desc_label,anc_CNV,desc_CNV,diff,type,new_type,site_pattern_i_segments))
        }
        changes <- changes + abs_diff
    }
    #print(record_df)
    print(anc_data)
    my_list<-list("changes"=changes,"record_df"=record_df,"anc_data_new"=anc_data,"segments"=site_pattern_i_segments)
    return(my_list)
}

m = as.data.frame(data[CHROM == chrN, 5:ncol(data)])
m[m<0] <- 0
p <- phyDat(m, type = "USER", levels = unique(sort(as.matrix(m))))
anc_mpr <- ancestral.pars(tree, p, type = "MPR")

n_colours<-max(data1[,c(5:50)])+1
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours<-getPalette(n_colours)

par(mfrow = c(1, 1))
plotAnc(tree, anc_mpr, 30, site.pattern = TRUE, cex.pie = 0.5)
plotAnc(tree, anc_mpr, 30, site.pattern = TRUE, cex.pie = 0.5, col = colours)

###use the ambiguity corrected changes to do plot
#all_changes <- sapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes0(tree, anc_mpr, i))
all_changes <- sapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes(tree, anc_mpr, i)$changes)
#anc_data1<-lapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes(tree, anc_mpr, i)$anc_data_new)

##############Get diff anc desc and all changes info(the order is important)
empty_list<-list()
df_undecomposed_change<-setNames(data.table(matrix(nrow=0,ncol=4)),c("CHROM","site_pattern","changes","segments"))
for (site_pattern_i in seq_len(max(attr(anc_mpr, "index")))){
  my_list_i<-count_changes(tree, anc_mpr, site_pattern_i)
  changes_i<-my_list_i$changes
  print(c("site_pattern","change",site_pattern_i,changes_i))
  df_changes_record<-my_list_i$record_df
  segments<-my_list_i$segments
  print(segments)
  empty_list<-append(empty_list,list(df_changes_record))
  df_undecomposed_change<-rbind(df_undecomposed_change,list(chrN,site_pattern_i,changes_i,segments))
}

#opt_csv0<-paste(chrN_mark,"site_patterns_changes_segments.csv",sep="_")
#fwrite(df_undecomposed_change,opt_csv0)
df_all_record<-rbindlist(empty_list)

#opt_csv01<-paste(chrN_mark,"all_changes_info.csv",sep="_")
#fwrite(df_all_record,opt_csv01)

###############

###07.01 plot tree use new CNV in which ambiguity has been removed
###Plot tree for each chromosome
###After running in python to get "All39chrs_changes_info_decomposed.csv" file
#
tree_pdf<-paste(chrN_mark,"trees.pdf",sep="_")
pdf(tree_pdf, height = 11.75, width = 8.25)
par(mfrow = c(1, 1), mar = c(4,3,3,2), oma = c(0,0,4,0))
df_changes_info<-fread("All39chrs_changes_info_decomposed.csv")
chr_x<-32
#16
df_changes_info_x<-df_changes_info[CHROM==chr_x]
par(mfrow = c(1, 1))
#seq_along(attr(anc_mpr, "index"))
for (site in seq_along(attr(anc_mpr, "index"))) {
  print(paste("site",site))
  location <- data[CHROM == chrN][site, sprintf("%s:%d-%d", CHROM, START, END)]
  site_pattern_i <- attr(anc_mpr, "index")[site]
  changes <- all_changes[site_pattern_i]
  df_changes_pattern_i<-df_changes_info_x[site_pattern==site_pattern_i]
  df_mark<-df_changes_pattern_i[!duplicated(df_changes_pattern_i),]
  plotAnc(tree, anc_mpr, site_pattern_i, site.pattern = TRUE, cex.pie = 0.5)
  #, col = colours
  for(i in 1:91){
    cnv_i_v<-anc_mpr[[i]][site_pattern_i,]
    cnv_i<-sum(cnv_i_v*(seq_along(cnv_i_v)-1))
    #bg_color<-colours[cnv_i+1]
    ####text all nodes with original Copy Number
    nodelabels(cnv_i,i,frame = "none",adj = 0.5,font=1,cex=0.8)
  }
  nChanges<-nrow(df_mark)
  print(paste("nChanges:",nChanges))
  if(nChanges>0){
    for (i in 1:nrow(df_mark)){
      nodelabels(df_mark[i,anc_cnv],df_mark[i,anc_index],frame = "c",adj = 0.5,font=1,cex=0.8,bg="red")
      nodelabels(df_mark[i,desc_cnv],df_mark[i,desc_index],frame = "c",adj = 0.5,font=1,cex=0.8,bg="red")
    }
  }
  #plotAnc(tree, anc_data1[[site_pattern]], site_pattern, site.pattern = TRUE, cex.pie = 0.5, col = colours)
  mtext(sprintf("Segment %d = %s; Changes = %.2f", site, location, changes), outer = TRUE)
}
dev.off()
#

########segment with number of changes
df_segment_change<-setNames(data.table(matrix(nrow=0,ncol=5)),c("CHROM","segment_ID","START","END","change"))
total_segments<-length(attr(anc_mpr,"index"))
for (i in 1:total_segments){
  start<-data1[dataset_segment_id==i,START]
  end<-data1[dataset_segment_id==i,END]
  change_index=attributes(anc_mpr)$index[i]
  print(c(i,change_index,all_changes[change_index]))
  df_segment_change<-rbind(df_segment_change,list(chrN,i,start,end,all_changes[change_index]))
}
#opt_csv1<-paste(chrN_mark,"segments_changes.csv",sep="_")
#fwrite(df_segment_change,opt_csv1)
endlist<-list("df_undecomposed_change"=df_undecomposed_change,
              "df_segment_change"=df_segment_change,
              "df_all_record"=df_all_record)
return(endlist)
}

##########
empty_list11<-list()
empty_list12<-list()
empty_list13<-list()
all_chroms<-unique(data[,.(CHROM)])
n_chroms<-nrow(all_chroms)
for (i in 1:n_chroms){
  chrN<-all_chroms[i]$CHROM
  final_list<-get_chrN_changes_info(data,chrN)
  df_undecomposed_change1<-final_list$df_undecomposed_change
  df_segment_change1<-final_list$df_segment_change
  df_changes_record1<-final_list$df_all_record
  empty_list11<-append(empty_list11,list(df_changes_record1))
  empty_list12<-append(empty_list12,list(df_segment_change1))
  empty_list13<-append(empty_list13,list(df_undecomposed_change1))
}

df_all_records00<-rbindlist(empty_list11)
opt_csv09<-"All39chrs_changes_info.csv"
fwrite(df_all_records00,opt_csv09)

df_all_records01<-rbindlist(empty_list12)
opt_csv10<-"All39chrs_segments_changes.csv"
fwrite(df_all_records01,opt_csv10)

df_all_records02<-rbindlist(empty_list13)
opt_csv11<-"All39chrs_site_pattern_changes_segments.csv"
fwrite(df_all_records02,opt_csv11)

rm(list=ls())
####2022.06.28
####generate anc node with desc
####node with parent
df_final<-setNames(data.table(matrix(nrow=0,ncol = 5)),c("anc_node","desc_nodes","tip_nodes","tip_labels","CNV_status"))

for(anc_node_i in attr(anc_mpr,"names")[47:91]){
  
  all_desc<-getDescendants(tree,anc_node_i)
  tip_desc<-all_desc[all_desc<47]
  tip_nodes_labels<-tree$tip.label[tip_desc]
  desc_nodes_i<-paste(all_desc,collapse = "|")
  tip_nodes_i<-paste(tip_desc,collapse = "|")
  tip_nodes_labels_i<-paste(tip_nodes_labels,collapse = "|")
  print(c(anc_node_i,desc_nodes_i))
  df_final<-rbind(df_final,list(anc_node_i,desc_nodes_i,tip_nodes_i,tip_nodes_labels_i),fill=T)
} 
#opt_csv2<-paste(chrN_mark,"anc_desc_nodes_info.csv",sep="\t")
fwrite(df_final,"anc_desc_nodes_info.csv")

df_final2<-setNames(data.table(matrix(nrow=0,ncol=3)),c("node","node_label","node_parent"))
for (node_index_i in 2:91){
  if (node_index_i<=46){
    node_label<-tree$tip.label[node_index_i] 
  }else{
    node_label<-node_index_i
  }
  if (node_index_i==47){
    parent_node_i<-"N"
  }else{
    parent_node_i<-getParent(tree,node_index_i)
  }
  df_final2<-rbind(df_final2,list(node_index_i,node_label,parent_node_i))
}
#opt_csv3<-paste(chrN_mark,"nodes_with_parent.csv",sep='\t')
fwrite(df_final2,"nodes_with_parent.csv")

####
####

# #install.packages("data.tree")
# dtree<-as.Node(tree)
# for (leaf_path_i in attr(dtree$leaves,"names")){
#   all_parent_nodes_with_leaf_name<-strsplit(leaf_path_i,"\\.")[[1]]
#   leaf_i_name<-tail(all_parent_nodes_with_leaf_name,n=1)
#   leaf_i_index<-match(leaf_i_name,attr(anc_mpr,"names"))
# }
# 
# ####2022.06.29
# strsplit(attr(dtree$leaves,"names")[1],"\\.")
# #leaf node label
# tail(strsplit(attr(dtree$leaves,"names")[1],"\\.")[[1]],n=1)


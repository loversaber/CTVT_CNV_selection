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


###########Function to plot tree after running python script
plot_chrN_tree<-function(data,chrN){
data1<-data[CHROM==chrN]
#as.character(args[1])
print("CHROM:")
chrN_mark<-paste("Chr",chrN,sep="")
print(chrN_mark)

fix_vector<-function(v){
  v_cnv<-sum(v*(seq_along(v)-1))
  print(v_cnv)
  max_05<-max(which(v<1 & v>0))
  min_05<-min(which(v<1 & v>0))
  cnv<-sum(v*(seq_along(v)-1))
  #eg. 1.5 -> 2; 3.5 -> 3
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
#plotAnc(tree, anc_mpr, 30, site.pattern = TRUE, cex.pie = 0.5, col = colours)

###use the ambiguity corrected changes to do plot
#all_changes <- sapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes0(tree, anc_mpr, i))
all_changes <- sapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes(tree, anc_mpr, i)$changes)
#anc_data1<-lapply(seq_len(max(attr(anc_mpr, "index"))), function(i) count_changes(tree, anc_mpr, i)$anc_data_new)


###07.01 plot tree use new CNV in which ambiguity has been removed
###Plot tree for each chromosome
###After running in python to get "All39chrs_changes_info_decomposed.csv" file
#
tree_pdf<-paste(chrN_mark,"trees.pdf",sep="_")
pdf(tree_pdf, height = 11.75, width = 8.25)
par(mfrow = c(1, 1), mar = c(4,3,3,2), oma = c(0,0,4,0))
df_changes_info<-fread("All39chrs_changes_info_decomposed.csv")
chr_x<-chrN
#16
df_changes_info_x<-df_changes_info[CHROM==chr_x]
par(mfrow = c(1, 1))
#seq_along(attr(anc_mpr, "index"))
for (site in seq_along(attr(anc_mpr, "index"))) {
  print(paste("Segment:",site))
  location <- data[CHROM == chrN][site, sprintf("%s:%d-%d", CHROM, START, END)]
  site_pattern_i <- attr(anc_mpr, "index")[site]
  changes <- all_changes[site_pattern_i]
  df_changes_pattern_i<-df_changes_info_x[site_pattern==site_pattern_i & dataset_segment_id==site]
  df_mark<-df_changes_pattern_i[!duplicated(df_changes_pattern_i),]
  plotAnc(tree, anc_mpr, site_pattern_i, site.pattern = TRUE, cex.pie = 0.5)
  #, col = colours
  for(i in 1:91){
    cnv_i_v<-anc_mpr[[i]][site_pattern_i,]
    cnv_i<-sum(cnv_i_v*(seq_along(cnv_i_v)-1))
    cnv_i<-round(cnv_i,2)
    #bg_color<-colours[cnv_i+1]
    ####text all nodes with original Copy Number
    if (((i %in% df_mark[,anc_index])==FALSE)&((i %in% df_mark[,desc_index])==FALSE)){
      nodelabels(cnv_i,i,frame = "none",adj = 0.5,font=1,cex=0.6)
      }
  }
  nChanges<-nrow(df_mark)
  print(paste("nChanges:",nChanges))
  if(nChanges>0){
    for (i in 1:nrow(df_mark)){
      nodelabels(df_mark[i,anc_cnv],df_mark[i,anc_index],frame = "c",adj = 0.5,font=1,cex=0.6,bg="red")
      nodelabels(df_mark[i,desc_cnv],df_mark[i,desc_index],frame = "c",adj = 0.5,font=1,cex=0.6,bg="red")
    }
  }
  #plotAnc(tree, anc_data1[[site_pattern]], site_pattern, site.pattern = TRUE, cex.pie = 0.5, col = colours)
  mtext(sprintf("Segment %d = %s; Changes = %.2f", site, location, changes), outer = TRUE)
}
dev.off()

}

##########Run for all 39 chromosomes
all_chroms<-unique(data[,.(CHROM)])
n_chroms<-nrow(all_chroms)
for (i in 13){
  chrN<-13
  final_list<-plot_chrN_tree(data,chrN)

}

rm(list=ls())

#################CHR11
chrN<-11
final_list<-plot_chrN_tree(data,chrN)

#######CHR38
chrN<-38
final_list<-plot_chrN_tree(data,chrN)

####CHR34 CTNND2
chrN<-34
final_list<-plot_chrN_tree(data,chrN)

####CHRX TKTL1
chrN<-'X'
final_list<-plot_chrN_tree(data,chrN)

####CHR26 LIF
chrN<-26 
final_list<-plot_chrN_tree(data,chrN)

####CHR33 CD47
chrN<-33
final_list<-plot_chrN_tree(data,chrN)

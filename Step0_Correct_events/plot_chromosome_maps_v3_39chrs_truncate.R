library(readxl)
library(data.table)
library(rectanglePacking)

## set input path(s)
setwd("~/Documents/work2022/work05/work_0705_chr1//")


get_dims<-function(fil,chrN){
  cnvtable<-as.data.table(read.csv(fil,sep="\t",header=T,check.names = T))
  cnvtable<-cnvtable[CHROM==as.integer(chrN)]
  cnvtable[,CHROM:=as.integer(CHROM)]
  cnvtable[,EventType := ifelse(grepl("GAIN",EVENT.TYPE),"gain","loss")]
  cnvtable[, oldID := CNV.EVENT.ID]
  cnvtable[, ID := paste("CTVT", CHROM, CNV.EVENT.ID, sep=".")]
  chrlengths0<-as.data.table(read.csv("dog_chr_length_noX.txt",sep=" ",header = T))
  chrlengths<-chrlengths0[CHROM==chrN]
  chrlengths[,CHROM:=as.integer(CHROM)]
  chrlengths[, GenomeStart := 0]
  chrlengths[, GenomeEnd := cumsum(as.numeric(LENGTH))]
  #chrlengths[, GenomeStart := 1 + c(0, shift(GenomeEnd)[2:.N])]
  chrlengths[, GenomeOffset := GenomeStart]
  
  cnvtable[chrlengths, GenomeStartPos := START + GenomeOffset, on = "CHROM"]
  cnvtable[chrlengths, GenomeEndPos := END + GenomeOffset, on = "CHROM"]
  
  cnv_table<-cnvtable
  chrom_<-as.integer(chrN)
  for (event_type_ in cnv_table[, sort(unique(EventType))]) {
    chrom_table <- cnv_table[CHROM == chrom_ & EventType == event_type_]
    rectangle_packing_cpp(chrom_table)
    cnv_table[chrom_table, layer := i.layer, on = c("CHROM", "EventType", "ID")]
  }
  ymin <- cnv_table[EventType=="loss", c(-max(layer))]
  ymax <- cnv_table[EventType=="gain", c(max(layer))]
  if(ymax==-Inf){
    ymax <- 1
  }else if(ymax > 15){
    ymax <- 16
  }else{
    ymax <- ymax
  }
  if(ymin==-Inf){
    ymin <- -1
  }else if(ymin < -15){
    ymin <- -14
  }else{
    ymin <- ymin
  }
  list(xlim = cnv_table[, c(min(GenomeStartPos), max(GenomeEndPos))],
         ylim = c(ymin, ymax))
}

plot_map<-function(fil,chrN,title_name,xmin_all,xmax_all){
ll1<-get_dims(fil,chrN)

ymin<-min(c(ll1$ylim))
ymax<-max(c(ll1$ylim))
  
cnvtable<-as.data.table(read.csv(fil,header = T,sep="\t",check.names = T))
cnvtable<-cnvtable[CHROM!="X"][,CHROM:=as.integer(CHROM)]
cnvtable<-cnvtable[CHROM==as.integer(chrN)]
chrlengths<-as.data.table(read.csv("dog_chr_length_noX.txt",sep=" ",header = T))
chrlengths[, GenomeStart := 0]
chrlengths[, GenomeEnd := cumsum(as.numeric(LENGTH))]
#chrlengths[, GenomeStart := 1 + c(0, shift(GenomeEnd)[2:.N])]
chrlengths[, GenomeOffset := GenomeStart ]

cnvtable[chrlengths, GenomeStartPos := START + GenomeOffset, on = "CHROM"]
cnvtable[chrlengths, GenomeEndPos := END + GenomeOffset, on = "CHROM"]


cnvtable[,EventType := ifelse(grepl("GAIN",EVENT.TYPE),"gain","loss")]
cnvtable[, oldID := CNV.EVENT.ID]
cnvtable[, ID := paste("CTVT", CHROM, CNV.EVENT.ID, sep=".")]

# Assign plot layers to each CNV
for (chrom_ in cnvtable[, sort(unique(CHROM))]) {
  for (event_type_ in cnvtable[, sort(unique(EventType))]) {
    chrom_table <- cnvtable[CHROM == chrom_ & EventType == event_type_]
    rectangle_packing_cpp(chrom_table)
    cnvtable[chrom_table, layer := i.layer, on = c("CHROM", "EventType", "ID")]
  }
}


get_plot_dimensions <- function(cnv_table) {
  ymin <- cnv_table[EventType=="loss", c(-max(layer))]
  ymax <- cnv_table[EventType=="gain", c(max(layer))]
  if(ymax==-Inf){
    ymax <- 1
  }else if(ymax > 15){
    ymax <- 16
  }else{
    ymax <- ymax
  }
  if(ymin==-Inf){
    ymin <- -1
  }else if(ymin < -15){
    ymin <- -14
  }else{
    ymin <- ymin
  }
  list(xlim = cnv_table[, c(min(GenomeStartPos), max(GenomeEndPos))],
       ylim = c(ymin, ymax))
}

plot_chrom_map <- function(cnv_table, chr_lengths, plot_options = NULL) {
  default_plot_options <- list(
    chromosome_height = 0.9,
    cnv_height = 0.7,
    col = "cornflowerblue",
    border = NA,
    title = "Chromosome Map",
    show_cnv_id = FALSE
  )
  
  if (!is.null(plot_options)) {
    stopifnot(is.list(plot_options))
    for (option in names(plot_options)) {
      default_plot_options[[option]] <- plot_options[[option]]
    }
  }
  
  plot_options <- default_plot_options
  
  dims <- get_plot_dimensions(cnv_table)
  if ("ylim" %in% names(plot_options)) {
    dims[["ylim"]] <- plot_options[["ylim"]]
  }
  ymax_new <- ceiling(ymax/5)*5
  ymin_new <- -ceiling(ymin/5)*5
  plot(NA, xlim = c(xmin_all,xmax_all), ylim = c(ymin_new,ymax_new), main = plot_options$title,
       ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
  
  # Useful to know the extents of the plotting window
  window_xrange <- par()$usr[1:2]
  window_yrange <- par()$usr[3:4]
  
  # Sort out the y-axis
  y_axis_limits <- round(2*dims$ylim, -1)/2
  y_axis_major_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 5)
  y_axis_minor_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 1)
  axis(2, at = y_axis_major_ticks, labels = abs(y_axis_major_ticks), las = 2)
  print(par()$tck)
  axis(2, at = y_axis_minor_ticks, labels = NA, tck=-0.01)
  
  mtext(c("Loss depth", "Gain depth"), 2, line=3, at = window_yrange/2, adj=c(0.5, 0.5))
  #mtext(chr_lengths[CHROM == "35", CHROM], 1, line = 0,
        #at = chr_lengths[CHROM == "35", (GenomeStart + GenomeEnd) / 2],
        #adj = 0.5)
  
  # Plot background rectangles for chromosomes
  rect(window_xrange[1], window_yrange[1], window_xrange[2], window_yrange[2],
       col = "grey95", border = NA, lwd = 2)
  chr_lengths[, rect(GenomeStart, window_yrange[1], GenomeEnd, window_yrange[2],
                     col = c("grey90", "grey95"), border = "NA")]
  GenomeStart1<-chr_lengths[CHROM == chrN, GenomeStart]
  GenomeEnd1<-chr_lengths[CHROM == chrN, GenomeEnd]
  polygon(c(GenomeStart1,GenomeStart1,GenomeEnd1,GenomeEnd1),c(-0.42,0.42,0.42,-0.42),col="grey",border="black")
  
  cnv_table[EventType == "gain",
            rect(xleft = GenomeStartPos,
                 ybottom = layer - plot_options$cnv_height/2,
                 xright = GenomeEndPos,
                 ytop = layer + plot_options$cnv_height/2,
                 col = plot_options$col, border = plot_options$border, lwd = 1)]
  cnv_table[EventType == "loss",
            rect(xleft = GenomeStartPos,
                 ybottom = -(layer - plot_options$cnv_height/2),
                 xright = GenomeEndPos,
                 ytop = -(layer + plot_options$cnv_height/2),
                 col = plot_options$col, border = plot_options$border, lwd = 1)]
  if (plot_options$show_cnv_id) {
    cnv_table[, text(x = (GenomeStartPos + GenomeEndPos)/2, y = ifelse(EventType=="gain", 1, -1) * (layer),
                     labels = oldID, adj = c(0.5, 0.5), cex = 0.6)]
  }
  
}

plot_chrom_map(cnvtable, chrlengths,
               plot_options = list(col = "cornflowerblue",
                                   ylim = c(ymin-1, ymax+1),
                                   title = title_name))
}

######################################Run Successfully
df_all<-fread("All_events_without_wrong_chrom_map.csv")
df_anc<-fread("Anc_chrom_map.csv")
opt_file<-"All_test.pdf"
pdf(opt_file,width=10,height=6)
par(mfrow = c(2, 1), mar = c(2,3,1,2), oma = c(0,0,2,0))
chr_i<-12
ll0<-get_dims("Anc_chrom_map.csv",chr_i)
ll1<-get_dims("All_events_without_wrong_chrom_map.csv",chr_i)

xmin_all<-min(c(ll1$xlim,ll0$xlim))
xmax_all<-max(c(ll1$xlim,ll0$xlim))

chr_i_mark<-paste("Chr",as.character(chr_i),sep="")
plot_map("Anc_chrom_map.csv",chr_i,chr_i_mark,xmin_all,xmax_all)
if(any(df_all[,.(CHROM)]==chr_i)){
  plot_map("All_events_without_wrong_chrom_map.csv",chr_i,"",xmin_all,xmax_all)
}else{
  plot(x=10,y=10,pch=".",col="white",xaxt="n",yaxt="n",axes=FALSE)
  text(x=10,y=10,"No Change for Ancestor")
}
print(chr_i)
dev.off()
##################################

df_all<-fread("All_events_without_wrong_chrom_map.csv")
df_anc<-fread("Anc_chrom_map.csv")

unique_chrs<-unique(df_all[,.(CHROM)])

opt_file<-"All_39chrs_chrom_map_with_Anc_truncated.pdf"
pdf(opt_file,width=10,height=6)
par(mfrow = c(2, 1), mar = c(2,3,1,2), oma = c(0,0,2,0))
for (chr_i in 1:38){
  ll0<-get_dims("Anc_chrom_map.csv",chr_i)
  ll1<-get_dims("All_events_without_wrong_chrom_map.csv",chr_i)
  
  xmin_all<-min(c(ll1$xlim,ll0$xlim))
  xmax_all<-max(c(ll1$xlim,ll0$xlim))
  
  chr_i_mark<-paste("Chr",as.character(chr_i),sep="")
  plot_map("Anc_chrom_map.csv",chr_i,chr_i_mark,xmin_all,xmax_all)
  if(any(df_all[,.(CHROM)]==chr_i)){
   plot_map("All_events_without_wrong_chrom_map.csv",chr_i,"",xmin_all,xmax_all)
  }else{
    plot(x=10,y=10,pch=".",col="white",xaxt="n",yaxt="n",axes=FALSE)
    text(x=10,y=10,"No Change for Ancestor")
  }
  print(chr_i)
  
}
dev.off()

rm(list=ls())

get_dims_X<-function(fil,chrN){
  cnvtable0<-as.data.table(read.csv(fil,sep="\t",header=T,check.names = T))
  cnvtable<-cnvtable0[CHROM==chrN]
  cnvtable[,EventType := ifelse(grepl("GAIN",EVENT.TYPE),"gain","loss")]
  cnvtable[, oldID := CNV.EVENT.ID]
  cnvtable[, ID := paste("CTVT", CHROM, CNV.EVENT.ID, sep=".")]
  chrlengths<-as.data.table(read.csv("dog_chr_length_X.txt",sep=" ",header = T))
  chrlengths[, GenomeStart := 0]
  chrlengths[, GenomeEnd := cumsum(as.numeric(LENGTH))]
  #chrlengths[, GenomeStart := 1 + c(0, shift(GenomeEnd)[2:.N])]
  chrlengths[, GenomeOffset := GenomeStart ]
  
  cnvtable[chrlengths, GenomeStartPos := START + GenomeOffset, on = "CHROM"]
  cnvtable[chrlengths, GenomeEndPos := END + GenomeOffset, on = "CHROM"]
  
  cnv_table<-cnvtable
  chrom_<-chrN
  for (event_type_ in cnv_table[, sort(unique(EventType))]) {
    chrom_table <- cnv_table[CHROM == chrom_ & EventType == event_type_]
    rectangle_packing_cpp(chrom_table)
    cnv_table[chrom_table, layer := i.layer, on = c("CHROM", "EventType", "ID")]
  }
  ymin <- cnv_table[EventType=="loss", c(-max(layer))]
  ymax <- cnv_table[EventType=="gain", c(max(layer))]
  list(xlim = cnv_table[, c(min(GenomeStartPos), max(GenomeEndPos))],
       ylim = c(ymin, ymax))
}

plot_map_X<-function(fil,chrN,title_name){
  ll1<-get_dims_X(fil,chrN)
  
  xmin_all<-min(c(ll1$xlim))
  xmax_all<-max(c(ll1$xlim))
  ymin<-min(c(ll1$ylim))
  ymax<-max(c(ll1$ylim))
  
  cnvtable0<-as.data.table(read.csv(fil,header = T,sep="\t",check.names = T))
  cnvtable<-cnvtable0[CHROM==chrN]
  chrlengths<-as.data.table(read.csv("dog_chr_length_X.txt",sep=" ",header = T))
  chrlengths[, GenomeStart := 0]
  chrlengths[, GenomeEnd := cumsum(as.numeric(LENGTH))]
  chrlengths[, GenomeOffset := GenomeStart]
  
  cnvtable[chrlengths, GenomeStartPos := START + GenomeOffset, on = "CHROM"]
  cnvtable[chrlengths, GenomeEndPos := END + GenomeOffset, on = "CHROM"]
  
  
  cnvtable[,EventType := ifelse(grepl("GAIN",EVENT.TYPE),"gain","loss")]
  cnvtable[, oldID := CNV.EVENT.ID]
  cnvtable[, ID := paste("CTVT", CHROM, CNV.EVENT.ID, sep=".")]
  
  # Assign plot layers to each CNV
  for (chrom_ in cnvtable[, sort(unique(CHROM))]) {
    for (event_type_ in cnvtable[, sort(unique(EventType))]) {
      chrom_table <- cnvtable[CHROM == chrom_ & EventType == event_type_]
      rectangle_packing_cpp(chrom_table)
      cnvtable[chrom_table, layer := i.layer, on = c("CHROM", "EventType", "ID")]
    }
  }
  
  
  get_plot_dimensions <- function(cnv_table) {
    ymin <- cnv_table[EventType=="loss", c(-max(layer))]
    ymax <- cnv_table[EventType=="gain", c(max(layer))]
    list(xlim = cnv_table[, c(min(GenomeStartPos), max(GenomeEndPos))],
         ylim = c(ymin, ymax))
  }
  
  plot_chrom_map <- function(cnv_table, chr_lengths, plot_options = NULL) {
    default_plot_options <- list(
      chromosome_height = 0.9,
      cnv_height = 0.7,
      col = "cornflowerblue",
      border = NA,
      title = "Chromosome Map",
      show_cnv_id = FALSE
    )
    
    if (!is.null(plot_options)) {
      stopifnot(is.list(plot_options))
      for (option in names(plot_options)) {
        default_plot_options[[option]] <- plot_options[[option]]
      }
    }
    
    plot_options <- default_plot_options
    
    dims <- get_plot_dimensions(cnv_table)
    if ("ylim" %in% names(plot_options)) {
      dims[["ylim"]] <- plot_options[["ylim"]]
    }
    plot(NA, xlim = c(xmin_all,xmax_all), ylim = c(ymin,ymax), main = plot_options$title,
         ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
    
    # Useful to know the extents of the plotting window
    window_xrange <- par()$usr[1:2]
    window_yrange <- par()$usr[3:4]
    
    # Sort out the y-axis
    y_axis_limits <- round(2*c(ymin,ymax), -1)/2
    y_axis_major_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 5)
    y_axis_minor_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 1)
    axis(2, at = y_axis_major_ticks, labels = abs(y_axis_major_ticks), las = 2)
    print(par()$tck)
    axis(2, at = y_axis_minor_ticks, labels = NA, tck=-0.01)
    
    mtext(c("Loss depth", "Gain depth"), 2, line=3, at = window_yrange/2, adj=c(0.5, 0.5))
    #mtext(chr_lengths[CHROM == "35", CHROM], 1, line = 0,
    #at = chr_lengths[CHROM == "35", (GenomeStart + GenomeEnd) / 2],
    #adj = 0.5)
    
    # Plot background rectangles for chromosomes
    rect(window_xrange[1], window_yrange[1], window_xrange[2], window_yrange[2],
         col = "grey95", border = NA, lwd = 2)
    chr_lengths[, rect(GenomeStart, window_yrange[1], GenomeEnd, window_yrange[2],
                       col = c("grey90", "grey95"), border = "NA")]
    GenomeStart1<-chr_lengths[CHROM == chrN, GenomeStart]
    GenomeEnd1<-chr_lengths[CHROM == chrN, GenomeEnd]
    polygon(c(GenomeStart1,GenomeStart1,GenomeEnd1,GenomeEnd1),c(-0.42,0.42,0.42,-0.42),col="grey",border="black")
    
    cnv_table[EventType == "gain",
              rect(xleft = GenomeStartPos,
                   ybottom = layer - plot_options$cnv_height/2,
                   xright = GenomeEndPos,
                   ytop = layer + plot_options$cnv_height/2,
                   col = plot_options$col, border = plot_options$border, lwd = 1)]
    cnv_table[EventType == "loss",
              rect(xleft = GenomeStartPos,
                   ybottom = -(layer - plot_options$cnv_height/2),
                   xright = GenomeEndPos,
                   ytop = -(layer + plot_options$cnv_height/2),
                   col = plot_options$col, border = plot_options$border, lwd = 1)]
    if (plot_options$show_cnv_id) {
      cnv_table[, text(x = (GenomeStartPos + GenomeEndPos)/2, y = ifelse(EventType=="gain", 1, -1) * (layer),
                       labels = oldID, adj = c(0.5, 0.5), cex = 0.6)]
    }
    
  }
  
  plot_chrom_map(cnvtable, chrlengths,
                 plot_options = list(col = "cornflowerblue",
                                     ylim = c(ymin-1, ymax+1),
                                     title = title_name))
}

opt_file<-"All_39chrs_chrom_map_with_Anc.pdf"
pdf(opt_file,width=10,height=6)
par(mfrow = c(2, 1), mar = c(2,3,1,2), oma = c(0,0,2,0))
chr_i<-"X"
plot_map_X("Anc_chrom_map.csv","X","ChrX")
if(any(df_all[,.(CHROM)]==chr_i)){
  plot_map_X("All_events_without_wrong_chrom_map.csv",chr_i,"")
}else{
  plot(x=10,y=10,pch=".",col="white",xaxt="n",yaxt="n",axes=FALSE)
  text(x=10,y=10,"No Change for Ancestor")
}
dev.off()
## clean up environment
rm(list=ls())

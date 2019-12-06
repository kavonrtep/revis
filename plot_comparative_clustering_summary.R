#!/usr/bin/env Rscript
library(optparse)



get_color = function(classification){
  ## 20 of unique colors, first is black
  unique_colors = c("#000000",
                    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075'
                    )
  ## rest wil be grey:
  grey_color = "#a9a9a9"
  unique_repeats = names(sort(table(classification), decreasing = TRUE))
  color_table = unique_colors[1:min(20,length(unique_repeats))]
  names(color_table) = unique_repeats[1:min(20,length(unique_repeats))]
  color = rep(grey_color, length(classification))
  names(color) = classification
  for (ac in names(color_table)){
    color[names(color) %in% ac] = color_table[ac]
  }
  return(color)
}


make_legend = function(color){
  ## simplify description:
  names(color) = gsub(".+/","",names(color))
  description = sapply(split(names(color), color), function(x) paste(unique(x), collapse=" "))
  legend_info = list(name = gsub("All", "NA", description), color = names(description))
  return(legend_info)
}

plot_rect_map = function(read_counts,cluster_annotation, output_file,Xcoef=1,Ycoef=1){
  ## read_counts : correspond to COMPARATIVE_ANALYSIS_COUNTS.csv
  ## cluster annotation : CLUSTER_TABLE.csv
  counts = read.table(read_counts,header=TRUE,as.is=TRUE)
  ## find which line is header
  header_line = grep("Cluster.*Supercluster.*Size", readLines(cluster_annotation))
  annot = read.table(cluster_annotation, sep="\t",header=TRUE,as.is=TRUE, skip = header_line - 1)
  ## remove counts which are not in annotation - only clusters in annot will be plotted!
  counts = counts[annot$Cluster,]
  N = nrow(annot)

  color_auto = get_color(annot$Automatic_classification)
  legend_info = make_legend(color_auto)
  manual = FALSE
  if (!is.null(annot$Manual_classification)){
    ## column with manual annotation exist - check if correct
    if (any(annot$Manual_classification %in% "")){
      message("manual classification is not complete, skipping")
    }else{
      color_manual = get_color(annot$Manual_classification)
      legend_info_manual = make_legend(color_manual)
      manual = TRUE
    }
  }
  M = as.matrix(counts[1:N,-(1:2)])
  rownames(M) = paste0("CL",rownames(M))
  Mn1=(M)/apply(M,1,max)
  Mn2=M/max(M)
  ord1 = hclust(dist(Mn1),method = "ward.D")$order
  ord2 = hclust(dist(t(Mn2)))$order
  # in inches
  wdth = (1.5 + N*0.03 ) * Xcoef
  hgt = (2 + ncol(M)*0.15) * Ycoef
  ptsize = round((wdth*hgt)^(1/4))*5

  pdf(output_file, width=wdth,height=hgt, pointsize = ptsize)  # was 50
  ploting_area_width = 3 + log10(N)*3
  ploting_area_sides = 1

  layout(matrix(c(4,2,3,4,1,3),ncol=3,byrow = TRUE),
         width=c(ploting_area_sides,ploting_area_width,ploting_area_sides),
         height=c(3,ncol(M)*0.5))
  par(xaxs='i', yaxs = 'i')
  par(las=2,mar=c(4,0,0,0),cex.axis=0.5)
  rectMap(Mn2[ord1,ord2],scale.by='none',col=color_auto[ord1], grid=TRUE)
  par(las=2,mar=c(1,0,1,0), mgp = c(2,1,0))
  barplot(annot$Size[ord1], col = 1)
  mtext(side = 2, "Cluster size", las = 3, line = 2, cex = 0.5)
  par(mfrow=c(1,1))
  plot.new()
  legend("topleft", col= legend_info$color, legend=legend_info$name, pch=15)
  if(manual){
    message("manual annotation - not implemented yet")
  }


  st = dev.off()
}

rectMap=function(x,scale.by='row',col=1,xlab="",ylab="",grid=TRUE,axis_pos=c(1,4),cexx=NULL,cexy=NULL){
  if (scale.by=='row'){
                                        #x=(x)/rowSums(x)
    x=(x)/apply(x,1,max)
  }
  if (scale.by=='column'){
    x=t(t(x)/apply(x,2,max))
  }
  nc=ncol(x)
  nr=nrow(x)
  coords=expand.grid(1:nr,1:nc)
  plot(coords[,1],coords[,2],type='n',axes=F,xlim=range(coords[,1])+c(-.5,.5),ylim=range(coords[,2])+c(-.5,.5),xlab=xlab,ylab=ylab)
  axis(axis_pos[1],at=1:nr,labels=rownames(x),lty=0,tick=FALSE,line=0,cex.axis=0.5/log10(nr))
  axis(axis_pos[2],at=1:nc,labels=colnames(x),lty=0,tick=FALSE,las=2,line=0 ,hadj=0, cex.axis=0.7)
  axis(2,at=1:nc,labels=colnames(x),lty=0,tick=FALSE,las=2,line=0 ,hadj=1, cex.axis=0.7)

  mtext(side = 1, "Cluster id", las=1, line = 3, cex = 0.5)
  line = 1.5 + log10(nr)
  mtext(side = 2, "Proportions of individual samples", las =0, line = line, cex = 0.5)
  s=(c(x)^(1/2))/2  # to get it proportional
  w = c(x)/2
  rect(coords[,1]-0.5,coords[,2]-s,coords[,1]+0.5,coords[,2]+s,col=col,border=NA)
  if (grid){
    abline(v=0:(nr)+.5,h=0:(nc)+.5,lty=2,col="#60606030")
  }
  box(col="#60606030",lty=2)
}

option_list <- list( 
  make_option(c("-c", "--cluster_table"), default=NA, type = "character",
              help="file from RepeatExplorer clustering - CLUSTER_TABLE.csv"),

  make_option(c("-m", "--comparative_counts"),default = NA,type = "character",
              help="file from RepeatExplorer2 output - COMPARATIVE_ANALYSIS_COUNTS.csv"),

  make_option(c("-o", "--output"), type="character",
              default="comparative_analysis_summary.pdf",
              help="File name for output figures (pdf document)")
)
opt = parse_args(OptionParser(option_list = option_list))
print(opt)
if (any(is.na(c(opt$cluster_table, opt$comparative_counts)))){
  message("\nBoth files: CLUSTER_TABLE.csv and COMPARATIVE_ANALYSIS_COUNTS.csv must be provided\n")
  q()
}


plot_rect_map(opt$comparative_counts, opt$cluster_table, opt$output)

## testing:
cluster_annotation = "example_data/example1_CLUSTER_TABLE.csv"
read_counts = "example_data/example1_COMPARATIVE_ANALYSIS_COUNTS.csv"
output_file = "figs/test.pdf"

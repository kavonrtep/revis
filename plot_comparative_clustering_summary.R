#!/usr/bin/env Rscript
library(optparse)

get_color = function(classification){
  ## 20 of unique colors, first is black
  unique_colors = c("#000000",
                    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075'
                    )[1:opt$number_of_colors]
  Ncol = length(unique_colors)
  ## rest wil be grey:
  grey_color = "#a9a9a9"
  ## unique repeats without All
  unique_repeats = names(sort(table(classification[!classification %in% "All"]), decreasing = TRUE))
  color_table = unique_colors[1:min(Ncol,length(unique_repeats))]
  names(color_table) = unique_repeats[1:min(Ncol,length(unique_repeats))]
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
  description = sapply(split(names(color), color), function(x) paste(unique(x), collapse=";"))
  description = gsub(".+;.+", "Other", description)
  description = gsub("All", "Other", description)
  if ("Other" %in% description & length(description) > 1){
    description = c(description[! description %in% "Other"], description[description %in% "Other"])
  }

  legend_info = list(name = gsub("All", "NA", description), color = names(description))

  return(legend_info)

}

plot_rect_map = function(read_counts,cluster_annotation, output_file,GS, Xcoef=1,Ycoef=1){
  ## read_counts : correspond to COMPARATIVE_ANALYSIS_COUNTS.csv
  ## cluster annotation : CLUSTER_TABLE.csv
  counts = read.table(read_counts,header=TRUE,as.is=TRUE)
  input_read_counts = unlist(read.table(read_counts, nrows = 1, comment.char = "",sep="\t")[-(1:2)])


  ## find which line is header
  header_line = grep("Cluster.*Supercluster.*Size", readLines(cluster_annotation))
  annot = read.table(cluster_annotation, sep="\t",header=TRUE,as.is=TRUE, skip = header_line - 1)
  ## remove counts which are not in annotation - only clusters in annot will be plotted!
  counts = counts[annot$Cluster,]
  N = nrow(annot)

  counts_automatic = counts
  annot_automatic = annot
  input_read_counts_automatic = input_read_counts
  # remove organelar and contamination if required make count correction
  if (opt$nuclear_only){
    exclude=grep("contamination|organelle",annot$Automatic_annotation)
    if (length(exclude)>0){
      counts_automatic = counts[-exclude, , drop=FALSE]
      annot_automatic = annot[-exclude, ,drop=FALSE]
      input_read_counts_automatic = input_read_counts - colSums(counts[exclude,-c(1:2) , drop=FALSE])
    }
  }
  color_auto = get_color(annot_automatic$Automatic_annotation)
  legend_info = make_legend(color_auto)
  params = list(Automatic_classification = list(
                  color =  color_auto,
                  legend = legend_info,
                  counts = counts_automatic,
                  annot = annot_automatic,
                  input_read_counts = input_read_counts_automatic
                )
                )

  if (!is.null(annot$Final_annotation)){
    
    ## column with manual annotation exist - check if correct
    if (any(annot$Final_annotation %in% "" | any(is.na(annot$Final_annotation)))){
      message("Final annotation is not complete, skipping")
    }else{

      counts_manual = counts
      annot_manual = annot
      input_read_counts_manual = input_read_counts
      ## correction must be done idependetly in case manual and automatic classification differ in annotation
      if (opt$nuclear_only){
        exclude=grep("contamination|organelle",annot$Final_annotation)
        if (length(exclude)>0){
          counts_manual = counts[-exclude, , drop=FALSE]
          annot_manual = annot[-exclude, ,drop=FALSE]
          input_read_counts_manual = input_read_counts - colSums(counts[exclude,-c(1:2) , drop=FALSE])
        }
      }
      ## append
      color_manual = get_color(annot_manual$Final_annotation)
      legend_info_manual = make_legend(color_manual)

      params$Manual_classification = list(
        color =  color_manual,
        legend = legend_info_manual,
        counts = counts_manual,
        annot = annot_manual,
        input_read_counts = input_read_counts_manual

      )
    }
  }

  ## set size of pdf output
  wdth = (3 + N*0.03 ) * Xcoef
  hgt = (2.2 + ncol(counts)*0.15) * Ycoef
  ptsize = round((wdth*hgt)^(1/4))*5

  pdf(output_file, width=wdth,height=hgt, pointsize = ptsize)  # was 50
  for (j in seq_along(params)){
    Nclust = nrow(params[[j]]$annot)
    ##prepare matrix for plotting
    M = as.matrix(params[[j]]$counts[1:Nclust,-(1:2)])
    rownames(M) = paste0("CL",rownames(M))
    Mn1=(M)/apply(M,1,max)
    Mn2=M/max(M)
    ord1 = hclust(dist(Mn1),method = "ward.D")$order
    ord2 = hclust(dist(t(Mn2)))$order

    ploting_area_width = 3 + log10(Nclust)*3
    ploting_area_sides = 1.5
    legend_width = 3
    title_height = 0.5
    layout(matrix(c(0,0,0,3,0,2,0,3,0,1,0,3),ncol=4,byrow = TRUE),
           width=c(ploting_area_sides,ploting_area_width,ploting_area_sides, legend_width),
           height=c(title_height, 3,ncol(M)*0.5))
    par(xaxs='i', yaxs = 'i')
    par(las=2,mar=c(4,0,0,0),cex.axis=0.5)
    if (any(is.na(GS))){
      rectMap(Mn2[ord1,ord2],scale.by='row',col=params[[j]]$color[ord1], grid=TRUE)
    }else{
      Mn3 = t(t(M) * (GS[colnames(M),] / params[[j]]$input_read_counts))[ord1,ord2]
      print(Mn3)
      Mn3 = Mn3/max(Mn3)
      print(Mn3)

      rectMap(Mn3,scale.by='none',col=params[[j]]$color[ord1], grid=TRUE)
    }
    par(las=2,mar=c(1,0,1,0), mgp = c(2,1,0))
    barplot(params[[j]]$annot$Size[ord1], col = 1)
    mtext(side = 2, "Cluster size", las = 3, line = 2, cex = 0.5)
    mtext(side=3, names(params)[j], las=0, line=1)
    plot.new()
    legend("topleft", col= params[[j]]$legend$color, legend=params[[j]]$legend$name, pch=15, cex=0.7, bty="n", pt.cex=1)
  }

  st = dev.off()
}

rectMap=function(x,scale.by='row',col=1,xlab="",ylab="",grid=TRUE,axis_pos=c(1,4),cexx=NULL,cexy=NULL){
  if (scale.by=='row'){
                                        #x=(x)/rowSums(x)
    x=(x)/apply(x,1,sum)
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
  #mtext(side = 2, "Proportions of individual samples", las =0, line = line, cex = 0.5)
  s=x/2
  w = c(x)/2
  rect(coords[,1]-0.5,coords[,2]-s,coords[,1]+0.5,coords[,2]+s,col=col,border=NA)
  if (grid){
    abline(v=0:(nr)+.5,h=0:(nc)+.5,lty=2,col="#60606030",lwd=0.2)
  }
  box(col="#60606030",lty=2, lwd=0.2)
}

option_list <- list( 
  make_option(c("-c", "--cluster_table"), default=NA, type = "character",
              help="file from RepeatExplorer2 clustering - CLUSTER_TABLE.csv"),

  make_option(c("-m", "--comparative_counts"),default = NA,type = "character",
              help="file from RepeatExplorer2 output - COMPARATIVE_ANALYSIS_COUNTS.csv"),

  make_option(c("-o", "--output"), type="character",
              default="comparative_analysis_summary.pdf",
              help="File name for output figures (pdf document)"),
  make_option(c("-N", "--number_of_colors"), type="integer", default=10,
              help="number of unqique color used from plotting"),

  make_option(c("-g", "--genome_size"),default = NA,type = "character",
              help="file from genome sizes of species provided in tab delimited file in the format:

                   species_code1   GenomeSize1
                   species_code2   GenomeSize2
                   species_code3   GenomeSize3
                   species_code4   GenomeSize4

                provide the same codes for species as in file COMPARATIVE_ANALYSIS_COUNTS.csv. If the genome
                sizes are provided, --nuclear_only is used by default. 
    "),
  make_option(c("-n", "--nuclear_only"), default = FALSE, type="logical",
              action = "store_true",
              help="remove all non-nuclear sequences (organelle and contamination). ")
)


opt = parse_args(OptionParser(option_list = option_list))

if (any(is.na(c(opt$cluster_table, opt$comparative_counts)))){
  message("\nBoth files: CLUSTER_TABLE.csv and COMPARATIVE_ANALYSIS_COUNTS.csv must be provided\n")
  q()
}

if (!opt$number_of_colors %in% 1:20){
  message("number of color must be in range 1..20")
  stop()
}

if (!is.na(opt$genome_size)){
  GS = read.table(opt$genome_size, header=FALSE, as.is=TRUE, row.names = 1)
  opt$nuclear_only=TRUE
}else{
  GS = NA
}

plot_rect_map(opt$comparative_counts, opt$cluster_table, opt$output, GS)

## testing:
cluster_annotation = "example_data/example1_CLUSTER_TABLE.csv"
read_counts = "example_data/example1_COMPARATIVE_ANALYSIS_COUNTS.csv"
output_file = "figs/test.pdf"
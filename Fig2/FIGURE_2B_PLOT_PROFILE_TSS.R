#####################################################################################
#      FIGURE_2B_PLOT_PROFILE_TSS
#####################################################################################
#
# Author(s): David Depierre, Cuvier Lab
# dav.depierre@gmail.com
# whoami: depierre
#
# Date: 2023-01-23
# Last update: 2023-01-23
#
#####################################################################################
#
# Project: PapolG / Z1 / RBM26 / prompt
#
# Description: script that do avg prof around TSS
#
# R version 3.4.2 (2017-09-28)
#
#########################################################################################################################
## LIBRARY
#########################################################################################################################


library("htmlwidgets")
library("viridisLite")
library("seqplots")


library(dplyr) # used for recover names(coverages)
library(gsubfn) # used in bam2coverage
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
# library(seqplots)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(rtracklayer)
library(devtools)



#########################################################################################################################
## fucntion
#########################################################################################################################

split_GRanges_inList <- function(GR, NnamesSplit, Nsplit = NULL){
  # namesSplit is either a character vector (then  decreasinglyordered splitted with Nsplit)
  # either a list of character vectors (then a list of granges is created according to list of names )
  namesSplit = get(NnamesSplit)
  if(is.numeric(namesSplit)){
    namesSplit = names(namesSplit[order(namesSplit, decreasing=T)])
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.character(namesSplit)){
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.list(namesSplit)){
    GRList = list()
    GRList = unlist(lapply(namesSplit, function(subnames){GRList = c(GRList, GR[names(GR) %in% subnames])}))
    names(GRList) =  paste0("GR_",names(namesSplit))
  }
  return(GRList)
}


#####################################################################################-

seqPlotSDoutliers_scaleFact <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,anchor=10000,sd=c(F,10),err=T,type="pf",smooth=T,spar=0.7, KR = F,gnme=NA,ignore.strand=F, scalingF = c(1,1), colvec = c("black","firebrick2")){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }

  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type, xanchored=anchor)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      Value_per_bin = c(gpsa.mtx)
      means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      means[which(is.nan(means))] <- 0 # change NA in 0
      if(sd[1] %in% T){
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd[2]] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        if(smooth){
          means = smooth.spline(1:(length(means)), means, spar=spar)$y
        }
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
      }
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      if(my.bw == bw.n[[scalingF[1]]]){
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means*scalingF[2]      #16.33333*1.761243 # change the means vector from getPlotSetArray object + scale fact for comparison

      }else{
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      }
    }

  }
  file.remove(paste0(o.tmp))
  par(mfrow=c(1,1))
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend_ext = T)
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend=F)
  # par(mfrow=c(2,2))
  # plot(density(Value_per_bin))

}


#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

tmp <- create(paste0(workdir,"PROJET_MTREC_IGH/FIGURE/AVG_RPOFILE/TMPgetPlotSetArray"))
tmp <- paste0(workdir,"PROJET_MTREC_IGH/FIGURE/AVG_RPOFILE/TMPgetPlotSetArray")

# BIGWIG file of H3K27me3
bwdir = paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/BIGWIG/")

BW_Zj1_shZ1_PolII = paste0(bwdir, "Zj1_shZ1_PolII.RPGC.sinorm.bw")
BW_Zj2_shZ1_Z1 = paste0(bwdir, "Zj2_shZ1_Z1.RPGC.sinorm.bw")
BW_Zj5_shZ1_PapolG = paste0(bwdir, "Zj5_shZ1_PapolG.RPGC.sinorm.bw")
BW_Zj6_shZ1_RBM26 = paste0(bwdir, "Zj6_shZ1_RBM26.RPGC.sinorm.bw")
BW_Zj7_shC_PolII = paste0(bwdir, "Zj7_shC_PolII.RPGC.sinorm.bw")
BW_Zj8_shC_Z1 = paste0(bwdir, "Zj8_shC_Z1.RPGC.sinorm.bw")
BW_Zj11_shC_PapolG = paste0(bwdir, "Zj11_shC_PapolG.RPGC.sinorm.bw")
BW_Zj12_shC_RBM26 = paste0(bwdir, "Zj12_shC_RBM26.RPGC.sinorm.bw")


# GENES

GRCh38_gene = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/Homo_sapiens.GRCh38.77.chromosomes.gene.names.rds"))
GRCh38_gene_TSS = resize(GRCh38_gene, 1, fix="start")
names(GRCh38_gene_TSS) = paste0(names(GRCh38_gene_TSS), ".1")

GNref_active = scan(paste0(workdir, "PROJET_MTREC_IGH/DATA/active_gene_union_Z1_MTR4_WT.txt"), what="character")
ensembl2entrez = read.table(paste0(workdir,"PROJET_MTREC_IGH/DATA/ensembl_to_entrez.txt"), header=TRUE)
GNref_active = as.character(ensembl2entrez[ ensembl2entrez$hg38_genes.hgnc_symbol %in% GNref_active, 1])

GRCh38_gene_TSS_active = GRCh38_gene_TSS[names(GRCh38_gene_TSS) %in% GNref_active]

#####################################################################
# PLOT
#####################################################################
GR_list_toPlot = c("GRCh38_gene_TSS_active")

for(GR in GR_list_toPlot){
  pdf(paste0(workdir,"PROJET_MTREC_IGH/FIGURE/AVG_RPOFILE/AVG_PROF_TSS_CHIPSEQ_RPGC_" ,GR,"_SD3_5000bp_smooth20.pdf"))
  par(lwd=2)
  seqPlotSDoutliers_scaleFact(c(BW_Zj8_shC_Z1),tmp,GR,c(-0.5,5),c(5000,5000),bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), colvec = c("#04be9eff"))
  
  seqPlotSDoutliers_scaleFact(c(BW_Zj7_shC_PolII),tmp,GR,c(-0.5,25),c(5000,5000),bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), colvec = c("#f5a623ff"))

  seqPlotSDoutliers_scaleFact(c(BW_Zj11_shC_PapolG),tmp,GR,c(-0.5,8),c(5000,5000),bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), colvec = c("#74a808ff"))

  seqPlotSDoutliers_scaleFact(c(BW_Zj12_shC_RBM26),tmp,GR,c(-0.5,4),c(5000,5000),bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), colvec = c("#bd10e0ff"))
  dev.off()
}




#end

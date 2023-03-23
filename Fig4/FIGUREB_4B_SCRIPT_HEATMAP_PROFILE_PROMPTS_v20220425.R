#####################################################################################
#      FIGUREB_4B_SCRIPT_HEATMAP_PROFILE_PROMPTS_v20220425
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
# Description: script that do heatmap around TSS
#
# R version 3.4.2 (2017-09-28)
#
#########################################################################################################################
## LIBRARY
#########################################################################################################################

biocLite("gridBase")

library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(gplots)
'%ni%' = Negate('%in%')
require(BiocGenerics)
require(parallel)
library(gsubfn)
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(rtracklayer)
library(nucleR)
library(gridBase)
library(dplyr)


#####################################################################################################################################
# LOAD FUNCTON
#####################################################################################################################################


source(paste0(workdir, "functionR/Script_HEATMAP_profile.R"))
# from https://github.com/ddepierre/analysis_scripts_H3K36_H3K27me3/blob/main/FIGURES/R_PLOT_FUNCTION/Script_HEATMAP_profile.R

#################################################################################################################################

## FUNCTION TO COMPUTE DIFFERENTIAL MATRIX


matMeans <- function(X,Y){ mean(c(X,Y)) }

matReplaceNA = function(M){
  Mv = c(M)
  Mv[which(is.na(Mv))] = 0
  Mv = matrix(Mv, ncol=1000)
  rownames(Mv) = rownames(M)
  return(Mv)
}


computeZscore_PROFMAT = function(PROF_KD, PROF_CTRL){
  PROF_KD = matReplaceNA(PROF_KD)
  PROF_CTRL = matReplaceNA(PROF_CTRL)
  ## Norm on 500bp between -5kb and -4.5kb
  # PROF_KD = PROF_KD*sum(c(PROF_CTRL[,1:1000]))/sum(c(PROF_KD[,1:1000]))
  ZSCORE = (PROF_KD-PROF_CTRL[rownames(PROF_KD),])/sqrt(matrix(mapply(matMeans, PROF_KD, PROF_CTRL), ncol=1000))
  ZSCORE = matReplaceNA(ZSCORE)
  rownames(ZSCORE) = rownames(PROF_KD)
  return(ZSCORE)
}


#####################################################################################################################################
#### LOAD DATA
#####################################################################################################################################
#### LOAD PROFILE_MATRIX

################################################################################################################################################################ RNA-seq
workdir = "/home/depierre/Bureau/work_cperrois/" # CUVIER10
workdir = "/home/cperrois/work/" # GENOTOUL

profdir = paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/PROFILE_MATRIX_SENS_ANTISENS/")
# Load data

WT_JENSEN_SRR3757038_hg38_ANTISENS = readRDS(paste0(profdir, "WT_JENSEN_SRR3757038_hg38_profmat_ANTISENS.RDS"))
WT_JENSEN_SRR3757040_hg38_ANTISENS = readRDS(paste0(profdir, "WT_JENSEN_SRR3757040_hg38_profmat_ANTISENS.RDS"))
Z1_JENSEN_SRR3757088_hg38_ANTISENS = readRDS(paste0(profdir, "Z1_JENSEN_SRR3757088_hg38_profmat_ANTISENS.RDS"))
Z1_JENSEN_SRR3757090_hg38_ANTISENS = readRDS(paste0(profdir, "Z1_JENSEN_SRR3757090_hg38_profmat_ANTISENS.RDS"))


#### LOAD QUANTIF CHIPSEQ TSS pm500bp


ZSCORE_TSS_shC_vs_shZ1_PapolG = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/QUANTIF/TSS/ZSCORE_TSS_shC_vs_shZ1_PapolG.RDS"))
ZSCORE_TSS_shC_vs_shZ1_RBM26 = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/QUANTIF/TSS/ZSCORE_TSS_shC_vs_shZ1_RBM26.RDS"))
ZSCORE_TSS_shC_vs_shZ1_Z1 = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/QUANTIF/TSS/ZSCORE_TSS_shC_vs_shZ1_Z1.RDS"))

names(ZSCORE_TSS_shC_vs_shZ1_PapolG) = unlist(lapply(strsplit(names(ZSCORE_TSS_shC_vs_shZ1_PapolG), "\\."), `[[`, 1))
names(ZSCORE_TSS_shC_vs_shZ1_RBM26) = unlist(lapply(strsplit(names(ZSCORE_TSS_shC_vs_shZ1_RBM26), "\\."), `[[`, 1))
names(ZSCORE_TSS_shC_vs_shZ1_Z1) = unlist(lapply(strsplit(names(ZSCORE_TSS_shC_vs_shZ1_Z1), "\\."), `[[`, 1))



#####################################################################################################################################
#### PREPARE DATA TO PLOT
#####################################################################################################################################

################################################################################################################################################################ ChIP-seq
# GET PROFMAT ZSCORE

ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS = computeZscore_PROFMAT(Z1_JENSEN_SRR3757088_hg38_ANTISENS, WT_JENSEN_SRR3757038_hg38_ANTISENS)
ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS = computeZscore_PROFMAT(Z1_JENSEN_SRR3757090_hg38_ANTISENS, WT_JENSEN_SRR3757040_hg38_ANTISENS)

saveRDS(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS, paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS.RDS"))
saveRDS(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS, paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS.RDS"))

ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS.RDS"))
ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS.RDS"))

rownames(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS) = paste0(rownames(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS),".1")
rownames(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS) = paste0(rownames(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS),".1")

rownames(WT_JENSEN_SRR3757038_hg38_ANTISENS) = paste0(rownames(WT_JENSEN_SRR3757038_hg38_ANTISENS),".1")

ZSCORE_TSS_shC_vs_shZ1_PapolG = sort(ZSCORE_TSS_shC_vs_shZ1_PapolG, decreasing=F)
ZSCORE_TSS_shC_vs_shZ1_RBM26 = sort(ZSCORE_TSS_shC_vs_shZ1_RBM26, decreasing=F)
ZSCORE_TSS_shC_vs_shZ1_Z1 = sort(ZSCORE_TSS_shC_vs_shZ1_Z1, decreasing=F)


#####################################################################################################################################
#### PLOT HEATMAP
#####################################################################################################################################


# average replicates
ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS = list(ZSCORE_PROFMAT_Z1_JENSEN_SRR3757088_vs_WT_JENSEN_SRR3757038_hg38_ANTISENS,ZSCORE_PROFMAT_Z1_JENSEN_SRR3757090_vs_WT_JENSEN_SRR3757040_hg38_ANTISENS) %>% simplify2array %>% apply(.,1:2,function(x){mean(x,na.rm=T)})
saveRDS(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS, paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS.RDS"))

# smooth
ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth = t(apply(matReplaceNA(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS), 1, function(X){smooth.spline(X, spar = 0.4)$y}))
saveRDS(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth, paste0(profdir, "ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth.RDS"))



pdf(paste0(workdir, "PROJET_MTREC_IGH/FIGURE/HEATMAP_PROFILE/PROMPTS/ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth_sortby_ZscoreQuantifKD_WT.pdf"))
rangeheatmap = c(1:1000)
heatMatrixMat(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth[names(ZSCORE_TSS_shC_vs_shZ1_Z1),rangeheatmap],RangeValue=c(-1,1), main = "ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_Z1", colorSet = c("#ffffff","#ffffff","#000000"))
heatMatrixMat(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth[names(ZSCORE_TSS_shC_vs_shZ1_PapolG),rangeheatmap],RangeValue=c(-1,1), main = "ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_PapolG", colorSet = c("#ffffff","#ffffff","#000000"))
heatMatrixMat(ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth[names(ZSCORE_TSS_shC_vs_shZ1_RBM26),rangeheatmap],RangeValue=c(-1,1), main = "ZSCORE_PROFMAT_Z1_JENSEN_88_90_vs_WT_JENSEN_38_40_hg38_ANTISENS_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_RBM26", colorSet = c("#ffffff","#ffffff","#000000"))

dev.off()



###END

#####################################################################################
#      FIGURE_3C_SCRIPT_HEATMAP_PROFILE_CHISEQ_TSS
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


#####################################################################################################################################
#### LOAD DATA
#####################################################################################################################################
#### LOAD PROFILE_MATRIX

profdir = paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PROFILE_MATRIX/")
PROFMAT_Zj7_shC_PolII = readRDS(paste0(profdir, "Zj7_shC_PolII.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj8_shC_Z1 = readRDS(paste0(profdir, "Zj8_shC_Z1.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj9_shC_MTR4 = readRDS(paste0(profdir, "Zj9_shC_MTR4.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj10_shC_Rad21 = readRDS(paste0(profdir, "Zj10_shC_Rad21.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj11_shC_PapolG = readRDS(paste0(profdir, "Zj11_shC_PapolG.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj12_shC_RBM26 = readRDS(paste0(profdir, "Zj12_shC_RBM26.RPGC.sinorm.profmat.RDS"))


PROFMAT_Zj1_shZ1_PolII = readRDS(paste0(profdir, "Zj1_shZ1_PolII.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj2_shZ1_Z1 = readRDS(paste0(profdir, "Zj2_shZ1_Z1.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj3_shZ1_MTR4 = readRDS(paste0(profdir, "Zj3_shZ1_MTR4.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj4_shZ1_Rad21 = readRDS(paste0(profdir, "Zj4_shZ1_Rad21.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj5_shZ1_PapolG = readRDS(paste0(profdir, "Zj5_shZ1_PapolG.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj6_shZ1_RBM26 = readRDS(paste0(profdir, "Zj6_shZ1_RBM26.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj7_shC_PolII = readRDS(paste0(profdir, "Zj7_shC_PolII.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj8_shC_Z1 = readRDS(paste0(profdir, "Zj8_shC_Z1.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj9_shC_MTR4 = readRDS(paste0(profdir, "Zj9_shC_MTR4.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj10_shC_Rad21 = readRDS(paste0(profdir, "Zj10_shC_Rad21.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj11_shC_PapolG = readRDS(paste0(profdir, "Zj11_shC_PapolG.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj12_shC_RBM26 = readRDS(paste0(profdir, "Zj12_shC_RBM26.RPGC.sinorm.profmat.RDS"))


#### LOAD QUANTIF CHIPSEQ TSS pm500bp

ZSCORE_TSS_shC_vs_shZ1_Z1 = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/QUANTIF/TSS/ZSCORE_TSS_shC_vs_shZ1_Z1.RDS"))



#####################################################################################################################################
#### PREPARE DATA TO PLOT
#####################################################################################################################################

GRCh38_gene = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/Homo_sapiens.GRCh38.77.chromosomes.gene.names.rds"))
GRCh38_gene_TSS = resize(GRCh38_gene, 1000, fix="start")
names(GRCh38_gene_TSS) = paste0(names(GRCh38_gene_TSS), ".1")

GNref_active = scan(paste0(workdir, "PROJET_MTREC_IGH/DATA/active_gene_union_Z1_MTR4_WT.txt"), what="character")
ensembl2entrez = read.table(paste0(workdir,"PROJET_MTREC_IGH/DATA/ensembl_to_entrez.txt"), header=TRUE)
GNref_active = as.character(ensembl2entrez[ ensembl2entrez$hg38_genes.hgnc_symbol %in% GNref_active, 1])

GRCh38_gene_TSS_active = GRCh38_gene_TSS[names(GRCh38_gene_TSS) %in% GNref_active]


################################################################################################################################################################ ChIP-seq
# GET PROFMAT ZSCORE


ZSCORE_PROFMAT_Zj1_shZ1_Zj7_shC_PolII = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj1_shZ1_Zj7_shC_PolII.RDS"))
ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1 = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1.RDS"))
ZSCORE_PROFMAT_Zj3_shZ1_Zj9_shC_MTR4 = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj3_shZ1_Zj9_shC_MTR4.RDS"))
ZSCORE_PROFMAT_Zj4_shZ1_Zj10_shC_Rad21 = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj4_shZ1_Zj10_shC_Rad21.RDS"))
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG.RDS"))
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26 = readRDS(paste0(profdir, "ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26.RDS"))

#####################################################################################################################################
#### PREPARE DATA
#####################################################################################################################################



ZSCORE_TSS_shC_vs_shZ1_Z1 = ZSCORE_TSS_shC_vs_shZ1_Z1[names(ZSCORE_TSS_shC_vs_shZ1_Z1) %in% GNref_active]

Q_filter = names(Q_Zj8_shC_Z1_TSS[Q_Zj8_shC_Z1_TSS>1 & Q_Zj11_shC_PapolG_TSS>1 & Q_Zj12_shC_RBM26_TSS>1])

ZSCORE_TSS_shC_vs_shZ1_Z1 = sort(ZSCORE_TSS_shC_vs_shZ1_Z1, decreasing=T)



ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f = ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1[rownames(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1) %in% GNref_active,] 
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f = ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG[rownames(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG) %in% GNref_active,] 
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f = ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26[rownames(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26) %in% GNref_active,] 

ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f = ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f[rownames(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f) %in% Q_filter,] 
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f = ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f[rownames(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f) %in% Q_filter,] 
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f = ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f[rownames(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f) %in% Q_filter,] 

ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f[is.infinite(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f)] <-0 
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f[is.infinite(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f)] <-0 
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f[is.infinite(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f)] <-0 

ZSCORE_TSS_shC_vs_shZ1_Z1 = ZSCORE_TSS_shC_vs_shZ1_Z1[names(ZSCORE_TSS_shC_vs_shZ1_Z1) %in% GNref_active]
ZSCORE_TSS_shC_vs_shZ1_Z1 = ZSCORE_TSS_shC_vs_shZ1_Z1[names(ZSCORE_TSS_shC_vs_shZ1_Z1) %in% Q_filter]


ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f_smth = t(apply(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_f, 1, function(X){smooth.spline(X, spar = 0.4)$y}))
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f_smth = t(apply(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_f, 1, function(X){smooth.spline(X, spar = 0.4)$y}))
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f_smth = t(apply(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_f, 1, function(X){smooth.spline(X, spar = 0.4)$y}))




########################################################################################################################
# filter by loss of Z1
########################################################################################################################

ZSCORE_TSS_shC_vs_shZ1_Z1_loss = ZSCORE_TSS_shC_vs_shZ1_Z1[ZSCORE_TSS_shC_vs_shZ1_Z1<0]

ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss = ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1[rownames(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1) %in% names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),] 
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss = ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG[rownames(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG) %in% names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),] 
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss = ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26[rownames(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26) %in% names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),] 



ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss[is.infinite(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss)] <-0 
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss[is.infinite(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss)] <-0 
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss[is.infinite(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss)] <-0 


ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss_smth = t(apply(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss, 1, function(X){smooth.spline(X, spar = 0.4)$y}))
ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss_smth = t(apply(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss, 1, function(X){smooth.spline(X, spar = 0.4)$y}))
ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss_smth = t(apply(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss, 1, function(X){smooth.spline(X, spar = 0.4)$y}))

ZSCORE_TSS_shC_vs_shZ1_Z1_loss = ZSCORE_TSS_shC_vs_shZ1_Z1_loss[names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss) %in% GNref_active]


#####################################################################################################################################
#### PLOT HEATMAP
#####################################################################################################################################



pdf(paste0(workdir, "PROJET_MTREC_IGH/FIGURE/HEATMAP_PROFILE/CHIPSEQ/ZSCORE_PROFMAT_smth_CHIPSEQ_shZ1_vs_shC_sortby_ZSCORE_TSS_shC_vs_shZ1_onlyZ1loss_Z1_20220425.pdf"))
rangeheatmap = c(1:1000)
heatMatrixMat(ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss_smth[names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),rangeheatmap],RangeValue=c(-3,3), main = "ZSCORE_PROFMAT_Zj2_shZ1_Zj8_shC_Z1_fZ1loss_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_Z1_loss")

heatMatrixMat(ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss_smth[names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),rangeheatmap],RangeValue=c(-3,3), main = "ZSCORE_PROFMAT_Zj5_shZ1_Zj11_shC_PapolG_fZ1loss_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_Z1_loss")

heatMatrixMat(ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss_smth[names(ZSCORE_TSS_shC_vs_shZ1_Z1_loss),rangeheatmap],RangeValue=c(-2.5,2.5), main = "ZSCORE_PROFMAT_Zj6_shZ1_Zj12_shC_RBM26_fZ1loss_smth", legend.name = "ZSCORE_TSS_shC_vs_shZ1_Z1_loss")

dev.off()



########################################################################################################################
# END
########################################################################################################################







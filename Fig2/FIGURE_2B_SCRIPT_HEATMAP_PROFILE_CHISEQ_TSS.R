#####################################################################################
#      FIGURE_2B_SCRIPT_HEATMAP_PROFILE_CHISEQ_TSS
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
PROFMAT_Zj11_shC_PapolG = readRDS(paste0(profdir, "Zj11_shC_PapolG.RPGC.sinorm.profmat.RDS"))
PROFMAT_Zj12_shC_RBM26 = readRDS(paste0(profdir, "Zj12_shC_RBM26.RPGC.sinorm.profmat.RDS"))


#### LOAD QUANTIF PROMPTS 2KB
qtdir = paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/")

Q_WT_JENSEN_SRR3757038_hg38_ANTISENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/WT_JENSEN_SRR3757038_hg38_ANTISENS_readsCounts_PROM_2KB.RDS"))
Q_WT_JENSEN_SRR3757038_hg38_SENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/WT_JENSEN_SRR3757038_hg38_SENS_readsCounts_PROM_2KB.RDS"))
Q_WT_JENSEN_SRR3757040_hg38_ANTISENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/WT_JENSEN_SRR3757040_hg38_ANTISENS_readsCounts_PROM_2KB.RDS"))
Q_WT_JENSEN_SRR3757040_hg38_SENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/WT_JENSEN_SRR3757040_hg38_SENS_readsCounts_PROM_2KB.RDS"))
Q_Z1_JENSEN_SRR3757088_hg38_ANTISENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/Z1_JENSEN_SRR3757088_hg38_ANTISENS_readsCounts_PROM_2KB.RDS"))
Q_Z1_JENSEN_SRR3757088_hg38_SENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/Z1_JENSEN_SRR3757088_hg38_SENS_readsCounts_PROM_2KB.RDS"))
Q_Z1_JENSEN_SRR3757090_hg38_ANTISENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/Z1_JENSEN_SRR3757090_hg38_ANTISENS_readsCounts_PROM_2KB.RDS"))
Q_Z1_JENSEN_SRR3757090_hg38_SENS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/QUANTIF/PROMOTER/Z1_JENSEN_SRR3757090_hg38_SENS_readsCounts_PROM_2KB.RDS"))

#### LOAD QUANTIF CHIPSEQ TSS pm500bp

Q_Zj7_shC_PolII_TSS = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/QUANTIF/TSS/Zj7_shC_PolII_nosinorm_readsCounts_TSS.RDS"))



#####################################################################################################################################
#### PREPARE DATA TO PLOT
#####################################################################################################################################
GNref_active = scan(paste0(workdir, "PROJET_MTREC_IGH/DATA/active_gene_union_Z1_MTR4_WT.txt"), what="character")
ensembl2entrez = read.table(paste0(workdir,"PROJET_MTREC_IGH/DATA/ensembl_to_entrez.txt"), header=TRUE)
GNref_active = as.character(ensembl2entrez[ ensembl2entrez$hg38_genes.hgnc_symbol %in% GNref_active, 1])

Q_Zj7_shC_PolII_TSS = Q_Zj7_shC_PolII_TSS[names(Q_Zj7_shC_PolII_TSS) %in% GNref_active]

Q_Zj7_shC_PolII_TSS = sort(Q_Zj7_shC_PolII_TSS, decreasing=T)


#####################################################################################################################################
#### PLOT HEATMAP
#####################################################################################################################################

pdf(paste0(workdir, "PROJET_MTREC_IGH/FIGURE/HEATMAP_PROFILE/CHIPSEQ/PROFMAT_CHIPSEQ_shC_sortby_Q_Zj7_shC_PolII_TSS.pdf"))
rangeheatmap = c(1:1000)
heatMatrixMat(PROFMAT_Zj7_shC_PolII[names(Q_Zj7_shC_PolII_TSS),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_Zj7_shC_PolII", legend.name = "Q_Zj7_shC_PolII_TSS")
heatMatrixMat(PROFMAT_Zj8_shC_Z1[names(Q_Zj7_shC_PolII_TSS),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_Zj8_shC_Z1", legend.name = "Q_Zj7_shC_PolII_TSS")
heatMatrixMat(PROFMAT_Zj11_shC_PapolG[names(Q_Zj7_shC_PolII_TSS),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_Zj11_shC_PapolG", legend.name = "Q_Zj7_shC_PolII_TSS")
heatMatrixMat(PROFMAT_Zj12_shC_RBM26[names(Q_Zj7_shC_PolII_TSS),rangeheatmap],winsorize=c(5,95), main = "PROFMAT_Zj12_shC_RBM26", legend.name = "Q_Zj7_shC_PolII_TSS")
dev.off()

###END

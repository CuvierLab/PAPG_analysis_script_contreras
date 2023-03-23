#####################################################################################
#      FIGURE_2C_VENNDIAGRAM
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
# Description: script that do venn diagram between TSS with peaks
#
# R version 3.4.2 (2017-09-28)
#
#########################################################################################################################
## LIBRARY
#########################################################################################################################

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
library(gplots)
library(Vennerable)
library(SuperExactTest)

#########################################################################################################################
## fucntion
#########################################################################################################################

plotVENNerable = function(data_to_venn, Nb_REF, outdir){
	vennNAME = paste(names(data_to_venn), collapse="_")
	Venn_ovlp = Venn(data_to_venn)
	color_Venn = VennThemes(compute.Venn(Venn_ovlp))
  #set1
  color_Venn[["Set"]][["Set1"]]$col = "#04be9e"
  color_Venn[["SetText"]][["Set1"]]$col = "#04be9e"
  color_Venn[["Face"]][["100"]]$fill = "#93d9cd"
  color_Venn[["Face"]][["100-1"]]$fill = "#93d9cd"
  #set2
  color_Venn[["Set"]][["Set2"]]$col = "#bd10e0"
  color_Venn[["SetText"]][["Set2"]]$col = "#bd10e0"
  color_Venn[["Face"]][["010"]]$fill = "#e29cf0"
  color_Venn[["Face"]][["010-1"]]$fill = "#e29cf0"
  #set3
  color_Venn[["Set"]][["Set3"]]$col = #74a808"
  color_Venn[["SetText"]][["Set3"]]$col = #74a808"
  color_Venn[["Face"]][["001"]]$fill = "#ceeb91"
  color_Venn[["Face"]][["001-1"]]$fill = "#ceeb91"
  #set1 n set2
  color_Venn[["Face"]][["110"]]$fill = "#a292e8"
  color_Venn[["Face"]][["110-1"]]$fill = "#a292e8"
  #set1 n set3
  color_Venn[["Face"]][["101"]]$fill = "#4ca696"
  color_Venn[["Face"]][["101-1"]]$fill = "#4ca696"
  #set2 n set3
	color_Venn[["Face"]][["011"]]$fill = "#a476ad"
  color_Venn[["Face"]][["011-1"]]$fill = "#a476ad"
  #TRIPLE_inter
	color_Venn[["Face"]][["111"]]$fill = "#8649b8"
  color_Venn[["Face"]][["111-1"]]$fill = "#8649b8"
	fishertest = fisher_namesList(data_to_venn,data_to_venn,Nb_REF)
	# Supertest for triple intersection
	ResSuperTest = supertest(data_to_venn, n =Nb_REF)
	SupetTestToPrint = as.matrix(ResSuperTest$P.value)
	colnames(SupetTestToPrint) = "intersect_pval"
	########### PLOT #########################
	pdf(paste0(outdir, "PLOTvenn_",vennNAME,".pdf"))
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE)
	par(mfrow=c(2,1))
	textplot(fishertest$p.value,  valign="top")
	textplot(fishertest$odds.ratio)
	textplot(SupetTestToPrint)
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE, show = list(SetLabels = FALSE, FaceText = ""))
	dev.off()
}



#########################################################################################################################
## load data
#########################################################################################################################


GRCh38_gene = readRDS(paste0(workdir, "PROJET_MTREC_IGH/DATA/RNASEQ_JENSEN/Homo_sapiens.GRCh38.77.chromosomes.gene.names.rds"))
GRCh38_gene_TSS = resize(GRCh38_gene, 1000, fix="start")
names(GRCh38_gene_TSS) = paste0(names(GRCh38_gene_TSS), ".1")

GNref_active = scan(paste0(workdir, "PROJET_MTREC_IGH/DATA/active_gene_union_Z1_MTR4_WT.txt"), what="character")
ensembl2entrez = read.table(paste0(workdir,"PROJET_MTREC_IGH/DATA/ensembl_to_entrez.txt"), header=TRUE)
GNref_active = as.character(ensembl2entrez[ ensembl2entrez$hg38_genes.hgnc_symbol %in% GNref_active, 1])

GRCh38_gene_TSS_active = GRCh38_gene_TSS[names(GRCh38_gene_TSS) %in% GNref_active]

shC_MTR4_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_MTR4_fdr4_pk.bed"))
shC_PapolG_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_PapolG_fdr4_pk.bed"))
shC_PolII_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_PolII_fdr4_pk.bed"))
shC_Rad21_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_Rad21_fdr4_pk.bed"))
shC_RBM26_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_RBM26_fdr4_pk.bed"))
shC_Z1_fdr4_pk = rtracklayer::import.bed(paste0(workdir, "PROJET_MTREC_IGH/DATA/CHIPSEQ/PEAK_CALLING/shC_Z1_fdr4_pk.bed"))



#########################################################################################################################
## filter active TSS
#########################################################################################################################

TSS_shC_MTR4_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_MTR4_fdr4_pk, ignore.strand=T, minoverlap=50)@from])
TSS_shC_PapolG_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_PapolG_fdr4_pk, ignore.strand=T, minoverlap=50)@from])
TSS_shC_PolII_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_PolII_fdr4_pk, ignore.strand=T, minoverlap=50)@from])
TSS_shC_Rad21_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_Rad21_fdr4_pk, ignore.strand=T, minoverlap=50)@from])
TSS_shC_RBM26_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_RBM26_fdr4_pk, ignore.strand=T, minoverlap=50)@from])
TSS_shC_Z1_fdr4_pk = names(GRCh38_gene_TSS_active[findOverlaps(GRCh38_gene_TSS_active, shC_Z1_fdr4_pk, ignore.strand=T, minoverlap=50)@from])



#########################################################################################################################
## plot
#########################################################################################################################


List1 = list(TSS_ACTIVE_shC_Z1_fdr4_pk = TSS_shC_Z1_fdr4_pk, TSS_ACTIVE_shC_RBM26_fdr4_pk = TSS_shC_RBM26_fdr4_pk, TSS_ACTIVE_shC_PapolG_fdr4_pk = TSS_shC_PapolG_fdr4_pk)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GRCh38_gene_TSS_active), outdir = outfig)

#########################################################################################################################
## END
#########################################################################################################################


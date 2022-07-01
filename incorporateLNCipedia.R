# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: Nov. 1st, 2018
# Update: June 2022 to incorporate GENCODE40

#Purpose: To create a nonredundant annotation set between LNCipedia and GENCODE

#Imports: GENCODE40 and LNCipedia5.2 (high confidence) gtf files
#Exports: 
 #venn diagram summary
 #gtf files to be further combined in bash script

# import packages -------------------------------------------------------------------------------
library(rtracklayer)
library(VennDiagram)
library(tidyverse)
library(gridExtra)

# import GENCODE40 annotation -------------------------------------------------------------------
gencode.v40.gtf <- "gencode.v40.primary_assembly.annotation.gtf"
gencode40 <- rtracklayer::import(gencode.v40.gtf, format = "gtf")
gencode40.df <- as.data.frame(gencode40)

gencode40.genes <- subset(gencode40.df, type == "gene") #list genes only in gencode annotation
gencode40.lnc.genes <- subset(gencode40.genes, gene_type == "lncRNA") #list genes only in gencode annotation

# import LNCipedia 5.2 (high confidence) annotation ---------------------------------------------
lncipedia.gtf <- "lncipedia_5_2_hc_hg38.formatted.gtf"
lncipedia5.2 <- rtracklayer::import(lncipedia.gtf, format = "gtf")
lncipedia5.2.df <- as.data.frame(lncipedia5.2)

# import 92 to 104 table, remove retired IDs
ens104 <- read.csv("ensembl104GeneUpdate.csv", header = TRUE)
geneIDs.104 <- ens104[grep("retired", ens104$New.stable.ID, invert = TRUE),]

# create list of intersection, unique, etc. -------------------------------------------------------
gencode40.ensBase <- gsub("\\.(.*)","" , gencode40.genes$gene_id)
gencode40.lnc.ensBase <- gsub("\\.(.*)","" , gencode40.lnc.genes$gene_id)
lncipedia.ensBase <- gsub("\\.(.*)","" , geneIDs.104$ensemblGeneID)

lncipediaIntersection <- gencode40.genes[gencode40.ensBase %in% lncipedia.ensBase,] #11356 lncipedia ENSG IDs in gencode 40

gencodeUnique <- gencode40.genes[!gencode40.ensBase %in% lncipedia.ensBase,] #50242 gencode genes (all biotypes) not in lncipedia
gencodeUnique.lncRNA <- gencode40.lnc.genes[!gencode40.lnc.ensBase %in% lncipedia.ensBase,] #6410 gencode genes (matching lncRNA biotype) not in lncipedia

gencodeTable <- table(gencode40.genes$gene_type)["lncRNA"]
lncipediaTable <- table(lncipediaIntersection$gene_type) #biotype of lncRNA in intersection of two annotations

gencodeKeep <- gencodeUnique #50242

#same list for transcripts
gencodeTranscriptTable <- table(subset(gencode40.df, type == "transcript")$transcript_type)
#lncipediaTranscriptIntersection <- subset(gencode40.df, type == "transcript")[gencode40.genes$transcript_id %in% geneIDs.104$ensemblGeneID,] #13206 lncipedia ENSG IDs in gencode 40
lncipediaTranscriptIntersection <- lncipediaIntersection[lncipediaIntersection$transcript_type %in% lncBiotypes,] #13206: same number 

# make venn diagram highlighting changes -----------------------------------------------------------
# calculate gene totals
# note these include an additional 8848 genes found via intersectBed
gencode.gene.total <- sum(gencodeTable) + 8848
lncipedia.gene.total <- length(unique(lncipedia5.2$gene_id))
intersection.gene.total <- sum(lncipediaTable) + 8848

dir.create("img")
dev.off()
tmp <- draw.pairwise.venn(gencode.gene.total,
                          lncipedia.gene.total,
                          intersection.gene.total, #additional 8848 genes from intersectBed
                          category = c("GENCODE 40", "lncipedia 5.2"),
                          fill = c("#00AFBB", "#E7B800"),
                          cex = c(4,4,4),
                          cat.cex = c(4,4),
                          alpha = rep(0.5, 2), cat.pos = c(0, -40), cat.dist = c(0.03, 0.10))


png("img/lncRNA_gene_venn.png", width = 800, height = 800)
grid.draw(tmp)
dev.off()



# tmp2 <- draw.pairwise.venn(sum(gencodeTranscriptTable),
#                           length(unique(lncipedia5.2$transcript_id)),
#                           sum(lncipediaTable),
#                           category = c("GENCODE 40", "lncipedia 5.2"),
#                           fill = c("#00AFBB", "#E7B800"),
#                           alpha = rep(0.5, 2), cat.pos = c(-105, 135), cat.dist = rep(0.08, 2))
# 
# png("img/lncRNA_transcript_venn.png")
# grid.draw(tmp)
# dev.off()

#create gencode gtf with lncipedia genes removed and export -----------------------------------------
gtfNew <- gencode40[gencode40$gene_id %in% gencodeUnique$gene_id,]

export(gtfNew, 'gencode.v40.primary_assembly.annotation.lncipediaRemoved.gtf')

#add gene level info to lncipedia annotation since it is exons only and export ----------------------
#see:
#https://bioinformatics.stackexchange.com/questions/4253/get-gene-lines-from-gtf-file
grl <- split(lncipedia5.2, lncipedia5.2$gene_id)  # This produces a GRangesList
grl <- endoapply(grl, function(x) {
  foo = x[1]
  foo$type = 'gene'
  start(foo) = min(start(x))
  end(foo) = max(end(x))
  return(c(foo, x))
})
gr <- unlist(grl)

export.gff2(subset(gr, type == "gene"), "lncipedia_5_2_hc_hg38.genes.gtf")

#add transcript level info to lncipedia annotation and export ----------------------------------------
trl <- split(lncipedia5.2, lncipedia5.2$transcript_id)  # This produces a GrangesList
trl <- endoapply(trl, function(x) {
  foo = x[1]
  foo$type = 'transcript'
  start(foo) = min(start(x))
  end(foo) = max(end(x))
  return(c(foo, x))
})
tr <- unlist(trl)

export.gff2(unlist(trl), "lncipedia_5_2_hc_hg38.transcripts.gtf")



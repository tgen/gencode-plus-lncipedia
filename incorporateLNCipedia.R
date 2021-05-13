# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: Nov. 1st, 2018

#Purpose: To create a nonredundant annotation set between LNCipedia and GENCODE

#Imports: GENCODE29 and LNCipedia5.2 (high confidence) gtf files
#Exports: 
 #venn diagram summary
 #gtf files to be further combined in bash script

# import packages -------------------------------------------------------------------------------
library(rtracklayer)
library(VennDiagram)
library(tidyverse)

# import GENCODE29 annotation -------------------------------------------------------------------
gencode.v29.gtf <- "gencode.v29.primary_assembly.annotation.gtf"
gencode29 <- rtracklayer::import(gencode.v29.gtf, format = "gtf")
gencode29.df <- as.data.frame(gencode29)

gencode29.genes <- subset(gencode29.df, type == "gene") #list genes only in gencode annotation

# import LNCipedia 5.2 (high confidence) annotation ---------------------------------------------
lncipedia.gtf <- "lncipedia_5_2_hc_hg38.formatted.gtf"
lncipedia5.2 <- rtracklayer::import(lncipedia.gtf, format = "gtf")
lncipedia5.2.df <- as.data.frame(lncipedia5.2)

# import LNCipedia 5.2 to Ensembl 92 table
geneIDs <- read.table("lncipedia_5_2_ensembl_92_genes.txt", header = TRUE)

# import 92 to 94 table, remove retired IDs
ens94 <- read.csv("ens94GeneUpdate.csv", header = FALSE)
ens94 <- ens94[grep("retired", ens94$V2, invert = TRUE),]
ens94$V2 <- as.character(ens94$V2)

# update to ensembl 94
geneIDs.94 <- left_join(geneIDs, ens94[,c(1:2)], by = c("ensemblGeneID" = "V1"))
geneIDs.94$V2 <- coalesce(geneIDs.94$V2, geneIDs.94$ensemblGeneID)
names(geneIDs.94) <- c("lncipediaGeneID", "ensembl92GeneID", "ensembl94GeneID")


# list lncRNA biotypes in GENCODE29 ---------------------------------------------------------------
#NOTE: not including TEC (to be experimentally confirmed) in list
lncBiotypes <- c("3prime_overlapping_ncRNA", "antisense",
                 "bidirectional_promoter_lncRNA", "lincRNA",
                 "non_coding", "processed_transcript",
                 "macro_lncRNA", "sense_intronic",
                 "sense_overlapping")

# create list of intersection, unique, etc. -------------------------------------------------------
#NOTE: this is restricted to above lncRNA biotypes
lncipediaIntersection <- gencode29.genes[gencode29.genes$gene_id %in% geneIDs.94$ensembl94GeneID,] #13206 lncipedia ENSG IDs in gencode 29
lncipediaIntersection <- lncipediaIntersection[lncipediaIntersection$gene_type %in% lncBiotypes,] #13206: same number 

gencodeUnique <- gencode29.genes[!gencode29.genes$gene_id %in% lncipediaIntersection$gene_id,] #gencode genes (matching lncRNA biotypes) not in lncipedia

lncipediaTable <- table(lncipediaIntersection$gene_type) #biotype of lncRNA in intersection of two annotations

gencodeTable <- table(gencode29.genes$gene_type)[lncBiotypes]

gencodeKeep <- gencodeUnique

#same list for transcripts
gencodeTranscriptTable <- table(subset(gencode29.df, type == "transcript")$transcript_type)
lncipediaTranscriptIntersection <- subset(gencode29.df, type == "transcript")[gencode29.genes$transcript_id %in% geneIDs.94$ensembl94GeneID,] #13206 lncipedia ENSG IDs in gencode 29
lncipediaTranscriptIntersection <- lncipediaIntersection[lncipediaIntersection$transcript_type %in% lncBiotypes,] #13206: same number 

# make venn diagram highlighting changes -----------------------------------------------------------
dir.create("img")

tmp <- draw.pairwise.venn(sum(gencodeTable),
                          length(unique(lncipedia5.2$gene_id)),
                          sum(lncipediaTable),
                          category = c("GENCODE 29", "lncipedia 5.2"),
                          fill = c("#00AFBB", "#E7B800"),
                          alpha = rep(0.5, 2), cat.pos = c(-105, 135), cat.dist = rep(0.08, 2))

png("img/lncRNA_gene_venn.png")
grid.draw(tmp)
dev.off()

tmp2 <- draw.pairwise.venn(sum(gencodeTranscriptTable),
                          length(unique(lncipedia5.2$transcript_id)),
                          sum(lncipediaTable),
                          category = c("GENCODE 29", "lncipedia 5.2"),
                          fill = c("#00AFBB", "#E7B800"),
                          alpha = rep(0.5, 2), cat.pos = c(-105, 135), cat.dist = rep(0.08, 2))

png("img/lncRNA_transcript_venn.png")
grid.draw(tmp)
dev.off()

#create gencode gtf with lncipedia genes removed and export -----------------------------------------
gtfNew <- gencode29[gencode29$gene_id %in% gencodeUnique$gene_id,]

export(gtfNew, 'gencode.v29.primary_assembly.annotation.lncipediaRemoved.gtf')

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



# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: July 19th, 2022

#Purpose: To create a gene annotation and transcript annotation text files with biotype information, etc.

#Imports: GENCODE40 and LNCipedia5.2 (high confidence) gtf files
# list of ens IDs in intersection
#Exports: 
#gene annotations
#transcript annotations

# import packages -------------------------------------------------------------------------------
library(rtracklayer) #for importing gtf
library(GenomicFeatures) #for creating TxDB
library(tidyverse)
#library(gridExtra)

# import combined annotation ---------------------------------------------
combinedAnno.gtf <- "gencode_v40.lncipedia_v5_2_hc.annotation.gtf"
combinedAnno <- rtracklayer::import(combinedAnno.gtf, format = "gtf")
combinedAnno.df <- as.data.frame(combinedAnno)

# import GENCODE40 annotation -------------------------------------------------------------------
gencode.v40.gtf <- "gencode.v40.primary_assembly.annotation.gtf"
gencode40 <- rtracklayer::import(gencode.v40.gtf, format = "gtf")
gencode40.df <- as.data.frame(gencode40)

gencode40.genes <- subset(gencode40.df, type == "gene") #list genes only in gencode annotation

gencode40.genes.anno <- gencode40.genes %>% dplyr::select(gene_id, gene_type, gene_name) %>%
  transmute(ensembl_id = gene_id, ensembl_gene_name = gene_name, ensembl_gene_type = gene_type)

#make tx2gene ---------------------------------------------
txdb <- makeTxDbFromGFF(combinedAnno.gtf, format = "gtf")
saveDb(x=txdb, file = "gencode_v40.lncipedia_v5_2_hc.annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME") 
head(tx2gene)
dim(tx2gene)
write_csv(tx2gene, "tx2gene.gencode_v40.lncipedia_v5_2_hc.csv")

#make .txt files with annotation features (gene, transcript, biotype, etc) ---------------------------------------------
combinedAnno.genes <- combinedAnno.df %>% dplyr::filter(type == "gene")
transcripts.anno <- combinedAnno.df %>% dplyr::filter(type == "transcript") %>%
  mutate(transcript_name = coalesce(transcript_name, transcript_id)) %>%
  mutate(transcript_type = replace_na(transcript_type, "lncipedia_lncRNA"))

genes.anno.tmp <- combinedAnno.genes %>% dplyr::select(gene_id, gene_type, gene_name) %>%
  mutate(gene_name = coalesce(gene_name, gene_id)) %>%
  mutate(gene_type = replace_na(gene_type, "lncipedia_lncRNA"))

#add in ensembl ID for matching lncipedia IDs -------------------------------------------------------------------
# import 92 to 106 table, remove retired IDs
ens106 <- read.csv("ensembl106GeneUpdate.csv", header = TRUE)
geneIDs.106 <- ens106[grep("retired", ens106$New.stable.ID, invert = TRUE),] %>%
  transmute(lncipediaGeneID = lncipediaGeneID,
            matching_ens_gene_id = gsub(" ", "", New.stable.ID),
            source = "lncipedia_map")

#import intersection/overlapping genes
intersection.map <- read_tsv("intersection.genes.nopc.txt", col_names = FALSE) %>%
  dplyr::select(X9, X18) %>%
  transmute(lncipediaGeneID = sapply(str_split(X18, "\\\""), "[", 2),
            matching_ens_gene_id = sapply(str_split(X9, "\\\""), "[", 2),
            source = "intersection"
            )

conversionTable <- rbind(geneIDs.106, intersection.map)

genes.anno <- left_join(genes.anno.tmp, conversionTable, by = c("gene_id" = "lncipediaGeneID")) %>%
  mutate(ensembl_id = if_else(!gene_type == "lncipedia_lncRNA", gene_id, matching_ens_gene_id)) %>%
  left_join(gencode40.genes.anno)

#tmp <- genes.anno %>% filter(gene_type == "lncipedia_lncRNA")
#table(tmp$ensembl_gene_type)


#add in biotype ----------------------------------------------------------------------------
#  biotype classes
#NOTE: see the following link
#https://uswest.ensembl.org/Help/Faq?id=468
#http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
pcBiotypes <- c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene",
                "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene",
                "IG_D_gene", "nonsense_mediated_decay", "non_stop_decay")
lncBiotypes <- c("lncRNA", "lncipedia_lncRNA", "processed_transcript", "retained_intron")
short.ncBiotypes <- c("snRNA", "snoRNA", "misc_RNA", "miRNA", "rRNA",
                      "Mt_tRNA", "Mt_rRNA", "ribozyme", "scaRNA",
                      "sRNA", "scRNA", "vault_RNA")
pseudogeneBiotypes <- c("processed_pseudogene", "transcribed_unprocessed_pseudogene",
                        "unprocessed_pseudogene", "transcribed_processed_pseudogene",
                        "transcribed_unitary_pseudogene", "rRNA_pseudogene",
                        "unitary_pseudogene", "polymorphic_pseudogene",
                        "pseudogene", "IG_V_pseudogene",
                        "translated_processed_pseudogene", "translated_unprocessed_pseudogene",
                        "TR_V_pseudogene",
                        "IG_C_pseudogene", "TR_J_pseudogene",
                        "IG_J_pseudogene", "IG_pseudogene")
TEC.Biotypes <- c("TEC")

genes.anno <- genes.anno %>% mutate(biotype_class = "unknown")
genes.anno$biotype_class[genes.anno$gene_type %in% pcBiotypes] <- "protein_coding"
genes.anno$biotype_class[genes.anno$gene_type %in% lncBiotypes] <- "short_noncoding"
genes.anno$biotype_class[genes.anno$gene_type %in% short.ncBiotypes] <- "long_noncoding"
genes.anno$biotype_class[genes.anno$gene_type %in% pseudogeneBiotypes] <- "pseudogene"
genes.anno$biotype_class[genes.anno$gene_type %in% TEC.Biotypes] <- "TEC"
unique(genes.anno$biotype_class) #checking

transcripts.anno <- transcripts.anno %>% mutate(biotype_class = "unknown")
transcripts.anno$biotype_class[transcripts.anno$transcript_type %in% pcBiotypes] <- "protein_coding"
transcripts.anno$biotype_class[transcripts.anno$transcript_type %in% lncBiotypes] <- "short_noncoding"
transcripts.anno$biotype_class[transcripts.anno$transcript_type %in% short.ncBiotypes] <- "long_noncoding"
transcripts.anno$biotype_class[transcripts.anno$transcript_type %in% pseudogeneBiotypes] <- "pseudogene"
transcripts.anno$biotype_class[transcripts.anno$transcript_type %in% TEC.Biotypes] <- "TEC"
unique(transcripts.anno$biotype_class) #checking

#output annotation description files ----------------------------------------------------------------
write_fst(genes.anno, "GRCh38_GENCODE40LNCipedia_geneAnnotations.fst")
write_fst(transcripts.anno, "GRCh38_GENCODE40LNCipedia_transcriptInfo.fst")


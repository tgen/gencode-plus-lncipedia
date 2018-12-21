#!/bin/bash
#
# Author: Elizabeth Hutchins
# Date: Nov. 1st, 2018
#
# Create nonredundant annotation combining GENCODE29 and LNCipedia 5.2 (high confidence)

#download gtf, fa file: GENCODE29
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
gunzip gencode.v29.primary_assembly.annotation.gtf.gz

#download gtf,fa file: lncipedia 5.2
wget https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.gtf
wget https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc.fasta

#download lncipedia 5.2/ensembl 92 conversion table
wget https://lncipedia.org/downloads/lncipedia_5_2/lncipedia_5_2_ensembl_92_genes.txt

#add quotes around gene names for compatibility
sed 's/_id /_id \"/g' lncipedia_5_2_hc_hg38.gtf > lncipedia_5_2_hc_hg38.formatted.gtf

sed -i 's/ ;/" ;/g' lncipedia_5_2_hc_hg38.formatted.gtf

for f in `seq 1 16`
do sed -i "s/_alias_${f} /_alias_${f} \"/g" lncipedia_5_2_hc_hg38.formatted.gtf
done

#remove headers
sed -i 1,2d lncipedia_5_2_hc_hg38.formatted.gtf
sed -i 1,5d gencode.v29.primary_assembly.annotation.gtf

#run Rscript to remove redundant lncRNAs from GENCODE and
#add gene and transcript rows to LNCipedia
Rscript incorporateLNCipedia.R

#extract only gene rows (exon rows will come from transcript file)
awk '{if($3=="gene") print}' lncipedia_5_2_hc_hg38.genes.gtf > lncipedia_5_2_hc_hg38.genesOnly.gtf

#combine annotation files
cat gencode.v29.primary_assembly.annotation.lncipediaRemoved.gtf lncipedia_5_2_hc_hg38.transcripts.gtf lncipedia_5_2_hc_hg38.genesOnly.gtf | grep -v "#" > gencode_v29.lncipedia_v5_2_hc.annotation.tmp
sort -k1,2n gencode_v29.lncipedia_v5_2_hc.annotation.tmp > gencode_v29.lncipedia_v5_2_hc.annotation.gtf
rm gencode_v29.lncipedia_v5_2_hc.annotation.tmp

#add headers
sed -i '1i##date 2018-11-14' gencode_v29.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##source-version rtracklayer 1.38.3' gencode_v29.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##format: gff-version 2' gencode_v29.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##contact: ehutchins@tgen.org' gencode_v29.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##provider: GENCODE, LNCipedia' gencode_v29.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##description: evidence-based annotation of the human genome (GRCh38), version 29 (Ensembl 94) with added genes from LNCipedia, version 5.2' gencode_v29.lncipedia_v5_2_hc.annotation.gtf


#create transcriptome fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

gffread -w gencode_v29.lncipedia_v5_2_hc.annotation.fa -g GRCh38.primary_assembly.genome.fa gencode_v29.lncipedia_v5_2_hc.annotation.gtf

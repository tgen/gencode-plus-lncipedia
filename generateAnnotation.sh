#!/bin/bash
#
# Author: Elizabeth Hutchins
# Date: Nov. 1st, 2018
# Update: June 2022 to incorporate GENCODE40
# Requirements: bedtools, gffread
#
# Create nonredundant annotation combining GENCODE40 and LNCipedia 5.2 (high confidence)

#download lncipedia 5.2/ensembl 92 conversion table
wget --no-check-certificate https://lncipedia.org/downloads/lncipedia_5_2/lncipedia_5_2_ensembl_92_genes.txt

# Run python script to generate ensembl92 to ensembl104 mappings
python3 ensembl_map_92_104.py

#make directory for annotation files
mkdir makeanno
cd makeanno

#download gtf, fa file: GENCODE40
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gtf.gz
gunzip gencode.v40.primary_assembly.annotation.gtf.gz

#download gtf,fa file: lncipedia 5.2
wget --no-check-certificate https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.gtf
wget --no-check-certificate https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc.fasta

#add quotes around gene names for compatibility
sed 's/_id /_id \"/g' lncipedia_5_2_hc_hg38.gtf > lncipedia_5_2_hc_hg38.formatted.gtf

sed -i 's/ ;/" ;/g' lncipedia_5_2_hc_hg38.formatted.gtf

for f in `seq 1 16`
do sed -i "s/_alias_${f} /_alias_${f} \"/g" lncipedia_5_2_hc_hg38.formatted.gtf
done

#remove headers
sed -i 1,2d lncipedia_5_2_hc_hg38.formatted.gtf
sed -i 1,5d gencode.v40.primary_assembly.annotation.gtf

#run Rscript to remove redundant lncRNAs from GENCODE and
#add gene and transcript rows to LNCipedia
Rscript incorporateLNCipedia.R

#create genes only file for gencode for intersectbed
#use intersectbed to find overlapping genes (reciprocal 70%)
#between the two annotations that aren't in the lncipedia conversion table
grep 'ENSEMBL  gene' gencode.v40.primary_assembly.annotation.lncipediaRemoved.gtf >  gencode.v40.primary_assembly.annotation.lncipediaRemoved.ENSEMBLgenes.gtf
grep 'HAVANA   gene' gencode.v40.primary_assembly.annotation.lncipediaRemoved.gtf >  gencode.v40.primary_assembly.annotation.lncipediaRemoved.HAVANAgenes.gtf
cat gencode.v40.primary_assembly.annotation.lncipediaRemoved.ENSEMBLgenes.gtf gencode.v40.primary_assembly.annotation.lncipediaRemoved.HAVANAgenes.gtf > gencode.v40.primary_assembly.annotation.lncipediaRemoved.genes.gtf

intersectBed -nonamecheck -wa -s -r -f 0.7 -a gencode.v40.primary_assembly.annotation.lncipediaRemoved.genes.gtf -b lncipedia_5_2_hc_hg38.genes.gtf > intersection.genes.gtf

#make list of ensembl IDs that intersect an remove from gencode annotation file
awk -F"\t" '{print $9}' intersection.genes.gtf | awk -F';' '{print $1}' | sed 's/gene_id //g' | sed 's/"//g' | sort | uniq > intersection.IDs.txt
grep -v -f intersection.IDs.txt gencode.v40.primary_assembly.annotation.lncipediaRemoved.gtf > gencode.v40.primary_assembly.annotation.lncipediaRemoved.nrd.gtf

#extract only gene rows (exon rows will come from transcript file)
awk '{if($3=="gene") print}' lncipedia_5_2_hc_hg38.genes.gtf > lncipedia_5_2_hc_hg38.genesOnly.gtf

#combine annotation files
cat gencode.v40.primary_assembly.annotation.lncipediaRemoved.nrd.gtf lncipedia_5_2_hc_hg38.transcripts.gtf lncipedia_5_2_hc_hg38.genesOnly.gtf | grep -v "#" > gencode_v40.lncipedia_v5_2_hc.annotation.tmp
sort -k1,2n gencode_v40.lncipedia_v5_2_hc.annotation.tmp > gencode_v40.lncipedia_v5_2_hc.annotation.gtf
rm gencode_v40.lncipedia_v5_2_hc.annotation.tmp

#add headers
sed -i '1i##date 2022-06-28' gencode_v40.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##source-version rtracklayer 1.38.3' gencode_v40.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##format: gff-version 2' gencode_v40.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##contact: ehutchins@tgen.org' gencode_v40.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##provider: GENCODE, LNCipedia' gencode_v40.lncipedia_v5_2_hc.annotation.gtf
sed -i '1i##description: evidence-based annotation of the human genome (GRCh38), version 40 (Ensembl 104) with added genes from LNCipedia, version 5.2' gencode_v40.lncipedia_v5_2_hc.annotation.gtf


#create transcriptome fasta
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

gffread -w gencode_v40.lncipedia_v5_2_hc.annotation.fa -g GRCh38.primary_assembly.genome.fa gencode_v40.lncipedia_v5_2_hc.annotation.gtf


#make decoys and salmon gentrome for indexing
#see: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
rm decoys.txt.bak

cat gencode_v40.lncipedia_v5_2_hc.annotation.fa GRCh38.primary_assembly.genome.fa > gentrome.fa
gzip gentrome.fa
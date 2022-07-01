Combines [GENCODE40](https://www.gencodegenes.org/human/release_40.html) and [lncipedia5.2](https://lncipedia.org/) gene annotations in a non-redundant manner.

Can be run with the following command:

```
bash generateAnnotation.sh
```

Dependencies:
 - [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml)

R packages:
 - [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
 - [tidyverse](https://www.tidyverse.org/)
 - [VennDiagram](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-35)



Note 1: gtf can be converted to bed with eautils
https://github.com/ExpressionAnalysis/ea-utils

```
ea-utils/clipper/gtf2bed gencode_v40.lncipedia_v5_2_hc.annotation.gtf > gencode_v40.lncipedia_v5_2_hc.annotation.bed12
```

or with ucsc tools
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
```
gtfToGenePred -genePredExt -geneNameAsName2 gencode_v40.lncipedia_v5_2_hc.annotation.gtf gene.tmp
genePredToBed gene.tmp gencode_v40.lncipedia_v5_2_hc.annotation.bed12
```


Note 2: salmon index can be created like so

```
salmon index -t gentrome.fa.gz -d decoys.txt -p 20 --gencode -i salmon_index/
```

Beñat Yañez
30/11/2023

------------------------------------------------------------------------

Main code used throught the Dissertation project.
The Dissertation was carried out for the Evolution course in the Ecology, Ecology, Evolution and Biodiversity MSc
at Edinburgh University.
Supervised by Simon Martin
------------------------------------------------------------------------

The bash (.sh) and python (.py) codes are used the Edinburgh University cluster system. 
The data files required to run them can be found in (/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data).
Much of the analysis uses code developed by Simon Martin found in (https://github.com/simonhmartin/genomics_general)
The R code all the data to be able to run it should be in this repository in the data folder.

------------------------------------------------------------------------

The results folder contains the ouput of the R code and the figures as pdfs.
Summary of each Script:

0fold_4fold_diversity3supergenes.sh
 - Run in the cluster
 - Just for the supergene region: Chromosome 15 from 5322257 to 6220875 and 6251636 to 7829094
 - Calculate 4fold and 0fold diversity at 100kb windows
 - Different diversity estimates for Mediterranean individuals that were Homozygote Chrysippus, Homozygote Meidterranean and Heterozygotes for the supergene

0fold_4fold_diversitygenome.sh
 - Run on cluster
 - Calculate 4fold and 0fold diversity at 100kb window
 - For the entire genome (One chromosome at a time)
 - Different diversity estimates for Mediterranean and African individuals

0fold_4fold_diversitysupergene.sh
 - Run on cluster
 - Calculate 4fold and 0fold diversity at 100kb window
 - Just for the supergene region: Chromosome 15 from 5322257 to 6220875 and 6251636 to 7829094
 - Different diversity estimates for Mediterranean and African individuals

0fold_4fold_genomes.sh
 - Run on cluster
 - Use the 0fold .bed file and 4fold .bed file to only select those sites from the .geno file
 - For the entire genome (One chromosome at a time)

0fold_4fold_ratio.R
 - Run on cluster
 - Take output of 0fold_4fold_diversitygenome.sh and calculate the ratio of 0 fold to 4 fold diversity
 - Modify imput to do the same for 0fold_4fold_diversity3supergenes.sh and 0fold_4fold_diversitysupergene.sh outputs

0fold_4fold_sites.sh
 - Run on cluster
 - Using Annotation file (.gff3), sequence (.fasta) and vcf file, identifies the 0 fold (any base change, modifies the amino acid) and 4 fold (any base change does not modify amino acid) and makes a .bed file of them

0fold_4fold_supergene.sh
 - Run on cluster
 - Use the 0fold .bed file and 4fold .bed file to only select those sites from the .geno file
 - Just for the supergene region: Chromosome 15 from 5322257 to 6220875 and 6251636 to 7829094

4DSitesdiversity.sh
 - Run on cluster
 - Uses a .geno file of all the 4 fold sites in the entire genome
 - Used if 0fold_4fold_diversitygenome.sh did not work

4DSitesdiversity_Supergene.sh
 - Run on cluster
 - Uses a .geno file of all the 4 fold sites in the supergene
 - Used if 0fold_4fold_diversitysupergene.sh did not work



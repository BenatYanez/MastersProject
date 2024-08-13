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

DiagnosticSNPs.sh
 - Run on cluster
 - Makes a subset of the SNPS from the PCR amplified region: Chromosome 15 from 6301177 to 6302042 
 - These sites are used in HardyWeinberg.R to identify diagnostic SNPs that can differentiate between 3 alleles (Chrysippus, Orientis-like and Mediterranean)

HardyWeinberg.R
 - Run natively
 - Using the genotypes from DiagnosticSNPs.sh, transforms them into bases (A,T,C,G)
 - Plots these bases along the DNA sequence
 - Using genotype frequencies from genome data and Sanger sequencing compare observed genotypes to expectations under Hardy-Weinberg using Chi-Squared test
 - Do this for the supergene (Two Alleles: Chrysippus homozygote, Mediterranean homozygote, Heterozygote)
 - Do this for the B locus (Three Alleles: Chrysippus homozygote, Mediterranean homozygote, Orientis-like homozygote, Chrysippus+Mediterranean heterozygote, Chrysippus+Orientis-like heterozygote, Mediterranean+Orientis-like heterozygote)

ImageTemperature.R
 - Run natively
 - Uses the Thermal Images in the data folder
 - Obtains the Temperature of the thorax from greyscale Thermal Images
 - Creates a dataframe with a time series of the temperature for these images, all runs are stored in the same dataframe

MakeBedFiles.R
 - Run on the cluster
 - Takes the bed files from 0fold_4fold_sites.sh and removes sites within the supergene to create a bed file for the genome
 - Takes the bed files from 0fold_4fold_sites.sh and only selects sites within the supergene to create a bed file for the supergene

PCA.R
 - Run natively
 - Uses the output from PCAInversionTypes.sh
 - Plots the first and second principal components against each other to identify Chrysippus homozygotes, Mediterranean Homozygotes and Heterozygotes for the supergene

PCAInversionTypes.py
 - Run on the cluster
 - Developed by Thomas Decroly
 - Runs a Principal Compenent Analysis on the vcfs of the supergene region of the Mediterranean population

PCAInversionTypes.sh
 - Run on the cluster
 - Runs PCAInversionTypes.py on the cluster

PlottingRatios.R
 - Run natively
 - Uses the diversity Ratios in the data folder created by 0fold_4fold_ratio.R
 - Plots the resulting ratio as a boxplot, separating results by geographic region (African, Mediterranean), genetic region (genome, supergene) and Sueprgene arrangement (Chrysippus homozygotes, Mediterranean Homozygotes and Heterozygotes)
 - Uses Welch t-tests to determine whether there are statistical differences between these

PowerAnalysis.R
 - Run natively
 - Part of the pilot study (During the Project Proposal phase)
 - Samples groups of different sizes under different heterozygote advantages and compares these to expectations under Hardy-Weinber to determine the power to detect Heterozygote advantage

SNPEffAnalysis.R
 - Run natively
 - Uses the output of UnifySnpEff.R
 - Calcuates the ratio of HIGH impact variants to callable sites in 100kb windows for genome, supergene for the Mediterranean and African population and the three supergene arrangements (Chrysippus homozygotes, Mediterranean Homozygotes and Heterozygotes)
 - Plots these as a boxplot
 - Uses Welch t-tests to determine whether there are statistical differences between these

StatisticalModels.R
- Run natively
- Uses the datframe made by ImageTemperature.R to analyse the results statistically
- Runs a quadratic and log transformed mixed effect model accounting for temporal autocorrelation using glmmTMB
- Does the same for a modified dataset where the starting point is temperature of 21ºC
- Plots the results

UnifySnpEff.R
 - Run on the cluster
 - Combines the files made by how_high.sh (One per chromosome) into a single file
 - Does it for both the supergene and Genome

how_high.py
 - Run on the cluster
 - Developed by Thomas Decroly
 - Uses SNPEff annotated vcfs to return the allele frequency of high impact variants
 - Ignores singletons
   
how_high.sh
 - Run on the cluster
 - Runs how_high.sh on the cluster for 100kb windows
 - Uses an annotated vcf (obtained by snpEff.sh) , an annotation file (.gff3) and a csv with the names of the samples being analysed

snpEff.sh
 - Run on the cluster
 - Uses SnpEff to annotate the genome vcf of each chromosome

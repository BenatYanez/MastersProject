------------------------------------------------------------------------

This is the  data used for the analysis that can be carried out natively
The files needed for analysis in the cluster are too large and should be requested from Simon Martin

------------------------------------------------------------------------

Explanation of the data

/SangerSequences
  - .ab1 files sent by Azenta Life Sciences.
  - They can be opened in SnapGene Viewer to visualize the chromatogram
  - Homozygotes will only have one peak per timepoint while for heterozygotes peaks will overlap
  - Used to genotype individuals
  - File name indicates which individual ID it belongs to
  - ID in the form of AA11BB222 where AA is the initial of who run the field trip to capture the individuals, 11 is the year the trip took place, BB is the site id (ie, AN is Andalusia), and 222 is a unique ID for each individual captured in each field trip (ie 001,002,003)

/ThermalImages
  - Contains several folder, each folder is for the day the temperature measurement experiments were done (3/06/2024-6/06/2024)
  - Each of these folders contains several folder that are numbered, the numbers correspond to a unique butterfly captured in Andalusia (ie, the 222 part of the ID name)
  - If the same butterfly was measured twice the same day the folder is names as X_2
  - Each folder contains jpgs in greyscale taken by a DT-980Y Thermal Imager

/diversity
  - .4D and .0D for 4fold and 0fold diversity estimates respectively
  - Data made by 0fold_4fold_diversitygenome.sh, 0fold_4fold_diversitysupergene.sh and 0fold_4fold_diversity3supergenes.sh
  - One file per chromosome 1-31 (Chromosome 1 and 21 do not exist), and the inversion (inversion.w100kb.m1kb.*D) 
  - One file of 4fold diversity for the entire genome (Autosome.w100kb.m1kb.4D)
  - File were diversity estimates are calculted for individuals belonging to different supergene arrangements rather than for all individuals (Separate.inversions.w100kb.m1kb.*D)
  - 1st collumn identifies the chromosome for that window. 2nd and 3rd collumn idenify the first and last site of the 100kb window respectively
  - 4th column the site at which the middle SNP is found (if there is 200 SNPs in a window, the site of the 100th SNP). 5th column is the total number of 4fold or 0fold SNPs within the window
  - The nect are the diversity estimates for the different populations (ie Mediterranean, African or 3 different supergene Arrangements)

Genotypes_DanausChrysippus.csv
  - Genotype of 93 individuals (from genome data, those catured in Andalusia and some captured previously)
  - Collums indicating the ID in the form of AA11BB222 where AA is the initial of who run the field trip to capture the individuals, 11 is the year the trip took place, BB is the site id (ie, AN is Andalusia), and 222 is a unique ID for each individual captured in each field trip (ie 001,002,003)
  - The Second column is the Genotype (Homozygote or Heterozygote)
  - The third and fourth column is which allele makes the genotype (Orientis-like, Mediterranean and Chrysippus). If the sample is homozygote only the third collumn is filled.

HighImpact100kbWindowsGenome/Supergene.tsv
  - Files obtained by how_high.sh
  - Show the number of HIGH impact sites in a 100kb window for all individuals in a sample (Mediterranean, African and 3 supergene arrangements) and the number of possible callable sites
- 1st collumn identifies the chromosome for that window. 2nd and 3rd collumn idenify the first and last site of the 100kb window respectively
- 4th collumn the number of HIGH impact sites in the Window
- 5th collumn the number of callable sites in the Window. (Some windows will not have any callable or HIGH impact sites)
- 6th collumn the population the results for that window belong to

Separate_InversionHet_Ratios_100kb/Separate_InversionHom1_Ratios_100kb/Separate_InversionHom2_Ratios_100kb/Inversion_Ratios_100kb/ZeroFour_Ratios_100kb.txt
  - Obtained with the 0fold_4fold_ratio.R script
  - The 0fold to 4fold ratios for 100kb windows for Heterozygotes , Chrysippus homozygotes and Mediterranean homozygotes for the supergene
  - Also The ratios for 100kb windows in the genome (ZeroFour_Ratios_100kb) and the Supergene (Inversion_Ratios_100kb), for Mediterranean and African samples
  - Same colums as /diversity but last columns represent the ratio rather than just diversity

Spain_samples.csv
  - Information on phenotype of the individuals captured in Andalusia the last week of May 2024.
  - 1st collumn ID in the form of AA11BB222 where AA is the initial of who run the field trip to capture the individuals, 11 is the year the trip took place, BB is the site id (ie, AN is Andalusia), and 222 is a unique ID for each individual captured in each field trip (ie 001,002,003).
  - 2nd collumn the sex of the butterfly (M=Male, F=Female)
  - 3rd collumn the morph of the forewing pale or dark
  - 4th and 5th collumn Colour of the thorax and abdomen (pale, dark, intermediate)
  - 6th collumn within Andalusia which town we captured the individuals in (Salobreña or Motril)
  - 7th collumn how deteriorated the wings were when we captured the individuals (Fresh=No deteriaration, Faded and Very Faded)
  - 8th collumn mean Length of both forewings from the middle of proximal end to the M2 vein end
  - 9th collumn  mean Height of both forewings from the end of 2A vein to end of R2 vein

dchry2.3.chr15.PCRAmplification.3Genotypes.vcf
  - vcf file for the genome samples
  - The region in this subset corresponds to the one amplified using PCR on the captured individuals
  - USed for genotyping


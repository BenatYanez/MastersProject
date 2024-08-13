#!/bin/bash

#Beñat Yañez
#From geno files (used by popgenWindows.py to obtain diversity) of each chromosome using the bed file with the 0fold and 4 fold sites (obtained with 0fold_4fold_sites.sh script),
#subset of the geno file to only include 0fold or 4fold sites (For the supergene)
#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N Degeneracy_inversion_subset                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h=c3
#$ -pe smp 1                                                                                  # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print
#Generate the inversion (supergene) geno file from a vcf
#python genomics_general/VCF_processing/parseVCF.py -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/vcfs/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.vcf.gz | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.Inversion.GQ20.DP8.geno.gz
#0fold
python genomics_general/mergeGeno.py -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.Inversion.GQ20.DP8.geno.gz -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.0Dsites.bed -f /data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.3/Dchry2.3.fa.fai --outputOnly 1 | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.0DSites.geno.gz
#4 fold
python genomics_general/mergeGeno.py -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.Inversion.GQ20.DP8.geno.gz -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.4Dsites.bed -f /data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.3/Dchry2.3.fa.fai --outputOnly 1 | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry.2.3.chr15.GQ20.DP8.Inversion.4DSites.geno.gz

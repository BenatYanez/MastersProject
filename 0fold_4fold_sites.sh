#!/bin/bash

#Beñat Yañez
#Obtain the degeneracy of each codon in the Danaus Chrysippus genome using the codingSiteTypes.py script by Simon Martin (https://github.com/simonhmartin/genomics_general?tab=readme-ov-file)
#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N Degeneracy                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h=c3
#$ -pe smp 2                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print

#Dchry2.3
python genomics_general/codingSiteTypes.py -a /data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.3/Dchry2.3.masked_annotation_transferred_from_Dchry2.2.tidy.gff3 -f gff3 -r /data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.3/Dchry2.3.fa --ignoreConflicts | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.coding_site_types.tsv.gz
#Obtain just the 4 fold and 0fold sites
gunzip -c /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.coding_site_types.tsv.gz | awk 'BEGIN {OFS="\t"}; (NR==1 || $5=="4") {print($1,$2)}' > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.4Dsites.bed
gunzip -c /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.coding_site_types.tsv.gz | awk 'BEGIN {OFS="\t"}; (NR==1 || $5=="0") {print($1,$2)}' > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/output.0Dsites.bed

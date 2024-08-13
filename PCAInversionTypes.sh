#!/bin/bash
#Run PCAInversionTypes.py on the cluster to obtain a PCA of the Supergene region for individuals of the Mediterranean population
#PCAInversionTypes.py developed by Thomas Decroly (https://github.com/thomdec/scripts)
# Grid Engine options (lines prefixed with #$)t
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h=c1
#$ -pe smp 8                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR prin
#$ -N PCA

python PCAInversionTypes.py

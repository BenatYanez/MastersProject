#!/bin/bash

#Beñat Yañez
#Obtain 0fold and 4 fold diversity in the supergene(in chromosome 15), for Mediterranean and African samples
#Uses the popgenWindows.py script to calculate diversity, made by Simon Martin (https://github.com/simonhmartin/genomics_general?tab=readme-ov-file)
#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N Diversity0_4Fold_Inversion                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h="c1|c2|c3"
#$ -pe smp 1                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print
Wsize=100
#0 fold
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Mediterranian -p African --popsFile Samples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.0DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/inversion.w"$Wsize"kb.m1kb.0D -f phased --writeFailedWindows
#4Fold
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Mediterranian -p African --popsFile Samples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry.2.3.chr15.GQ20.DP8.Inversion.4DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/inversion.w"$Wsize"kb.m1kb.4D -f phased --writeFailedWindows


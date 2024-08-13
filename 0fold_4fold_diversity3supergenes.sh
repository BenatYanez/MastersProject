#!/bin/bash

#Beñat Yañez
#Obtain the diversity of 0fold and 4fold sites for Mediterranean Homozygotes, Chrysippus homozygotes and the Heterozygotes
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
#0fold sites
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Hom1 -p Hom2 -p Het -p African --popsFile InverSamples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.0DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/Separate.inversions.w"$Wsize"kb.m1kb.0D -f phased --writeFailedWindows
#4Fold
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Hom1 -p Hom2 -p Het -p African --popsFile InverSamples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry.2.3.chr15.GQ20.DP8.Inversion.4DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/Separate.inversions.w"$Wsize"kb.m1kb.4D -f phased --writeFailedWindows

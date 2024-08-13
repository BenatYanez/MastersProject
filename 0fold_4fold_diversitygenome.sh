#!/bin/bash

#Beñat Yañez
#Obtain 0fold and 4 fold diversity in the genome(separated into chromosomes), for Mediterranean and African samples
#Uses the popgenWindows.py script to calculate diversity, made by Simon Martin (https://github.com/simonhmartin/genomics_general?tab=readme-ov-file)
#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N Diversity0_4Fold                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h="c1|c2|c3"
#$ -t 1-31
#$ -tc 15
#$ -pe smp 1                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print
num=$SGE_TASK_ID
Wsize=100
#0fold
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Mediterranian -p African --popsFile Samples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr"$num".GQ20.DP8.0DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/chr"$num".w"$Wsize"kb.m1kb.0D -f phased --writeFailedWindows
#4fold
python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Mediterranian -p African --popsFile Samples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr"$num".GQ20.DP8.4DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/chr"$num".w"$Wsize"kb.m1kb.4D -f phased --writeFailedWindows

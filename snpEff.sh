#!/bin/bash

#Beñat Yañez
#Use SNPEff to annotate the vcfs for each chromosome and the Supergene Region
#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N SnpeEff                                                                   # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -t 1-31
#$ -tc 10
#$ -l h="c1|c2|c3"
#$ -pe smp 1                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR prin
num=$SGE_TASK_ID
contig=$(awk -v lineid=$SGE_TASK_ID 'NR==lineid{print;exit}' contigs.txt)

#Supergene
#java -jar /ceph/users/byanez/.conda/envs/MastersProject/share/snpeff-5.2-0/snpEff.jar -c /ceph/users/byanez/.conda/envs/MastersProject/share/snpeff-5.2-0/snpEff.config Dchry2.3 /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/vcfs/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.vcf.gz > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SNPEFFEanotation.eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.ann.vcf

#Genome
java -jar /ceph/users/byanez/.conda/envs/MastersProject/share/snpeff-5.2-0/snpEff.jar -c /ceph/users/byanez/.conda/envs/MastersProject/share/snpeff-5.2-0/snpEff.config Dchry2.3 -noStats /data/martin/genomics/analyses/Danaus_popgen/europe/2_variants/filt_chromosomes/eu40.vs.dchry2.3.chr"$contig".GQ20.DP8.vcf.gz | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/SNPEFFEanotation.eu40.vs.dchry2.3.chr"$num".GQ20.DP8.ann.vcf.gz

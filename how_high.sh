##Need to do sconda MastersProject
#Run the how_high.py script developed by Thomas Decroly to obtain the number of High impact sites and Callable sites in windows of 100kb for the Supergene region for the 3 arrangements
#To look at the genome and compare African to Mediterranean samples instead substitute SnpEffSamplesSupergene.csv with SnpEffSamples.csv and change output name


# Grid Engine options (lines prefixed with #$)
#$ -N HowHIGH                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h="c1|c2|c3"
#$ -t 1-31
#$ -tc 10
#$ -pe smp 1                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print
num=$SGE_TASK_ID
python how_high.py -v /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/SNPEFFEanotation.eu40.vs.dchry2.3.chr"$num".GQ20.DP8.ann.vcf.gz -a /data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.3/Dchry2.3.masked_annotation_transferred_from_Dchry2.2.tidy.gff3 -p SnpEffSamplesSupergene.csv -w 100000 -o HighImpactRatio.100kbWindowsSupergene.chr"$num".tsv

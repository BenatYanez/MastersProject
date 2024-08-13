#Beñat Yañez
#If 0fold_4fold_supergene.sh does not for 4 fold sites
#From a vcf only containing  4D sites subset only the supergene region, transform this into a geno file using parseVCF.py and
#Obtain diversity of 4 fold sites from the supergene(inversion) using popgenWindows.py script by Simon Martin (https://github.com/simonhmartin/genomics_general?tab=readme-ov-file)

#Need to run "sconda MastersProject"

# Grid Engine options (lines prefixed with #$)
#$ -N FourFoldDiversity                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h="c1|c2|c3"
#$ -t 1
#$ -pe smp 5                                                                                   # Run array job on this sub-server
#$ -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/output/          # Folder for STDOUT print
#$ -e /data/martin/genomics/analyses/Danaus_popgen/Benat_project/scripts/error/           # Folder for STDERR print
Wsize=$1

bcftools view /data/martin/genomics/analyses/Danaus_popgen/europe/2_variants/eu40.vs.dchry2.3.autosome.DP8_100.GQ20.4dsites.vcf.gz -r chr15:5322257-6220875,chr15:6251636-7829094 -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/vcfs/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.4DSites.vcf.gz
python genomics_general/VCF_processing/parseVCF.py -i /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/vcfs/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.4DSites.vcf.gz | bgzip > /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.4DSites.geno.gz

python /ceph/users/smartin/Research/genomics_general/popgenWindows.py --windType coordinate --windSize "$Wsize"000 --minSites 1000 -p Mediterranian -p African --popsFile Samples.txt -g /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/geno/eu40.vs.dchry2.3.chr15.GQ20.DP8.Inversion.4DSites.geno.gz -o /data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/diversity/Inversion.w"${WSize}"kb.m1kb.4D -f phased --writeFailedWindows

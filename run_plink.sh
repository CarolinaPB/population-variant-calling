#!/bin/bash
#SBATCH --time=1:0:0
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --error=plink%j.txt
#SBATCH --job-name=plink
#SBATCH --exclude=fat001,fat002,fat101,fat100
#SBATCH --mem=20000

module load plink/1.9-180913

#FILE=/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/pop-var-calling/concat/allchr.vc.gz
# FILE=/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/pop-var-calling/concat/all_samples.vcf.gz
FILE=/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/pop-var-calling/concat/all_samples.vep.vcf.gz
SAMPLES_TO_KEEP=/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/pop-var-calling/concat/old_samples_plink_keep.txt

plink --vcf $FILE --pca --double-id --out chicken_old_samples --chr-set 38 --allow-extra-chr --keep $SAMPLES_TO_KEEP

# configfile: "config.yaml"
from snakemake.utils import makedirs

pipeline = "population-var-calling" # replace your pipeline's name

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
MAPPING_DIR = config["MAPPING_DIR"]
PREFIX = config["PREFIX"]
SPECIES = config["SPECIES"]

SAMPLES, = glob_wildcards(os.path.join(MAPPING_DIR, "{samples}.sorted.bam"))


localrules: create_bam_list, create_file_log    

rule all:
    input:
        files_log,
        f"{PREFIX}.eigenvec",
        f"results/variant_calling/{PREFIX}.vcf.stats.txt"

        
rule create_bam_list:
    output:
        temp("bam_list.txt")
    message:
        "Rule {rule} processing"
    params:
        bam_dir = MAPPING_DIR
    shell:
        "ls {params.bam_dir}/*.bam > {output}"


rule var_calling_freebayes:
    input:
        ref=ASSEMBLY,
        bam= rules.create_bam_list.output,
    output:
        "results/variant_calling/{prefix}.vcf.gz"
    params:
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        scripts_dir = os.path.join(workflow.basedir, "scripts")
    group:
        'group'
    shell:
        """
module load freebayes bcftools vcflib python/2.7.15 samtools
{params.scripts_dir}/freebayes-parallel.sh <({params.scripts_dir}/fasta_generate_regions.py {input.ref}.fai {params.chunksize}) 2 \
-f {input.ref} \
--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 \
-L {input.bam} | vcffilter -f 'QUAL > 20' {input} | bgzip -c > {output}
tabix -p vcf {output}
        """

rule run_vep:
    input:
        rules.var_calling_freebayes.output
    output:
        vcf = 'results/final_VCF/{prefix}.smoove.square.vep.vcf.gz',
        # warnings = "5_postprocessing/{prefix}.suqare.vep.vcf.gz_warnings.txt",
        summary = "'results/final_VCF/{prefix}.smoove.square.vep.vcf.gz_summary.html"
    message:
        'Rule {rule} processing'
    conda:
        "envs/vep_dependencies.yaml"
    params:
        species = SPECIES
    group:
        'group'
    shell:
        """
module load samtools
/cm/shared/apps/SHARED/ensembl-vep/vep -i {input} \
--format vcf \
--buffer_size 5000 \
--offline \
--dir /lustre/nobackup/SHARED/cache/ \
--species {params.species} \
--vcf \
--force_overwrite \
-o {output.vcf} \
--fork 1 \
--compress_output bgzip \
--canonical \
--gene_phenotype \
--regulatory \
--numbers \
--symbol
        """

rule bcftools_stats:
    input:
        # rules.concat_vcf.output.vcf
        rules.run_vep.output.vcf
    output:
        "results/variant_calling/{prefix}.vcf.stats.txt"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        """
        module load bcftools
        bcftools stats -s - {input} > {output}
        """

rule PCA:
    input:
        # rules.concat_vcf.output.vcf
        rules.run_vep.output.vcf
    output:
        "{prefix}.eigenvec"
    message:
        'Rule {rule} processing'
    params:
        prefix=PREFIX
    group:
        'group'
    shell:
        """
        module load R/3.6.2
        module load plink/1.9-180913

        plink --vcf {input} --pca --double-id --out {params.prefix} --chr-set 38 --allow-extra-chr
        Rscript scripts/PCA.Rscript {params.prefix}.eigenvec
        """
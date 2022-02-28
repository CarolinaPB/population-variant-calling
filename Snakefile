# configfile: "config.yaml"
from snakemake.utils import makedirs

pipeline = "population-var-calling"

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
MAPPING_DIR = config["MAPPING_DIR"]
PREFIX = config["PREFIX"]
SPECIES = config["SPECIES"]
NUM_CHRS = config["NUM_CHRS"]

SAMPLES, = glob_wildcards(os.path.join(MAPPING_DIR, "{samples}.sorted.bam"))


localrules: create_bam_list, create_file_log    

rule all:
    input:
        files_log,
        f"results/PCA/{PREFIX}.pdf",
        f"results/final_VCF/{PREFIX}.vep.vcf.stats"

        
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
        temp("results/variant_calling/{prefix}.vcf.gz"),
        temp("results/variant_calling/{prefix}.vcf.gz.tbi")
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
        vcf = 'results/final_VCF/{prefix}.vep.vcf.gz',
        warnings = "results/final_VCF/{prefix}.vcf.gz_warnings.txt",
        summary = "results/final_VCF/{prefix}.vep.vcf.gz_summary.html"
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

rule index_vcf:
    input:
        rules.run_vep.output.vcf
    output:
        'results/final_VCF/{prefix}.vep.vcf.gz.tbi'
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools
tabix -p vcf {input}
        """

rule bcftools_stats:
    input:
        vcf = rules.run_vep.output.vcf,
        idx = rules.index_vcf.output
    output:
        "results/final_VCF/{prefix}.vep.vcf.stats"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        """
        module load bcftools
        bcftools stats -s - {input.vcf} > {output}
        """


rule PCA:
    input:
        vcf = rules.run_vep.output.vcf,
    output:
        eigenvec = "results/PCA/{prefix}.eigenvec",
        eigenval = "results/PCA/{prefix}.eigenval",
    message:
        'Rule {rule} processing'
    params:
        prefix= os.path.join("results/PCA",PREFIX),
        num_chrs = NUM_CHRS
    group:
        'group'
    shell:
        """
        module load plink/1.9-180913
        plink --vcf {input.vcf} --pca --double-id --out {params.prefix} --chr-set {params.num_chrs} --allow-extra-chr --threads 8
        """

rule plot_PCA:
    input:
        eigenvec = rules.PCA.output.eigenvec,
        eigenval = rules.PCA.output.eigenval,
    output:
        "results/PCA/{prefix}.pdf"
    message:
        'Rule {rule} processing'
    params:
        rscript = os.path.join(workflow.basedir, "scripts/basic_PCA.R")
    group:
        'group'
    shell:
        """
module load R/3.6.2
echo $CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
Rscript {params.rscript} --eigenvec={input.eigenvec} --eigenval={input.eigenval} --output={output}
        """
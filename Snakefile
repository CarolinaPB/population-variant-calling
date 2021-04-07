configfile: "config.yaml"

from glob import glob
from Bio import SeqIO
pipeline = "population-var-calling" # replace your pipeline's name

workdir: config["OUTDIR"]

include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
BED = config["BED"]
MAPPING_DIR = config["MAPPING_DIR"]
PREFIX = config["PREFIX"]

SAMPLES, = glob_wildcards(os.path.join(MAPPING_DIR, "processed_reads/{samples}.sorted.bam"))

CHROMOSOMES_LARGE = []

if not os.path.isfile("small_chrs.txt"):
    with open("small_chrs.txt", "w") as small_chrs_file:
        with open(ASSEMBLY) as assembly:
            for record in SeqIO.FastaIO.FastaIterator(assembly): 
                if record.id.endswith(".1"):
                    small_chrs_file.write(record.id)
                else:
                    CHROMOSOMES_LARGE.append(record.id)
else:
    with open(ASSEMBLY) as assembly:
        for record in SeqIO.FastaIO.FastaIterator(assembly): 
            if record.id.endswith(".1"):
                pass
            else:
                CHROMOSOMES_LARGE.append(record.id)



localrules: create_bed_windows, create_bam_list, create_file_log    

rule all:
    input:
        files_log,
        expand("concat/allchr_{prefix}.vcf.gz",prefix=PREFIX),
        expand("concat/allchr_{prefix}_stats.txt", prefix = PREFIX),
        expand("{prefix}.eigenvec", prefix=PREFIX)
        

rule create_bed_windows:
    input:
        BED
    output:
        bed=temp("create_bed_windows.done") # make temp
    message:
        "Rule {rule} processing"
    params:
        # window_size = 10000000
        window_size = 5000000
    shell:
        """
        module load bedtools
        bedtools makewindows -b {input} -w {params.window_size} | awk '{{print $0 > "bedfiles/region_"$1"-"$2"-"$3".bed"}}'
        touch {output}
        """

# to run small bins in one go and big bins separately
#  bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 < {params.window_size}' | awk '{{print $0}}' >> {output.bed}
# bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 >= {params.window_size}' | awk '{{print $0 > "bedfiles/region_"$1"-"$2"-"$3".bed"}}'


rule create_bam_list:
    output:
        temp("bam_list.txt")
    message:
        "Rule {rule} processing"
    params:
        bam_dir = MAPPING_DIR
    shell:
        "ls {params.bam_dir}processed_reads/*.bam > {output}"


BEDFILES2, = glob_wildcards("bedfiles/{bed}.bed")
rule run_freebayes:
    input:
        regions = "bedfiles/{bed}.bed",
        reference=ASSEMBLY,
        bams= rules.create_bam_list.output,
        # idx_done = "index_RG.done"
    output:
        vcf="vcf/{bed}.vcf.gz"
        # idx=temp("vcf/{bed}.vcf.gz.tbi")
    message:
        "Rule {rule} processing"
    group:
        'group'
    shell:
        """
        module load freebayes vcflib bcftools samtools

        freebayes --use-best-n-alleles 4 \
        -f {input.reference} \
        --bam-list {input.bams} \
        --targets {input.regions} \
        --min-base-quality 20 \
        --min-mapping-quality 30 \
        --min-alternate-fraction 0.2 \
        --haplotype-length 0 | \
        vcffilter -f 'QUAL > 20' | bgzip -c > {output.vcf}

        tabix -p vcf {output.vcf}
        """

#remove rule and change paths of next rules (remove the "sorted" dir)
# rule sort_index_vcf:
#     input:
#         rules.run_freebayes.output.vcf
#     output:
#         vcf = temp("sorted/{bed}.sorted.vcf.gz"),
#         idx = temp("sorted/{bed}.sorted.vcf.gz.tbi")
#     message:
#         "Rule {rule} processing"
#     shell:
#         """
#         module load bcftools samtools
#         bcftools sort -Oz -m 2G {input} > {output.vcf}
#         tabix -p vcf {output.vcf}
#         """

rule concat_vcf:
    input:
        # vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
        vcf = expand("vcf/{bed}.vcf.gz", bed=BEDFILES2)
    output:
        vcf = "concat/allchr_{prefix}.vcf.gz",
        idx = "concat/allchr_{prefix}.vcf.gz.tbi"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        """
        module load bcftools;
        bcftools concat --threads 16 --allow-overlaps --remove-duplicates -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """

rule bcftools_stats:
    input:
        rules.concat_vcf.output.vcf
    output:
        "concat/allchr_{prefix}_stats.txt"
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
        rules.concat_vcf.output.vcf
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
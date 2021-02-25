configfile: "config.yaml"

from glob import glob
from Bio import SeqIO
pipeline = "population-var-calling" # replace your pipeline's name


include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
BED = config["BED"]
MAPPING_DIR = config["MAPPING_DIR"]

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



localrules: create_bed_windows, create_bam_list    

rule all:
    input:
        files_log,
        # "bedfiles/region_small_bins.bed",
        # expand("concat/{chromosome}.vcf.gz", chromosome=CHROMOSOMES_LARGE),
        # "concat/small_chrs.vcf.gz",
        # expand("bcftools_stats/stats_{chr}.txt", chr=CHROMOSOMES_LARGE),
        # "bcftools_stats/stats_small_chrs.txt"
        "concat/allchr.vcf.gz",
        "concat/allchr_stats.txt",

rule create_bed_windows:
    input:
        BED
    output:
        # bed = "bedfiles/region_small_bins.bed",
        bed="create_bed_windows.done"
    message:
        "Rule {rule} processing"
    params:
        # window_size = 10000000
        window_size = 5000000
    shell:
        """
        module load bedtools
        # bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 < {params.window_size}' | awk '{{print $0}}' >> {output.bed}
        # bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 >= {params.window_size}' | awk '{{print $0 > "bedfiles/region_"$1"-"$2"-"$3".bed"}}'
        bedtools makewindows -b {input} -w {params.window_size} | awk '{{print $0 > "bedfiles/region_"$1"-"$2"-"$3".bed"}}'
        touch {output}
        """


rule create_bam_list:
    output:
        "bam_list.txt"
    message:
        "Rule {rule} processing"
    params:
        bam_dir = MAPPING_DIR
    shell:
        "ls {params.bam_dir}/processed_reads/*.bam > {output}"


BEDFILES2, = glob_wildcards("bedfiles/{bed}.bed")
rule run_freebayes:
    input:
        regions = "bedfiles/{bed}.bed",
        reference=ASSEMBLY,
        bams= rules.create_bam_list.output,
        # idx_done = "index_RG.done"
    output:
        "vcf/{bed}.vcf.gz"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load freebayes vcflib bcftools

        freebayes --use-best-n-alleles 4 \
        -f {input.reference} \
        --bam-list {input.bams} \
        --targets {input.regions} \
        --min-base-quality 20 \
        --min-mapping-quality 30 \
        --min-alternate-fraction 0.2 \
        --haplotype-length 0 | \
        vcffilter -f 'QUAL > 20' | bgzip -c > {output}
        """


rule sort_index_vcf:
    input:
        rules.run_freebayes.output
    output:
        vcf = "sorted/{bed}.sorted.vcf.gz",
        idx = "sorted/{bed}.sorted.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load bcftools samtools
        bcftools sort -Oz -m 2G {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# rule concat_vcf_individual:
#     input: 
#         vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
#         # chr_list = create_perchr_vcf("chr_list.txt")
#     output:
#         vcf = "concat/{chromosome}.vcf.gz",
#         idx = "concat/{chromosome}.vcf.gz.tbi"
#     message:
#         "Rule {rule} processing"
#     group:
#         'group'
#     shell:
#         """
#         module load bcftools;
#         bcftools concat --threads 12 --allow-overlaps --remove-duplicates -r {wildcards.chromosome} -Oz -o {output.vcf} {input.vcf}
#         tabix -p vcf {output.vcf}
#         """

# rule concat_vcf_group:
#     input: 
#         vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
#         regions = "small_chrs.txt"
#     output:
#         vcf = "concat/small_chrs.vcf.gz",
#         idx = "concat/small_chrs.vcf.gz.tbi"
#     message:
#         "Rule {rule} processing"
#     group:
#         'group'
#     shell:
#         """
#         module load bcftools;
#         bcftools concat --threads 12 --allow-overlaps --remove-duplicates -R {input.regions} -Oz -o {output.vcf} {input.vcf}
#         tabix -p vcf {output.vcf}
#         """


# rule bcftools_stats:
#     input:
#         "concat/{chr}.vcf.gz"
#     output:
#         "bcftools_stats/stats_{chr}.txt"
#     message:
#         "Rule {rule} processing"
#     group:
#         'group'
#     shell:
#         """
#         module load bcftools
#         sample_list=$(bcftools query -l {input})
#         bcftools stats -S $sample_list > {output}
#         """

rule concat_vcf:
    input:
        vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
    output:
        vcf = "concat/allchr.vcf.gz",
        idx = "concat/allchr.vcf.gz.tbi"
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
        "concat/allchr_stats.txt"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        """
        module load bcftools
        sample_list=$(bcftools query -l {input})
        bcftools stats -S $sample_list > {output}
        """
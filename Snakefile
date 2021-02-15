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
        "bam_list.txt",
        "bedfiles/region_small_bins.bed",
        # "var_calling/vcf.vcf.gz",
        # "concat/complete_vcf.vcf.gz",
        # "vcf/{bed}.vcf.gz"
        # expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam"), samples=SAMPLES),
        # expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam.bai"), samples=SAMPLES),
        # expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG2.bam"), samples =SAMPLES),
        # "chr_list.txt",
        expand("concat/{chromosome}.vcf.gz", chromosome=CHROMOSOMES_LARGE),
        "concat/small_chrs.vcf.gz"

rule create_bed_windows:
    input:
        BED
    output:
        bed = "bedfiles/region_small_bins.bed",
    message:
        "Rule {rule} processing"
    params:
        # window_size = 10000000
        window_size = 5000000
    shell:
        """
        module load bedtools
        bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 < {params.window_size}' | awk '{{print $0}}' >> {output.bed}
        bedtools makewindows -b {input} -w {params.window_size} | awk '$3-$2 >= {params.window_size}' | awk '{{print $0 > "bedfiles/region_"$1"-"$2"-"$3".bed"}}'
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

# rule create_chr_list:
#     input:
#         ASSEMBLY
#     output:
#         "chr_list.txt"
#     message:
#         "Rule {rule} processing "
#     shell:
#         """
#         grep '>' {input} | awk '{{print substr($1,2); }}' > {output}
#         """


# rule add_readgroup_tag:
#     input:
#         os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.bam")
#     output:
#         os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam")
#     message:
#         "Rule {rule} processing"
#     group:
#         "RG"
#     shell:
#         "module load samtools && samtools addreplacerg --threads 12 -r ID:{wildcards.samples} -r SM:{wildcards.samples} -o {output} {input}"


# rule add_readgroup_tag:
#     input:
#         expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.bam"), samples=SAMPLES)
#     output:
#         expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam"), samples =SAMPLES)
#     message:
#         "Rule {rule} processing"
#     # resources:
#     #     time="3:0:0"
#     run:
#         processed_files = os.listdir(os.path.join(MAPPING_DIR,"processed_reads/"))
#         for file in input:
#             samp_id = file.rsplit('/',1)[1]
#             samp_id = samp_id.split(".")[0]
#             processed_reads_path = os.path.join(MAPPING_DIR,"processed_reads", samp_id + ".sorted.RG.bam")
#             if processed_reads_path not in processed_files:
#                 # processed_reads_path = os.path.join(MAPPING_DIR,"processed_reads", samp_id + ".sorted.RG.bam")
#                 shell("module load samtools && samtools addreplacerg --threads 12 -r ID:"+samp_id+" -r SM:"+samp_id+" -o "+processed_reads_path+ " "+ file)

# rule index_RG:
#     input:
#         expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam"), samples =SAMPLES)
#     output:
#         expand(os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam.bai"), samples =SAMPLES),
#         "index_RG.done"
#     message:
#         "Rule {rule} processing"
#     run:
#         for file in input:
#             shell("module load samtools && samtools index -@ 16 "+file)
        # shell("touch index_RG.done")

# rule index_RG:
#     input:
#         rules.add_readgroup_tag.output
#     output:
#         os.path.join(MAPPING_DIR,"processed_reads/{samples}.sorted.RG.bam.bai")
#     message:
#         "Rule {rule} processing"
#     group:
#         "RG"
#     shell:
#         "module load samtools && samtools index -@ 16 {input}"



BEDFILES2, = glob_wildcards("bedfiles/{bed}.bed")
# print(BEDFILES2)
rule run_freebayes:
    input:
        regions = "bedfiles/{bed}.bed",
        reference=ASSEMBLY,
        bams= "bam_list.txt",
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



# rule merge_vcf:
#     input:
#         expand("vcf/{bed}.vcf.gz", bed=BEDFILES2)
#     output:
#         "var_calling/vcf.vcf.gz"
#     message:
#         "Rule {rule} processing"
#     shell:
#         "module load bcftools && bcftools merge -Oz {input} > {output}"

rule sort_index_vcf:
    input:
        rules.run_freebayes.output
    output:
        vcf = "sorted/{bed}.sorted.vcf.gz",
        idx = "sorted/{bed}.sorted.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    shell:
        """"
        module load bcftools samtools
        bcftools sort -Oz -m 2G {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# rule index_vcf:
#     input:
#         rules.sort_vcf.output
#     output:
#          "sorted/{bed}.sorted.vcf.gz.tbi"
#     message:
#         "Rule {rule} processing"
#     shell:
#         "module load bcftools && tabix -p vcf {input}"

rule concat_vcf_individual:
    input: 
        vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
        # chr_list = create_perchr_vcf("chr_list.txt")
    output:
        vcf = "concat/{chromosome}.vcf.gz",
        idx = "concat/{chromosome}.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    # params:
    #     chr = chromosome
    shell:
        """
        module load bcftools;
        bcftools concat --threads 12 --allow-overlaps --remove-duplicates -r {wildcards.chromosome} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """

rule concat_vcf_group:
    input: 
        vcf = expand("sorted/{bed}.sorted.vcf.gz", bed = BEDFILES2),
        regions = "small_chrs.txt"
    output:
        vcf = "concat/small_chrs.vcf.gz",
        idx = "concat/small_chrs.vcf.gz.tbi"
    message:
        "Rule {rule} processing"
    # params:
    #     chr = chromosome
    shell:
        """
        module load bcftools;
        bcftools concat --threads 12 --allow-overlaps --remove-duplicates -R {input.regions} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """
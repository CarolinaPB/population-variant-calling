# Population level variant calling

## First follow the instructions here

[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT

This is a pipeline that takes reads aligned to a genome (in `.bam` format) and performs variant calling with `Freebayes`. It uses VEP to annotate the resulting VCF, calculates statistics, and calculates and plots a PCA.

#### Tools used


- [Freebayes](https://github.com/freebayes/freebayes) - variant calling using short reads
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) - vcf statistics
- [Plink](https://www.cog-genomics.org/plink/) - compute PCA
- R - Plot PCA

| ![DAG](https://github.com/CarolinaPB/pop-var-calling/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |

### Edit config.yaml with the paths to your files

```
ASSEMBLY: /path/to/fasta
MAPPING_DIR: /path/to/bams/dir
PREFIX: <prefix>
OUTDIR: /path/to/outdir
SPECIES: <species>
NUM_CHRS: <number of chromosomes>
```

- ASSEMBLY - path to genome fasta file
- MAPPING_DIR - path to directory with bam files to be used
  - the pipeline will use all bam files in the directory, if you want to use a subset of those, create a file named `bam_list.txt` and add the path to the bam files you want to use. One path per line.
- PREFIX -  prefix for the created files
- OUTDIR - directory where snakemake will run and where the results will be written to  
  If you want the results to be written to this directory (not to a new directory), open config.yaml and comment out `OUTDIR: /path/to/outdir`
- SPECIES - species name to be used for VEP
- NUM_CHRS - number of chromosomes for your species (necessary for plink). ex: 38

## RESULTS

The most important files are and directories are:  

- **<run_date>_files.txt** dated file with an overview of the files used to run the pipeline (for documentation purposes)
- **results** directory that contains
  - **variant_calling/final_VCF** directory with variant calling VCF files, as well as VCF stats
    - {prefix}.vep.vcf.gz - final VCF file
    - {prefix}.vep.vcf.gz.stats
  - **PCA** PCA results and plot
    - {prefix}.eigenvec and {prefix}.eigenval - file with PCA eigenvectors and eigenvalues, respectively
    - {prefix}.pdf - PCA plot

The VCF files has been filtered for `QUAL > 20`. Freebayes is ran with parameters `--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2`. For more details check the Snakefile.

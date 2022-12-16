# Workshop 7: Viral Multi-Omics Workflow

In this workshop, we will explain how to identify viral genomes (vMAGs) from metagenomic assemblies followed by dereplication, annotation, taxonomic classification, and mapping for coverage and relative abundance of recovered viral genomes. We will also cover how to make viral linkages using two main methods (CRISPR and consensus methods).

##Viral Recovery and Identification


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<your email>
#SBATCH --partition=wrighton-hi,wrighton-low

cd ../raw_reads 
# copy reads from ORG-DATA location 
cp /home/ORG-Data-2/Agribiome/CURE_metaG_March2021/15_S10_L003_R* .
# rename reads to sample name
mv 15_S10_L003_R1_001.fastq.gz WCRC_304_R1.fastq.gz
mv 15_S10_L003_R2_001.fastq.gz WCRC_304_R2.fastq.gz
```

# Workshop 7: Viral Multi-Omics Workflow

In this workshop, we will explain how to identify viral genomes (vMAGs) from metagenomic assemblies followed by dereplication, annotation, taxonomic classification, and mapping for coverage and relative abundance of recovered viral genomes. We will also cover how to make viral linkages using two main methods (CRISPR and consensus methods).

## Viral Recovery and Identification

To identify viral genomes within assemblies, we use a VirSorter2. Per their paper, VirSorter2 is “a DNA and RNA virus identification tool that leverages genome-informed database advances across a collection of customized automatic classifiers to improve the accuracy and range of virus sequence detection”. 

Before using this tool, we need to do two things: (1) make sure assemblies are properly renamed and (2) pull only scaffolds of a certain length. 

Important: VirSorter2 can also be used with metaT data to identify RNA viruses, however today we will focus on DNA viruses.

## Step 1: Rename scaffolds in assemblies to include sample names
This step helps in downstream analyses as identified viral genomes will already be named with the sample that they originated from. This has likely been done in previous steps, but if not, use a sed command to do so:

```
sed ‘s/>/>[new-name]/’ [assembly.fa] > [named-assembly.fa]
```

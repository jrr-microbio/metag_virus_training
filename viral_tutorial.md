# Workshop 7: Viral Multi-Omics Workflow

In this workshop, we will explain how to identify viral genomes (vMAGs) from metagenomic assemblies followed by dereplication, annotation, taxonomic classification, and mapping for coverage and relative abundance of recovered viral genomes. We will also cover how to make viral linkages using two main methods (CRISPR and consensus methods).

## Viral Recovery and Identification

To identify viral genomes within assemblies, we use [VirSorter2](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y). Per their paper, VirSorter2 is “a DNA and RNA virus identification tool that leverages genome-informed database advances across a collection of customized automatic classifiers to improve the accuracy and range of virus sequence detection”. 

Before using this tool, we need to do two things: (1) make sure assemblies are properly renamed and (2) pull only scaffolds of a certain length. 

Important: VirSorter2 can also be used with metaT data to identify RNA viruses, however today we will focus on DNA viruses.

### Step 1: Rename scaffolds in assemblies to include sample names
This step helps in downstream analyses as identified viral genomes will already be named with the sample that they originated from. This has likely been done in previous steps, but if not, use a sed command to do so:

```
sed ‘s/>/>new-name/’ assembly.fa > named-assembly.fa
```

Sed commands are very useful and have many options. Here, we are using it to rename the contigs within the assembly. We first call sed and then provide the specific options which are within the ‘ ‘ and delimited by /. The s option (substitution) tells sed that we would like to find [x] and replace it with [z]. In this example scenario, the command would be: sed ‘s/x/z/’. Finally, we tell sed which file to act upon and specify that we want to write this output to a new file using >.  

* *Important: make sure to always replace the > before your new name when re-naming. The carrot before each scaffold name is important for many programs - don’t forget it! 

### Step 2: Pull contigs 10kb in length 
Viral genomes can be challenging to confidently identify from contigs, so to increase the likelihood that genomes we recover are viral we want to run VirSorter2 on relatively long contigs. Therefore, we’ll first pull contigs that are at least 10kbp from our assemblies using pullseq:

```
pullseq.py -i contigs.fa -m 10000 -o contigs_10000.fa
```

In this command, we give pullseq our input assembly (-i contigs.fa) and require a minimum length of 10,000bp (-m 10000) and then write this to a new file (-o contigs_10000.fa). This is the file that you will provide VirSorter in the next step.
Note: For best data management practice, delete the subsetted-assembly file (contigs_10000.fa) once you have run the following steps to identify viral genomes. You can re-pull contigs at any time if needed, and the output file will be the same.

### Step 3: Virsorter2, CheckV and DRAM-v for viral recovery and identification

To confidently identify contigs of viral nature we will use information provided by VirSorter2, CheckV, and DRAM-v through a series of steps. 

We follow the protocol developed by the Sullivan Lab et al. found [here](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=5). Following this SOP, you will first run VirSorter2, followed by CheckV and VirSorter again before finally running DRAM-v. Briefly, VirSorter2 will first identify contigs of possible viral origins and CheckV will then identify and trim the sequences for possible host contamination. We then run VirSorter2 again on the curated sequences to generate an ‘affi-contigs.tab’ file that will feed into DRAM-v and is important for AMG identification. This protocol is designed to check all recovered viral genomes as thoroughly as possible so that you can be sure these are not contaminated with host (bacterial & archaeal) DNA and have enough evidence that the contigs are viral in origin. Therefore, it’s important to also follow the manual curation steps at the end of the SOP which are guided by the results of DRAM-v. 

```
Overview of Sullivan Lab SOP: 
Step 1: Run VirSorter2 → identify viral genomes
Step 2: CheckV → identify and remove possible host contamination 
Step 3 VirSorter2 again → run new, trimmed sequences through virsorter for affi-tab file
Step 4: Run DRAM-v → annotation of viral sequences, needed for manual curation
Step 5: Manual curation → investigate suspect viral genomes 
```

Unlike sed and pullseq commands which can (generally) be run directly on the command line, the above commands for VirSorter2, CheckV and DRAM-v should always be run as bash scripts and submitted to slurm using sbatch [file].sh. Note: To run DRAM-v on the server, the source command is: 

```
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.0
```

#### Stop & Think:
How might size and quality of metagenomes influence the recovery of vMAGs in this step?
Though we recommend using 10kb contigs, in some scenarios you may also want to consider including contigs >5kb 

### Step 4: Cluster vMAGs 95/85
In an effort to standardize viral genome identification, the viral community got together to establish a consensus on best practices in a paper titled [Minimum Information about an Uncultivated Virus Genome (MIUViG)](https://www.nature.com/articles/nbt.4306). Within those rules, and much like you do for MAGs, we cluster the viral genomes to remove any duplicates. Viral clustering is done at 95% ANI across 85% of the shortest contig that is being compared - in other words, to cluster, a viral genome must be 95% similar across 85% of its genome to be considered the same viral population (i.e., vMAG). To do this, we will use some additional features that are included as part of [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) software. 

First, create a blast+ database. This is part of the blastn package already installed on the server so you can just directly call this:

```
makeblastdb -in final-viral-combined-for-dramv_nobadchars.fa -dbtype nucl -out my_db
```

Note: If you are running VirSorter across multiple samples, you want to concatenate all final VirSorter2 outputs and use them as input for this first makeblastdb step.

Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:

```
blastn -query final-viral-combined-for-dramv_nobadchars.fa -db my_db -outfmt '6 std qlen slen' -max_target_seqs 10000 -o my_blast.tsv -num_threads 10
```

Notes: If you have more than 10,000 sequences, increase the -max_target_seqs number to more slightly above that number. “-outfmt 6” is a specific type of format that blastn can output. The headers for the output file will have this [format](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). The following flag “std” is the default order of outfmt 6 output. It then also specifies to give the query length and sequence length in addition to other values. 
 
This next part is done by some custom scripts that are provided as part of the CheckV download (but that aren’t installed on the server). You have to go in and download them from [here](https://bitbucket.org/berkeleylab/checkv/downloads/). Once you unzip that file, go into the “scripts” folder, and the two you need are anicalc.py and aniclust.py. Go ahead and copy those into a directory on the server where you have your unclustered viral genomes. 

After copying them, make sure that permissions are set to read/execute (run: chmod 777 anicalc.py ; chmod 777 aniclust.py) and then calculate pairwise ANI by combining local alignments between sequence pairs:

```
python anicalc.py -i my_blast.tsv -o my_ani.tsv
```

Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):

```
python aniclust.py --fna my_seqs.fna --ani my_ani.tsv --out my_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
```

That spits out a file that has all of your cluster representatives in the first column, and everything that clustered into that cluster on the 2nd column. We extract the IDs of the representatives…

```
awk '{print $1}' my_clusters.tsv > clustered_headers.txt
```

...and to get the final sequence file we use a pullseq script available on the server:

```
pullseq_header_name.py -i final-viral-combined-for-dramv_nobadchars.fa -n clustered_headers.txt -o final_95-85_clustered_vMAGs.fasta -e F
```

The result is a final file called final_95-85_clustered_vMAGs.fasta that is your final clustered viral database.

### Step 5: Rerun DRAM-v and delete old unclustered DRAM-v
Now that you have your final clustered database, we want to run DRAM-v once again so that we can get all of our viral annotations as well as auxiliary metabolic genes (AMGs). We run DRAM-v similar to how the Sullivan lab SOP ran it, but we add the Uniref flag to get as much information as possible. 

```
#!/bin/bash
#SBATCH --nodes=1                  # Number of requested nodes
#SBATCH --ntasks=10               # Cores per node requested
#SBATCH --time=7-50:00:00        # Timelimit
#SBATCH --mem=170gb                # Job memory (1024gb=1tb)
#SBATCH --mail-type=BEGIN,END,FAIL         
#SBATCH --mail-user=jrodram@colostate.edu
#SBATCH --partition=wrighton-hi
#SBATCH --nodelist=zenith

source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.0

#annotate
DRAM-v.py annotate -i final_95-85_clustered_vMAGs.fasta -v viral-affi-contigs-for-dramv_nobadchars.tab -o ./output_clustered/uniref_clustered_vMAGs_DRAMv --min_contig_size 0 --threads 15  --use_uniref &> log_dramv_clustered_viruses_uniref.txt

#distill
DRAM-v.py distill -i annotations.tsv -o ./distill_clustered_vMAGs_output/
```

#### Stop and think:
Once this DRAM-v run finishes, go back and delete the older DRAM-v run, as you will not need this anymore. Additionally,  the &> log_dramv_clustered_viruses_uniref.txt file in the last command is a log output that is from the linux ecosystem. DRAM-v now outputs a log file so you do not necessarily need that - however I like to include just in case.

With the final distill command, you get an output table called “amg_summary.tsv” that has all of the AMG information you will need for AMG analyses. This will give you AMGs that have their own categories (as labeled in the DRAM manuscript). We usually take AMGs that are categories 1-3, and remove those that are categories 4-5. However, another thing to watch out for is transposons. For these putative AMGs, make sure that they are not labeled with the “T” flag and “F” flags, as these are less likely to be viral - and even more so the “genome” itself is possibly not viral. If you have multiple AMGs per viral genome, it is critical that you go to the annotations and manually confirm again that it is, in fact, a viral genome. Usually when there are more than 3 AMGs within a single viral genome, I consider this slightly suspect. You can manually confirm these are viral by looking for viral-like genes within the annotations like the Sullivan Lab SOP mentions. I mention this because sometimes high confidence ‘keep’ labeled viral genomes end up in the post-qc steps and have no indication of being viral.

### Step 6: Viral Taxonomy
Similar to bacteria and archaea, we often want to know whether the viral genomes we have identified are novel. However, since viruses do not have universal marker genes, we cannot use normal methods like comparing a specific gene (unless it is a broadly conserved gene for a specific viral group). In virus-land, we identify taxonomy by assessing protein content similarity between different viruses. To do this, we use a tool called vContact2. The input is the .faa proteins file from the DRAM-v output.

Before we get started, we need to subset/parse out our sequences a little. Rory and Josué worked on this script that can help manually parse vContact2 output based on scaffold headers. First we change the scaffold headers to a line that the script can parse downstream:

```
sed ‘s/>/>this_study_/’ DRAMv_genes.faa > renamed_DRAMv_genes.faa
```

After this is done, we can source vContact2 and run the commands:

```
source /opt/Miniconda2/miniconda2/bin/activate vContact2
```

This next command preps the files for vContact2 input. It basically indexes all your gene IDs per genome and writes them into proper format for vContact2 use:

```
vcontact2_gene2genome -p renamed_DRAMv_genes.faa -o viral_genomes_g2g.csv -s 'Prodigal-FAA'
```

Once you have the gene2genome csv file, you can run vContact2. Don’t worry too much about these flags - these are default and recommended. The only thing that might change through time is “–db”. Right now, version 211 is the latest version of Viral Refseq. However, this gets updated somewhat often so keep an eye out.

```
#Run actual vContact2 command
vcontact2 -r renamed_DRAMv_genes.faa --rel-mode Diamond -p viral_genomes_g2g.csv --db ProkaryoticViralRefSeq211-Merged --pcs-mode MCL --vcs-mode ClusterONE --pc-evalue 0.0001 --reported-alignments 25 --max-overlap 0.8 --penalty 2.0 --haircut 0.1 --pc-inflation 2 --vc-inflation 2 --min-density 0.3 --min-size 2 --vc-overlap 0.9 --vc-penalty 2 --vc-haircut 0.55 --merge-method single --similarity match --seed-method nodes --sig 1 --max-sig 300 --mod-inflation 5 --mod-sig 1 --mod-shared-min 3 --link-sig 1 --link-prop 0.5 --c1-bin /opt/Miniconda2/miniconda2/envs/vContact2/bin/cluster_one-1.0.jar -o vContact2_output -t 15

#deactivate the source (this will also happen automatically if you just close your terminal)
source /opt/Miniconda2/miniconda2/bin/deactivate 
```

The vContact2 file contains every single viral genome you gave as input, as well as every single viral genome that was added into the network as part of the tool. The column “VC” corresponds to “Viral Cluster”. If you sort by this column, everything that shares the same “VC” means they have enough protein content similarity to be considered a viral “genera”. If the virus shares a cluster ID with only viruses from your site of study, that cluster is considered a “novel viral genera” since they did not cluster with any of viral sequences in the known database. If the virus clusters with a virus from the RefSeq / Viral databases, then it is a “known taxonomy” virus that shares taxonomy with whatever genome it clustered to in the RefSeq database.  If there are multiple possibilities within a single known taxonomy cluster, the recommendation is usually to go with whatever the “majority” dictates.

This is quite tedious to do manually, and so now that you have vContact2 run and the output file we can work on parsing it with the script. Open the output file “genome_by_genome_overview.csv”. Add in a second column right after the “Genome” column named “genome_category”. Then, filter everything in the “Genome” column that has “this_study” in its header name (these correspond to the sequences that we input into the database) and populate the cells in the new column as “this_study”. Now, filter everything that is NOT in your study, and populate the cells in the new column as “refseq_genome”. Once this is done, you can run the actual parser:

```
vContact2_parser_v2.py -i genome_by_genome_overview.csv -o parsed_genome_by_genome_overview.csv
```

The resulting parsed genome by genome overview will then contain the files that you need to establish taxonomy for your viruses. It populates a column called “auto type” that includes the different viral groups that are present in each viral cluster. This is much easier to get overall quantities than manually going into each cluster and verifying what viruses are in it.

Note: If you are trying to use this script to say more than just whether a virus has known taxonomy, the genome_category column can be leveraged to collate each cluster into “groups”. The way that you can use this is by adding more than just the two categories “this study” and “refseq_genome” into the possible options in the “genome_category” column. For example, if you have viruses from another ecosystem, and you are seeing if any of yours cluster to another system, add in a third category “ecosystem_2” to “genome_category”. The parsed output file will then tell you when a genome from “this study” will cluster to “ecosystem_2” or a “refseq_genome”.

### Step 7: mapping reads to viruses
Similar to MAGs, we can also map metagenomic reads to viral genomes to determine coverage and relative abundance. This follows a very similar workflow as mapping to MAGs or genes, however there are some key differences. Here we include an example for mapping with bbmap, but the same workflow can be followed using bowtie (as was shown with mapping to MAGs). 

First, map trimmed reads to a concatenated file of all viral contigs using bbmap: 

```
bbmap.sh in1=R1_All_trimmed.fastq in2=R2_All_trimmed.fastq ref=FINAL_all78metaG_4893viral_curated_95-90.fna threads=30 out=mapped.sam 
```

The above command should always be run as a bash script and submitted to slurm. The above command tells bbmap that we would like to map trimmed reads (in1= and in2=) to our reference file of concatenated viruses (ref=). Of course, mapped.sam is our output (out=) and I am requesting this job to be run on 30 processors (threads=).

As with any mapping job, we next want to convert our sam file to a bam file. Reminder: sam files are large and can be deleted after conversion to a bam file. Bam your sams! The ‘@’ symbol in the below command designates processors.Then, we filter for minID of 95% using the reformat command and finally sort our sam file again using samtools.

```
samtools view -@ 30 -bS mapped.sam > mapped.bam
reformat.sh in=mapped.bam out=mapped_95id.bam minidfilter=0.95 primaryonly=t pairedonly=t
samtools sort -@ 15 -o mapped_95id_sorted.bam mapped_95id.bam
```

Next, we’ll use two coverM commands to determine coverage which can be used to calculate relative abundance. Given vMAGs are viral contigs, coverM contig mode is applied with two commands. First, --min-covered-fraction 75 and next followed by -m reads_per_base to calculate coverage. Similar to requirements set for MAGs, here vMAGs must have a minimum covered fraction >75% to be considered present. Any sequences that do not meet this threshold should not be considered in the second output using reads_per_based. Coverage values can then be calculated from the reads per base output x 151 bp. It is also important to consider your depth cut off here, generally 1-3x is acceptable to be considered ‘present’.

Since you often map many reads instead of just one pair, it’s useful to use the linux wildcard “ * “ to tell coverM to consider all files that end with the same suffix, which in this case is “__95id_sorted.bam”. 

```
coverm contig --bam-files *_95id_sorted.bam -m covered_fraction -t 15 --min-covered-fraction 75 &> coverM_covered_fraction_75_clusteredviruses.txt

coverm contig --bam-files *_95id_sorted.bam -m reads_per_base -t 15 &> coverm_reads-per-base_clustered.txt
```

### Step 8: Make a master spreadsheet
It is important to start generating a viral supplementary file early, as all of these files will be required upon publication. This is also very helpful for analyses to have much of the important information in one place. At this point, I usually start making a spreadsheet in excel that includes the following information:

Tab 1: vMAG Summary: Virus genome ID, any viral ID renaming (sometimes this is helpful), genome average abundance, checkV scores and quality metrics from checkV output files, genome length, and viral taxonomy
Tab 2: vMAG annotations
Tab 3: AMG Summary/Distill output
Tab 4: vMAG read mapping output


### Step 9: Making host-viral linkages
Once you have a database of both recovered MAGs (host genomes) and recovered vMAGs (viral genomes), you can then attempt to make host-viral linkages. Since viruses depend on their hosts to proliferate, we can better understand the role that a viral population may have in the broader all microbial community if we understand which hosts they interact with. 

There are two main approaches that we use to make connections between a MAG and a vMAG. These are (1) CRISPR based linkages, (2) consensus method using tools that consider oligonucleotide frequency and sequence similarity between hosts and viruses. Generally CRISPR based linkages are considered to be the strongest way to make host-viral linkages. The consensus method does not depend on hosts containing a CRISPR array and instead uses two (maybe even 3) different tools to determine likely virus-host linkages. For the consensus approach, we recommend  VirHostMatcher and PHIST. Additional tools that can be used include software like WiSH and VirHostMatcher-Net.

### Option 1: CRISPR-based host-virus inkages
This method of making host-virus linkages depends upon your host MAG containing a CRISPR array. Unfortunately, CRISPR arrays are more commonly encoded and deployed in some ecosystems compared to others and therefore it may not be surprising if you don’t detect many (or any) CRISPR arrays in and across your MAGs. If this is the case, you’ll want to pair or substitute this approach with the consensus approaches (options 2 and 3).

Why is this approach the strongest way to make host-virus linkages? CRISPR-Cas systems are often thought of as bacterial and archaeal immune systems. CRISPR works by recording memories of viral interactions by integrating small pieces of viral DNA as spacers within the hosts’ CRISPR array. Within the array, spacers are interspaced with identical repeat sequences and flanked by Cas (CRISPR-associated) genes. These saved memories help to protect the host against recurrent invasion by the same viral population by more rapidly identifying and degrading the invading nucleic acids from the virus. For our use, we can extract these spacers (which are incorporated snippets of viral DNA) and match them to our vMAGs.

The basic overview to making host-viral linkages vira CRISPR is:
Identify CRISPR arrays within a MAG
Extract spacers 
BLASTn these spacers against our vMAG database 
Filter for ‘good’ matches to make high-confidence host-viral linkages 

The first thing you’ll need to do is identify CRISPR arrays in your MAGs. To do this, we will use Geneious. You’ll also need to install the CRISPR Recognition Tool (CRT) plugin for Geneious. 

Once you have this installed and Geneious opened, import your MAGs as individual fastas. Select one a MAG, and go to ‘Annotate & Predict’ on the top toolbar. Select ‘Find CRISPR loci’ from the dropdown menu. You will then see a display window with some requirements in identifying CRISPR arrays - make sure to select ‘more options’ and enter the following parameters: 
min number of repeats a CRISPR must contain: 4 
minimum length of a CRISPR’s repeated region: 19 
maximum length of a CRISPR’s repeated region: 55
minimum length of a CRISPR’s non-repeated region (or spacer region): 19 maximum length of a CRISPR’s non-repeated region (or spacer region): 48 
length of a search window used to discover CRISPR’s: 8

You will see that CRISPR arrays are identified in green and yellow highlighted sequences within contigs of the MAG fasta file. Explore this! 

Stop and think: 
CRISPR arrays are made of many short repeated sequences (quite literally the ‘repeat’ sequences that interspace the spacers) and therefore often break metagenomic assemblies. Because of this, you’ll see that many arrays are often found at the end of contigs and you should not assume that every highlighted ‘array’ is a distinct array. If you are interested in counting the number of distinct arrays within a MAG, you will want to count the number of different repeat sequences, as these should be nearly identical within a given array.

Next, we can extract these spacers by going to Tools → Extract Annotations → make sure the prompted sentence reads “Annotation type / Is / CRISPR Spacer”. This will create a new file of only spacers, which you can export via File → Export → Document and select as a fasta.

Note: you can do much of this in bulk - you do not need to do each MAG individually. Make sure that scaffolds within bins are properly renamed before beginning.

Now that we have a fasta file of our spacers, we need to BLASTn them to our vMAG database to make linkages. First, we’ll first make a nucleotide blastDB (-dbtype nucl) with a file of our concatenated viruses: 

```
makeblastdb -in concatenated-viruses.fna -out concatenated-viruses.fna_DB -dbtype nucl
```

Next, BLAST the spacers (-query) using the database of viruses we just made: 

```
blastn -query spacers.fa -db concatenated-viruses.fna_DB -outfmt 6 -out BLAST_spacers_output.txt -evalue 1e-5 -num_threads 5 -word_
size 7 -max_target_seqs 10000 -dust no
```

This BLAST command is specialized for short sequences by using the following flags: -word_size, -max_target_seqs and -dust. When blasting spacers at any time, make sure to include these flags. 

By specifying out format 6 (-outfmt), the headers in the blast output text file will be:
query sequence id → these are your spacers that hit 
subject sequence id → these are the virus the spacer hit to
percentage of identical matches
alignment length
number of mismatches	
number of gap openings	
start of alignment in query	
end of alignment in query
start of alignment in subject
end of alignment in subject
Bitscore

Each row in the output file lists each spacer (column 1) that has hit to a virus (column 2) with the next 9 columns providing details on how ‘good’ each hit is. We want to filter this to only high quality hits. Current standards allow 1 mismatch between spacers + virus with 0 gap openings. 

Once you have filtered the BLAST output, you now have a list of host spacers that link to a virus, and thus host-virus linkages!

### Option 2: VirHostMatcher
[VirHostMatcher](https://academic.oup.com/nar/article/45/1/39/2605663) uses oligonucleotide frequency between a viral and host genome to determine if they infect each other. 

To run VirHostMatcher, you create three different directories. One of them will contain all the viral genomes (i.e., “virus), one will contain all the bacterial / archaeal genomes (i.e., “host”), and the last will be an empty folder titled “output”. Note: For viral and host genomes, they need to be in a single .fasta file for each. So if you have 125 genomes, you will have 125 fasta files one with each genome.

Once you have your folders set up and your genomes in place, you can run the command:

```
python /opt/VirHostMatcher/vhm.py -v virus -b host -o output
```

In the output folder, you will then have a BUNCH of different files for all of the distance metrics that this method uses. You will only be using the file titled “d2star_k6.csv”. This file is a matrix where hosts are the columns, and viruses are the rows. The different scores that are within each matrix are the dissimilarities. The important thing about VirHostMatcher is that you ONLY consider a possible hit if it is < 0.25. ie., you need to parse out this file to only have the linkages that are under 0.25. After that, for each virus, you only choose the LOWEST d2* metric as the most likely hit. If a virus has multiple hits under 0.25, you only take the lowest hit. While theoretically possible that a virus infects multiple hosts, this tool is quite squishy so taking a stringent approach is best.

### Option 3: PHIST
[PHIST](https://academic.oup.com/bioinformatics/article/38/5/1447/6460800) is different from VirHostMatcher because instead of using oligonucleotide frequencies, it uses exact sequence similarity between a virus and its putative hosts. 

```
phist.py -k 25 -t 25 ./virus ./host common_kmers.csv predictions.csv
```

The output predictions.csv is a little more straightforward here! For this one, you can take a virus (column phage) and a host (column host) as being linked if the adj-pvalue is less than 0.05.

Once you have both tools run, you then subset / select the virus-host links that are confidently predicted by both VirHostMatcher ( d2* < 0.25) and PHIST (adj.-pvalue < 0.05) to the same hosts. Note: You can have more than 1 putative host with PHIST, but NOT with VirHostMatcher. You ONLY take the lowest d2* measurement.

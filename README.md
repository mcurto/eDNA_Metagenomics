# eDNA Metagenomics: Scripts and codes used in Curto, M., Veríssimo, A., Santos, C. D., Riccioni, G., Ribeiro, F., Jentoft, S., Alves, M. J., Gante, H. F. (in prep.) Improving whole biodiversity monitoring with eDNA metagenomics. 


In this manuscript, Illumina shot gun sequencing data was produced for environmental DNA (eDNA) water samples from the Ave River to evaluate eDNA metagenomics applicability in retrieving biological diversity information. In this scope, we evaluated first the ability of this method in describing the whole biological community and then the fish community. For the latter a comparison was made with eDNA metabarcoding using the 12S primers from Miya et al (2015). Moreover, we performed several tests comparing inter- and intraspecific identity distribution of several publicly available shotguns sequencing datasets with the reference genomes of Cyprinus carpio, Anguilla anguilla and Anguilla rostrata. Therefore, this repository is divided in three parts:

- [eDNA metagenomics shot gun sequencing data analysis](#edna-metagenomics-shotgun-sequencing-data-analysis)
- [Metabarcoding analysis](#metabarcoding-analysis)
- [Sliding window analysis to obtain identity distributions](#sliding-window-analysis-to-obtain-identity-distributions)



# eDNA metagenomics shotgun sequencing data analysis

The bioinformatic pipeline for read classification consisted in three steps:

- [Quality control and paired read merging](#quality-control-and-paired-read-merging)
- [Homology search with blast](#homology-search-with-blast)
- [Taxonomic classification](#taxonomic-classification)


## Quality control and paired read merging

Quality control was done with Trimmomatic (Bolger et al. 2014) and reads were merged with PEAR (Zhang et al. 2014). Trimmomatic trimmed the SmarterThruPlex adapters and low-quality nucleotides at the 3’ while pear run with default options. In both programs reads smaller than 50bp were excluded.

Both programs can be run for multiple samples using the shell script **Run_Trimomatic_and_PEAR.bash**. Three positional arguments are used for this script. 1st the directory containing the paired fastq files compressed as gunzip, 2nd the directory to save the processed fastq files, and 3rd the path containing Trimmomatic executable and directory with adapter sequences:

    bash Run_Trimomatic_and_PEAR.bash /path/to/input_fastq.gz/ /path/to/outpu_fastq.gz/ /path/to/Trimmomatic/

All outputs are saved in the gunzip compressed format. The resulting files were then converted into fasta farmat.

## Homology search with blast
The program blastn from ncbi-blast-2.12.0+ was used to find homology between the shot-gun reads and different reference databases: 1) The nucleotide database from genebank (nt); 2) Fish genomes of all species with resources available for all families present in Portugal plus the genomes from fish that are commonly eaten by humans (genome); 3)Fish transcriptomes of all species with resources available for all families present in Portugal plus the transcriptomes from fish that are commonly eaten by humans (transcriptomes). The program ran with the following code:

    blastn -query /path/to/Metagenomics_reads.fasta \
    -db /path/to/nt/database/nt -perc_identity 90 -qcov_hsp_perc 0.9 \
    -out /path/to/Blast_output.blastout -evalue 1e-10 \
    -outfmt "6 qseqid sseqid staxids pident qlen length evalue bitscore" \
    -num_threads 10

The R1 and R2 files from unmerged reads were blasted separately using the same parameters.

## Taxonomic classification

The nucleotide database was already prebuilt with taxonomic information which was not the case for the remaining ones. For that reason, their blast outputs had to be processed differently. In both cases, the output files were analyzed in four steps: i) combined unassembled paired reads; ii) filter the results based on e-value, match identity, and query coverage; iii) adding the taxonomic lineage to the matches; iv) summarize the final taxonomy per read. 

### nt database:

#### Combined unassembled paired reads:

The blast results from unassembled paired reads were merged with the script **MergeR1AndR2.py**. This script takes four positional arguments, in this order: 1<sup>st</sup> the path to the read 1 blast outputs, 2<sup>nd</sup> the path to the read 2 blast outputs, 3<sup>rd</sup> the path to the fasta files, 4<sup>th</sup>  the output directory.

    python3 MergeR1AndR2_v2.py /path/to/directory/blastout_Read1/ /path/to/directory/blastout_Read2/ /path/to/directory/fasta/ path/to/directory/output/

The resulting files were concatenated with the assembled output.

#### Filter results:

Blast results were filtered using the script **FilterBlast.py**. This takes five positional arguments: 1<sup>st</sup> the path to the a directory containing the blast outputs that should have the ".blastout" extension; 2<sup>nd</sup> the path to the directory where the filtered files should be saved; 3<sup>rd</sup> the minimum e-value to keep a match; 4<sup>th</sup> the minimum identity; 5th the minimum query coverage. Here an example with the parameters used in the paper:

    python3 FilterBlast.py /path/to/directory/Blast_outputs/ /path/to/directory/Blast_filtered_outputs/ 1e-10 99 0.9

This script will merge matches showing the same parameters.

#### Add taxonomic lineage to results:

Taxonomic lineage was added in three steps. First the taxonomic IDs were extracted with the script **ExtractTaxID_Filt.py**. Second the lineage was retrieved using taxonkit (Shen and Ren 2021). Third lineage information was added to the filtered blast outputs with the script **AddTaxonomyFiltFile_WithTaxIDInfo.py**

The script ExtractTaxID_Filt.py takes two positional arguments: 1<sup>st</sup> the directory containing the filtered blast outputs, and 2<sup>nd</sup> the path to the output file where the taxonomic IDs will be saved:

    python3 ExtractTaxID_Filt.py /path/to/directory/Blast_filtered_outputs/ /path/to/TaxID_file.txt

The resulting Tax ID file was used to extract the lineages using the taxonkit tool. This step requires the taxonomy database from genbank (https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) to be downloaded and decompressed.

Taxonkit can be run in the following way:

    taxonkit lineage --data-dir /path/to/taxdump/ /path/to/TaxID_file.txt -o /path/to/Lineage_file.txt

The lineages are then added to the filtered files with the script **TaxonomyFiltFile_WithTaxIDInfo.py**. This will be done by adding an extract column to the blast output with that information. This script takes three positional arguments: 1<sup>st</sup> the directory containing the filtered reads, 2<sup>nd</sup> the path to the file containing the lineages, 3<sup>rd</sup> the output directory:

    python3 AddTaxonomyFiltFile_WithTaxIDInfo.py /path/to/directory/Blast_filtered_outputs/ /path/to/Lineage_file.txt /path/to/directory/Blast_lineage_outputs/


#### Summarize taxonomy:

Taxonomy can be summarized by defining a final lineage per read with the script **SummariseTaxonomyPerRead.py**.

In case multiple best matches are found to a certain read, the script **SummariseTaxonomyPerRead.py** will find the most recent common ancestor among them. The resulting output will be a tab separated text file containing the following information per line: read name, match identity, read length, alignment length, e-vaule, bit score, and final lineage. The script needs two positional arguments: 1<sup>st</sup> the directory containing the lineage files, 2<sup>nd</sup> the directory to save the outputs per sample:

    python3 SummariseTaxonomyPerRead.py /path/to/directory/Blast_lineage_outputs/ /path/to/directory/Final_lineage_per_read/


### Genome and transcriptome databases:
The analysis for the genomes and transcriptomes databases only differed from the nt analysis in the step where lineage information is added. More precisely, this step is done before the filtering step with the script **AddLineage2Genomes.py**. For that two tab-separated tables need to be prepared:

1) A correspondence between accession number, taxon and TaxID as follows:

|   |   |   |
|---|---|---|
| accession number 1 | taxon 1 | TaxID 1 |
| accession number 2 | taxon 2 | TaxID 2 |

2) A correspondence between TaxID and lineage:

|   |   |
|---|---|
| TaxID 1 |  lineage1 |
| TaxID 2 |  lineage2 | 


Then the script **AddLineage2Genomes.py** is used to add both TaxID and lineage information to the filtered blast files. It takes four positional arguments: 1<sup>st</sup> the list containing the  correspondence between accession number, taxon and TaxID; 2<sup>nd</sup> the correspondence between TaxID and lineage; 3<sup>rd</sup> the directory where the filtered were saved; and 4<sup>th</sup> the directory to save the output files: 

    python3 AddLineage2Genomes.py /path/to/Accession_species_TaxID_correspondance.txt /path/to/TaxID_Lineages_correspondence.txt path/to/directory/Blast_outputs/ /path/to/directory/Blast_lineage_outputs/


### Merging and filtering results from different databases
The results from different databases were merged for fish taxa. This was done by:

**1)** merging and sorting the results from different databases with the following command:

    cat path/to/result/nt/database/* path/to/result/genome/database/* path/to/result/transcriptome/database/* | sort -k1 > path/to/all_blastout_results_concatenated.blastout

**2)** Filter for the best match and merge equally good matches using the similar criteria to the filtering step. The main difference is that the e-value is not used since it is database dependent. The script needs two positional arguments: 1<sup>st</sup> the directory containing the input files, 2<sup>nd</sup> the directory to save the output files:

    python3 FilterBlast_MergedPerRead.py path/to/concatenated/files/ path/to/concatenated/files_filtered/

**3)** Summarize the results per read as described before with the script **SummariseTaxonomyPerRead_BestDB.py**:

    python3 SummariseTaxonomyPerRead_BestDB.py path/to/concatenated/files_filtered/ path/to/concatenated/files_filtered_sumtax/

These results were then filtered to ensure that only fish matches were kept. This was done by checking if there was another sub-optimal match to a fish taxon. To do that, the fish matches were extracted from unfiltered files for all comparison using grep and the taxonomic information for all these matches was retrieved as described above (see Add taxonomic lineage to results section). The script **Positive_read_check.py** was used to filter the results. The script takes four arguments: 1<sup>st</sup> the merged output with the best matches across all databases, 2<sup>nd</sup> the blast output file with all matches for the same reads across all databases with lineage information included, 3<sup>rd</sup> the target taxonomic group, and 4<sup>th</sup> a prefix to save output files. In this case two output files are saved. One with the filtered data (.filt.blastout) and a second with the number of matches the specified taxon (.check_data.blastout). Here an example of a code:

    python3 Positive_read_check.py path/to/concatenated/files_filtered_sumtax.blastout path/to/all_matches_lineage.txt Actinopteri  path/to/output/directory/prefix    
 
    
# Metabarcoding analysis

Amplicon sequence variants (ASVs) for metabarcoding analysis were obtained with the R package DADA2 as described in the script **dada2_12S_Metabarcoding.R** and saved on a ASV table.

The sequences were extracted in the fasta format using the script **extract_ASVs.py** as follows:

    python3 extract_ASVs.py ASV_table.txt ASV_sequences.fasta

Reads taxonomic assignment was obtained as described for the shotgun sequencing analysis.


# Sliding window analysis to obtain identity distributions

Identity distributions resulting from the comparison of publicly available shotgun sequencing data were estimated within sliding windows running over the alignment between the shot gun reads and the reference genome. This was done using the script **Sliding_window_ID.py** that requires a vcf file containing the allele call from this comparison and the correspondent average read depth per position. Here we present the complete analytical pipeline showing how this was done using the comparison between the reads from Carassius auratus (DRR172221) with the Cyprinus carpio genome (GCF_018340385.1).


The downloaded shotgun sequencing files were quality controlled with trimmomatic. Here an example for the sample DRR172221:

    java -jar trimmomatic-0.39.jar PE -threads 6 -phred33 DRR172221_1.fastq.gz DRR172221_2.fastq.gz DRR172221_1.paired.fastq.gz SRR13304660_1.unpaired.fastq.gz DRR172221_2.paired.fastq.gz SRR13304660_2.unpaired.fastq.gz ILLUMINACLIP:All_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70

Then they were mapped to the reference genome with bwa (Li and Durbin 2009) and converted to bam files with samtools (Li et al 2009). Here an example with the mapping to the genome with reference GCF_018340385.1:

- Referene genome indexing

      bwa index GCF_018340385.1_ASM1834038v1_genomic.fna.gz

- Mapping to reference

      bwa mem -M -t 10 GCF_018340385.1_ASM1834038v1_genomic.fna.gz DRR172221_1.paired.fastq.g DRR172221_2.paired.fastq.g | samtools view -bS - > SRR13304660.bam

PCR duplicates were marked with Picard after sorting the bam files:

- Bam sorting 
    
      java -jar picard.jar SortSam I=DRR172221.bam O=DRR172221-sorted.bam SORT_ORDER=coordinate
 
- Marking duplicates

      java -jar picard.jar MarkDuplicates I=SRR13304660-sorted.bam O=DRR172221-sorted-md.bam M=DRR172221-md-metrics.txt


And the resulting bam files were indexed with samtools and used for variant calling with freebayes (Garrison and Marth 2012):

- Indexing:
    
      samtools index DRR172221-sorted-md.bam

- Variant calling:

      freebayes -f GCF_018340385.1_ASM1834038v1_genomic.fna DRR172221-sorted-md.bam > SRR13304660-sorted-md.freebayes.vcf

Indels were removed with vcf tools (Danecek et al. 2011):
    
    vcftools --vcf DRR172221-sorted-md.freebayes.vcf --remove-indels --recode --recode-INFO-all --out DRR172221-sorted-md.freebayes.snps-only.vcf

Depth per position is calculated with samtools:
    
    samtools depth -a DRR172221-sorted-md.bam > DRR172221-sorted-md.tsv

Then the scrip **Sliding_window_ID.py** is used to estimate the degree of identity per sliding window. The script requires five positional arguments:  1<sup>st</sup> sliding window size, 2<sup>nd</sup> vcf file, 3<sup>rd</sup> depth file, 4<sup>th</sup> minimum depth, 5th output file:

    python Sliding_window_ID.py 189 SRR13304660-sorted-md.freebayes.vcf DRR172221-sorted-md.tsv 3 DRR172221.nr_SNPs.txt


# References:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Zhang, J., Kobert, K., Flouri, T., & Stamatakis, A. (2014). PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), 614-620.

Shen, W., & Ren, H. (2021). TaxonKit: a practical and efficient NCBI taxonomy toolkit. Journal of Genetics and Genomics, 48(9), 844-850.

Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

Li H., Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.

Picard Toolkit. (2019). Picard toolkit. Broad Institute, Github Repository.

Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907.

Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156-2158.

Miya, Masaki, et al. "MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species." Royal Society open science 2.7 (2015): 150088.


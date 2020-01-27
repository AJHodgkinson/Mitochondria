# Mitochondrial Gene Expression and RNA Modification

In this area you will find instructions and code to be able to generate mitochondrial-encoded 
gene expression quantifications and RNA modification levels, as described in the following publications:

Ali, A.T., *et al*. 2019. Nuclear Genetic Regulation of the Human Mitochondrial Transcriptome. *eLife* **8**: e41927
Ali, A.T., *et al*. 2019. Nuclear genetic regulation of human m1A/G mitochondrial RNA modification. *BioRxiv*: 666339

## Dependencies

Perl Module Requirements: 

Parallel::Loops;
Getopt::Long;
Pod::Usage;

Software Requirements:

[trim_galore] (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
[STAR] (https://github.com/alexdobin/STAR)
[Samtools (1.0 or newer)] (http://www.htslib.org);
[Pigz] (http://zlib.net/pigz/);

Aliases:

SAMPLE_NAME = Name of sample used

## Aligning and Filtering Data

Trim Reads for adaptors, quality and polyA tails:

	mkdir SAMPLE_NAME.trim
  trim_galore -stringency 3 -q 20 -o SAMPLE_NAME.trimQ20 --phred33 --paired SAMPLE_NAME_1.fq.gz SAMPLE_NAME_2.fq.gz



Map RNA Sequencing Data to a reference genome:

	STAR --genomeDir STARgenome2pass --runThreadN threads --alignEndsType EndToEnd --readFilesIn fastq_1.fq.gz
	fastq_2.fq.gz --readFilesCommand zcat --outFileNamePrefix Outname --outFilterMismatchNmax 10


# Mitochondrial Gene Expression and RNA Modification

In this area you will find instructions and code to be able to generate mitochondrial-encoded 
gene expression quantifications and RNA modification levels, as described in the following publications:

Ali, A.T., *et al*. 2019. Nuclear Genetic Regulation of the Human Mitochondrial Transcriptome. *eLife* **8**: e41927

Ali, A.T., *et al*. 2019. Nuclear genetic regulation of human m1A/G mitochondrial RNA modification. *BioRxiv*: 666339

## Dependencies

**Perl Module Requirements:** 

Parallel::Loops;

Getopt::Long;

Pod::Usage;

**Software Requirements:**

[trim_galore] (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/);

[Prinseq] (http://prinseq.sourceforge.net);

[STAR] (https://github.com/alexdobin/STAR);

[Samtools (1.0 or newer)] (http://www.htslib.org);

[Pigz] (http://zlib.net/pigz/);

**Aliases:**

SAMPLE_NAME = Name of sample used;

REFERENCE_GENOME = Reference genome to align to;

READ_LENGTH = Length of reads;

MISMATCH = Number of mismatches to allow per read pair;

BAM_FILE = Final BAM file for processing;

## Aligning and Filtering Data

Trim Reads for adaptors, quality and polyA tails:

	mkdir SAMPLE_NAME.trim
 
	trim_galore -stringency 3 -q 20 -o SAMPLE_NAME.trimQ20 --phred33 --paired SAMPLE_NAME_1.fq.gz SAMPLE_NAME_2.fq.gz

	perl prinseq-lite.pl -fastq SAMPLE_NAME.trim/SAMPLE_NAME_1_val_1.fq -fastq2 SAMPLE_NAME.trim/SAMPLE_NAME_2_val_2.fq -out_good SAMPLE_NAME.trim/test -trim_tail_left 5 -trim_tail_right 5 -min_len 20

Map RNA Sequencing Data to a reference genome:

	STAR --genomeDir REFERENCE_GENOME --runThreadN 10 --readFilesIn SAMPLE_NAME.trim/SAMPLE_NAME_1_val_1.fq.gz SAMPLE_NAME.trim/SAMPLE_NAME_2_val_2.fq --readFilesCommand zcat --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang (READ_LENGTH - 1) --outFilterMismatchNmax MISMATCH --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:SRR603365 PU:Illumina PL:Illumina LB:SAMPLE_NAME SM:SAMPLE_NAME CN:Seq_centre --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000 --twopassMode Basic --outFileNamePrefix SAMPLE_NAME.STAR.pass2.

	samtools index SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.bam

Keep only properly paired reads:

	samtools view -@ 10 -f 0x0002 -b -o SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.PP.bam SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.bam

	samtools index SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.PP.bam

Keep uniquely mapped reads

	samtools view -h SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.PP.bam | grep -P "NH:i:1\t|^@" | samtools view -bS - > SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

	samtools index SAMPLE_NAME.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam
	
## Collect count data for modification inference:

	perl pileupAlleleExtractor_mito.pl --Bam BAM_FILE --Out SAMPLE_NAME

### USAGE

perl pileupAlleleExtractor_mito.pl [options] --Fastq1 <forward_fastq_file> --Out <out_name>

**Options:**

--Bam <bam_file> : File containing aligned data (Required)

--Out <out_name> : Out name for file stub (Required)

--MinQ <quality_threshold> : Phred score cut-off to include nucleatide in count (Optional)

--RefFasta <path_to_ref_fasta> : Specify path to reference fasta sequence (all in one file) (Required if using STAR to map, or building a reference genome with Tophat2).

--SamtoolsPath <path_to_samtools> : Use to specify path to samtools, default assumes it is in current path (Optional).

--noBaq:  Use this flag with no arguments to turn of BAQ tool in Samtools

**Perl Module Requirements:**

Getopt::Long
Pod::Usage

### OUTFILE

This tool will create an outfile that can be used to count alternative alleles at each site.

Column Headers: Chromosome, Reference Position, Reference allele, Nucleotides(Forward strand upper case, reverese strand lower case):Nucleotide Counts.

	chrM	55	T	A,C,G,T,a,c,g,t,*,-:0,0,1,431,0,0,0,23,0,0
	chrM	56	A	A,C,G,T,a,c,g,t,*,-:453,0,0,0,19,0,0,0,0,0
	chrM	57	T	A,C,G,T,a,c,g,t,*,-:0,0,0,468,0,0,0,22,0,0



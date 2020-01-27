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
[Prinseq] (http://prinseq.sourceforge.net)
[STAR] (https://github.com/alexdobin/STAR)
[Samtools (1.0 or newer)] (http://www.htslib.org);
[Pigz] (http://zlib.net/pigz/);

Aliases:

SAMPLE_NAME = Name of sample used
REFERENCE_GENOME = Reference genome to align to
READ_LENGTH = Length of reads
MISMATCH = Number of mismatches to allow per read pair

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

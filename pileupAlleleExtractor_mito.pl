
#####HISTORIC#####
#Made by Alan Hodgkinson
#Name: pileupAlleleExtractor_mito.pl
#Script to create countfile for mito positions from BAM input.
#=========
#04/09/2015 : First draft of the script (AJH)
#=========



###SOFTWARE REQUIREMENTS (best in current path):

#Samtools (1.1 or higher)

###Use diagnostics:
use strict;
use Getopt::Long;
use Pod::Usage;

###Defaults and Inputs

my $bam_file = undef;
my $out_id = undef;
my $minQ = 13;
my $ref_fasta = "/away/ahodgkinson/Reference/human_g1k_v37chr.fasta";
my $samtools = "samtools";
my $nobaq = 0;

my $cmd_help;

GetOptions (
	    'Bam=s' => \$bam_file,
	    'Out=s' => \$out_id,
	    'MinQ=i' => \$minQ,
	    'RefFasta=s' => \$ref_fasta,
	    'SamtoolsPath=s' => \$samtools,
	    'noBaq' => \$nobaq,
	    'help' => \$cmd_help
) or exit(1);

pod2usage(1) if $cmd_help;

##File and software checks/warnings:

if (not defined $bam_file) {
  die "***You must provide a BAM file for the generation of a countfile.  Use the --Bam flag\n";
}
if (not defined $out_id) {
  die "***You must provide an out ID for the generation of a countfile.  Use the --Out flag\n";
}
if ($ref_fasta eq "/away/ahodgkinson/Reference/human_g1k_v37chr.fasta") {
  warn "--Default path to fasta file selected: /away/ahodgkinson/Reference/human_g1k_v37chr.fasta\n";
}
if ($samtools eq "samtools") {
  warn "--Default path to samtools selected\n";
}

my $outfile = $out_id.".MTcountfile.".$minQ.".txt";
if (-e $outfile) {
  die "***$outfile already exists: please rename the output with --Out or rename/move/delete the existing file\n";
}

##Main Program

my $r = int(rand(10000000));
my $tempfile = "temp_pileup".$r.".pu";

if (-e $tempfile) {
  die "***$tempfile already exists: Madness has happened and the random number generator has created a 1 in 10 million event just for you! Just re-run the program without changing anything and it should be fine next time (Oh, and try playing the lottery!)\n";
}

my $command = "$samtools mpileup -Q $minQ -d 10000000 -f $ref_fasta -r chrM $bam_file";
if ($nobaq) {
  $command = "$samtools mpileup -B -Q $minQ -d 10000000 -f $ref_fasta -r chrM $bam_file";
}

system ("$command > $tempfile");

open(INPUT, "$tempfile") || die "can't open $tempfile to read: $!\n";
open(OUTFILE, ">$outfile") || die "can't write to $outfile: $!\n";

my $prior = "";

while(my $line = <INPUT>){
  chomp $line;
  my @line_content = split(/[\s\t]/,$line);
  my $chr = $line_content[0];
  my $pos = $line_content[1];
  my $ref = uc($line_content[2]);
  my $dp = $line_content[3];
  my $bases = ($line_content[4]);
  my $score = $line_content[5];
  
  #Remove positional characters
  $bases =~ s/\$//g;  $bases =~ s/\^.//g; #these are start and end of read tag - as detailed above, mapping quality score also removed (following ^)
  my @bases = split(//,$bases);
  
  my %allele_count;
  
  #Store information on whether indels are on forward or reverse strand - the information for deletions is on the previous line (letter, size etc) - on the present line, a deletion is just *, so need to store whether it is on a rev or for strand from previous line
  my $fr_indels;
  for (my $pp=0;$pp<@bases;$pp++) {
    if ($bases[$pp] =~ /\d/) {
      if ($bases[$pp-1] eq "-") {
	if ($bases[$pp+1] =~ /[A-Z]/i) {
	  $fr_indels .= "F";
	}
	if ($bases[$pp+1] =~ /[a-z]/i) {
	  $fr_indels .= "R";
	}
      }
    }
  }
  
  #Treat the indels - i.e. revome the information for them from the strand.  As above this info is actually a shadow for the following position, so should be removed (not sure about insertions here, but we ignore them).
  my $indel_status = 0;
  my $indelSize = 0;
  foreach my $base (@bases){
    if($base eq '+' || $base eq '-'){
      $indel_status = 1;
      $base = "!";
    }
    elsif($indel_status && $base =~ m/\d/){
      if($indelSize == 0){
	$indelSize = $base;
      }
      else{
	$indelSize .=$base;
      }
      $base = "!";
    }
    elsif($indel_status && $indelSize >0 && $base !~ m/\d/){
      $base = "!"; #Transform in string again after and remove those - i.e. any information about indels
      $indelSize = $indelSize -1;
    }
  }
  
  my $bases2 = join('',@bases);
  $bases2 =~ s/!//g;
  
  #Phase2: REMOVE THE BASES WITH QUALITY LESS THAN $minQ
  @bases = split(//,$bases2);
  my @score = split(//,$score);	
  
  my @nuc = ('A','C','G','T','a','c','g','t','*','-');
  foreach my $t(@nuc){
    $allele_count{$t}=0;
  }
  
  my $indel_position = 0;
  my @indel_check = split (//,$prior);

  for (my $i=0;$i<scalar(@score);$i++){  #loop through bases and record allele called
    if ($bases[$i] eq '*') {  #change deletion symbol if on reverse strand
      if (length($prior)>0) { #if indel information is recorded
	if ($indel_check[$indel_position] eq "R") {
	  $bases[$i] = '-';
	}
	$indel_position++;
      }
    }
    
    if(ord($score[$i])-33<$minQ){}
    else{
      if( $bases[$i] =~ m/[\.,]/){
	if ($bases[$i] eq ".") {
	  if(!defined $allele_count{$ref}){
	    $allele_count{$ref} = 0;
	  }
	  $allele_count{$ref} +=1; 
	}
	if ($bases[$i] eq ",") {
	  if(!defined $allele_count{lc($ref)}){
	    $allele_count{lc($ref)} = 0;
	  }
	  $allele_count{lc($ref)} +=1; 
	}
      }
      elsif($bases[$i] =~ m/[\<\>]/){}
      else{
	if(!defined $allele_count{$bases[$i]}){
	  $allele_count{$bases[$i]} = 0;
	}
	$allele_count{$bases[$i]} +=1;
      }
    }
  }
  
  print OUTFILE "$chr\t$pos\t$ref\t"; #print data along with base calls
  my $index_tmp = 0;
  foreach my $keys (@nuc){
    if($index_tmp==0){
      print OUTFILE $keys;
    }
    else{
      print OUTFILE ",".$keys;
    }
    $index_tmp+=1;
  }
  print OUTFILE ":";
  $index_tmp = 0;
  foreach my $keys (@nuc){
    if($index_tmp==0){
      print OUTFILE $allele_count{$keys};
    }
    else{
      print OUTFILE ",".$allele_count{$keys};
    }
    $index_tmp+=1;
  }
  print OUTFILE "\n";
  $prior = $fr_indels; #Store previous lines indel information
}

close (INPUT);
close (OUTFILE);

system ("rm $tempfile");






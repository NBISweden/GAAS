#!/usr/bin/env perl
# Takes a fasta-file (usually a genomic or transcriptomic assembly) and checks for
# potential problems as well as calculates a few basic statistics.
#
# NBIS 2018
# jacques.dainat@nbis.se

use warnings;
use strict;
use Statistics::R;
use POSIX qw(strftime);
use File::Basename;
use Pod::Usage;
use Getopt::Long;
use Carp;
use IO::File;
use Bio::SeqIO;
use GAAS::GAAS;

my $header = get_gaas_header();
my $opt_infile;
my $opt_size=1000;
my $opt_dirRes;
my $opt_help = 0;


if ( !GetOptions( 'f|infile=s' => \$opt_infile,
                  's|size=i' => \$opt_size,
                  'o|out|output=s' => \$opt_dirRes,
                  'h|help!'        => \$opt_help ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( !( defined($opt_infile) ) ) {
    pod2usage( {
           -message => "$header\nMust specify at leat one parameter:\nfasta_purify.pl -f input.fa [-o Ouput_directory] ",
           -verbose => 0,
           -exitval => 1 } );
}


my ($ouputName,$path,$ext) = fileparse($opt_infile,qr/\.[^.]*/);
my $outstream = IO::File->new();
my $outstreamError;
my $outstreamFix;


if ( defined($opt_dirRes) ) {
	if (-d $opt_dirRes) {
		print "$opt_dirRes output directory already exits !\n";exit;
	}
	mkdir $opt_dirRes;

    $outstream->open( "$opt_dirRes/$ouputName"."_report.txt", 'w' ) or
        croak(sprintf( "Can not open '%s' for writing %s", "$opt_dirRes/$ouputName"."_report.txt", $! ));

    open(my $fh, '>', "$opt_dirRes/$ouputName"."_purified.fa") or die "Could not open file '$opt_dirRes/$ouputName"."_purified.fa' $!";
    $outstreamFix = Bio::SeqIO->new(-fh => $fh , -width  => 80, -format => 'fasta');

		$outstreamError = IO::File->new();
		$outstreamError->open( "$opt_dirRes/$ouputName"."_problem_found.log", 'w' ) or
				croak(sprintf( "Can not open '%s' for writing %s", "$opt_dirRes/$ouputName"."_problem_found.log", $! ));
}
else {
    $outstreamFix = Bio::SeqIO->new(-fh => \*STDOUT, -width  => 80, -format => 'fasta');

		$outstream->fdopen( fileno(STDERR), 'w' ) or
				        croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

# INPUT FILE
my $inseq = Bio::SeqIO->new(-file   => "<$opt_infile", -format => 'fasta');

my $nb_seq=0;
my $nb_seq_under_size=0;
my $nb_seq_with_IUPAC=0;
my $nb_seq_lowerCase=0;
my $nb_seq_Nstart=0;
my $nb_seq_NstartToShort=0; # too shor after removing leading Ns
my $nb_seq_Nend=0;
my $nb_seq_NendtToShort=0; # too shor after removing trailing Ns
my $nb_seq_NstartEnd=0;
my $nb_seq_pureN=0;

my $total_lowerCase=0;
my $total_Nremoved=0;
my $total_IUPAC=0;

my $size_in=0;
my $size_out=0;

my $size_internal_N = 10000;
my $nb_long_internal_N = 0;

#Calculate the statistics from the entries in the hash
while( my $seqObj = $inseq->next_seq() ) {
  $nb_seq++;
  $seqObj->desc("");
  my $sequence = $seqObj->seq();
  $size_in += length($sequence);

	# ---------- 1 FILTER SIZE --------------
	if (! is_long_enough($sequence,$opt_size) ){
		if ($outstreamError){
			print $outstreamError "seq_id:".$seqObj->id."\ttoo short\t".length($sequence)."\n";
		}
		next;
	}

 # ---------- 2 CLEAN Ns from extremities --------------
	# Pure Ns
	if ($sequence =~ /^N+$/){
		$nb_seq_pureN++;
		print $outstreamError "seq_id:".$seqObj->id."\tPure Ns sequence\tsize:".length($sequence)."\n";
		$nb_seq_under_size++;
		$total_Nremoved += length($sequence);
		next;
	}
	# LEADING Ns
	my $start_end_n=0;
  if ($sequence =~ /^N/){
		$nb_seq_Nstart++;
		$start_end_n=1;
		my $original_size = length($sequence);
		# remove leading Ns
		$sequence =~ s/^N+//g;
		$seqObj->seq($sequence);
		$total_Nremoved += $original_size-length($sequence);

		if ($outstreamError){
			print $outstreamError "seq_id:".$seqObj->id."\tNstart\tsize:".length($sequence)."\n";
		}

		if (! is_long_enough($sequence,$opt_size) ){
			if ($outstreamError){
				print $outstreamError "seq_id:".$seqObj->id."\ttoo short after removing leading Ns\tsize:".length($sequence)."\tquantity Ns removed:".$original_size-length($sequence)."\n";
			}
			$nb_seq_NstartToShort++;
			next;
		}
  }
	# TRAILING Ns
  if ($sequence =~ /N$/){
		$nb_seq_Nend++;
		$start_end_n=1;
		my $original_size = length($sequence);
		# remove trailing Ns
		$sequence =~ s/N+$//g;
		$seqObj->seq($sequence);
		$total_Nremoved += $original_size-length($sequence);

		if ($outstreamError){
			print $outstreamError "seq_id:".$seqObj->id."\tNend\tsize:".length($sequence)."\n";
		}

		if (! is_long_enough($sequence,$opt_size) ){
			if ($outstreamError){
				print $outstreamError "seq_id:".$seqObj->id."\ttoo short after removing trailing Ns\tsize:".length($sequence)."\tquantity Ns removed:".$original_size-length($sequence)."\n";
			}
			$nb_seq_NendtToShort++;
			next;
		}
  }

	if($start_end_n){$nb_seq_NstartEnd++;}

	# ---------- 3 DETECT long internal NNNNN --------------
  my @internal_ns = $sequence =~ /[ATGCURYSWKMBDHVatgcuryswkmbdhv]([Nn]+)[ATGCURYSWKMBDHVatgcuryswkmbdhv]/g;
	foreach my $internal_n ( @internal_ns ){
		if(length($internal_n) > $size_internal_N){
			$nb_long_internal_N++;
		}
	}

	# ---------- 4 CLEAN lower cases nucleotides --------------

	my $lowerCaseCount += ($sequence =~ tr/atgcunryswkmbdhv/atgcunryswkmbdhv/);
  $total_lowerCase += $lowerCaseCount;
  if($lowerCaseCount > 0){
		if ($outstreamError){
			print $outstreamError "seq_id:".$seqObj->id."\tcondains lowercase nucleotides\tquantity:".$lowerCaseCount."\n";
		}
		$nb_seq_lowerCase++;
		#Fix the lowercase issue
		$sequence =~ tr/atgcunryswkmbdhv/ATGCUNRYSWKMBDHV/;
		$seqObj->seq($sequence);
  }

  # ---------- IUPAC? --------------
  my $nb_IUPAC_here += ($sequence =~ tr/RYSWKMBDHVryswkmbdhv/RYSWKMBDHVryswkmbdhv/);
  $total_IUPAC += $nb_IUPAC_here;
  if($nb_IUPAC_here > 0){
	  if ($outstreamError){
      print $outstreamError "seq_id:".$seqObj->id."\tIUPAC found\tquantity:".$nb_IUPAC_here."\n";
		}
	  $nb_seq_with_IUPAC++;
  }

  $size_out += length($sequence);
  print $outstreamFix->write_seq( $seqObj );
}


my $seq_over_size = $nb_seq - $nb_seq_under_size ;
########################
#print out information #
my $stingToPrint=undef;
my $date = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
my $lineA= "-" x 80;
my $lineB= "-" x 78;
$lineB = "|".$lineB."|\n";
$stingToPrint .= $lineA."\n";
$stingToPrint .= "|".sizedPrint(basename($opt_infile),78)."|\n";
$stingToPrint .= "|".sizedPrint("Analysis launched the $date",78)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Comment",47)."|".sizedPrint("Number",15).
"|".sizedPrint("Task",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total sequence in",47,"L")."|".sizedPrint("$nb_seq",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total nucleotide in",47,"L")."|".sizedPrint("$size_in",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences pure Ns",47,"L")."|".sizedPrint("$nb_seq_pureN",15).
"|".sizedPrint("remove",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with leading and/or",47,"L")."|".sizedPrint("",15).
"|".sizedPrint("",14)."|\n";
$stingToPrint .= "|".sizedPrint("trailing Ns (pure N sequences not included)",47,"L")."|".sizedPrint("$nb_seq_NstartEnd",15).
"|".sizedPrint("trim",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total nb Ns removed (including pure)",47,"L")."|".sizedPrint("$total_Nremoved",15).
"|".sizedPrint("trim",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total nb of long internal N-regions >$size_internal_N",47,"L")."|".sizedPrint("",15).
"|".sizedPrint("",14)."|\n";
$stingToPrint .= "|".sizedPrint("This is problematic for Genemark",47,"L")."|".sizedPrint("$nb_long_internal_N",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("After removing leading and trailing Ns",47,"L")."|".sizedPrint("",15).
"|".sizedPrint("",14)."|\n";
$stingToPrint .= "|".sizedPrint("Nb of removed sequences (<$opt_size)",47,"L")."|".sizedPrint("$nb_seq_under_size",15).
"|".sizedPrint("remove",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with lowercase nucleotides",47,"L")."|".sizedPrint("$nb_seq_lowerCase",15).
"|".sizedPrint("fix",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of lowercase nucleotides",47,"L")."|".sizedPrint("$total_lowerCase",15).
"|".sizedPrint("fix",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with IUPAC nucleotides",47,"L")."|".sizedPrint("$nb_seq_with_IUPAC",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of IUPAC nucleotides",47,"L")."|".sizedPrint("$total_IUPAC",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total sequence out",47,"L")."|".sizedPrint("$seq_over_size",15).
"|".sizedPrint("inform",14)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Total nucleotide out",47,"L")."|".sizedPrint("$size_out",15).
"|".sizedPrint("inform",14)."|\n$lineB";
print $outstream "$stingToPrint";

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

sub is_long_enough{
	my ($sequence, $opt_size) = @_;

	if(length($sequence) < $opt_size){
		$nb_seq_under_size++;
		return undef;
	}
	return 1;
}

__END__
=head1 NAME

gaas_fasta_purify.pl

=head1 DESCRIPTION

This script aims to fix and inform the user about common problems found in an assembly,
in order to avoid problems during an annotation project.

	* lowercase nucleotides: There are considered as repeat by annotation tool,
	while most assemblers consider them as region with low coverage. Fixed by the tool.
	* Leading and trailing Ns sequences are forbidden for submission to ENA. Trimmed by the tool.
	* Long internal N regions (> 10000) is problematic in order to run Genemark. The tool inform the user.
	* IUPAC code might be problematic for some tools. e.g, an extra option is mandatory for MAKER. The tool inform the user.
	* Short sequences can be useless and time consuming for annotation tools. The tool remove them. Default 1000 bp.
  * Long sequence in single line fasta might raise issues some tools. The tool fold them (80 characters).

=head1 SYNOPSIS

    ./gaas_fasta_purify.pl --f=infile [--output=Directory]
    ./gaas_fasta_purify.pl --help

=head1 OPTIONS

=over 8

=item B<--f>, B<--infile> or B<-f>

Input fasta file containing DNA sequences.

=item B<--size> or B<-s>

Integer. Filter the sequence shorter to this size (in bp). Default: 1000

=item  B<--out>, B<--output> or B<-o>

[OPTIONAL] Output directory where diffrent output files will be written. If no output is specified, the result will written to STDOUT.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut

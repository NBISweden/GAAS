#!/usr/bin/env perl
# Takes a fasta-file (usually a genomic or transcriptomic assembly) and checks for
# potential problems as well as calculates a few basic statistics.
#
# NBIS 2018
# jacques.dainat@nbis.se

use warnings;
use strict;
use Try::Tiny;
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
my $nb_seq = 0;
my %seq_with_N_at_start_or_end;
my $nbseq_withLowerCase=0;
my $total_lowerCaseCount;
my $Ncount=0;
my $pureNseq=0;
my $totalcount=0;
my $totalcountOver1000=0;
my $totalcountOver10000=0;
my $gccount=0;
my $total_noNs=0;
my $nb_U =0;
my $nb_seq_with_U = 0;
my $nb_IUPAC = 0;
my $nb_seq_with_IUPAC = 0;
my @sequencelength=();
my @sequencelengthOver1000=();
my @sequencelengthOver10000=();

my $opt_infile;
my $opt_dirRes;
my $opt_help = 0;


if ( !GetOptions( 'f|infile=s' => \$opt_infile,
                  'o|out|output=s' => \$opt_dirRes,
                  'h|help!'         => \$opt_help ) )
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
           -message => "$header\nMust specify at leat one parameter:\nfasta_statistics.pl -f input.fa [-o Ouput_directory] ",
           -verbose => 0,
           -exitval => 1 } );
}


my ($ouputName,$path,$ext) = fileparse($opt_infile,qr/\.[^.]*/);
my $outstream = IO::File->new();

if ( defined($opt_dirRes) ) {
	if (-d $opt_dirRes) {
 		print "$opt_dirRes output directory already exits !\n";exit;
		}
	mkdir $opt_dirRes;

    $outstream->open( "$opt_dirRes/$ouputName"."_stat.txt", 'w' ) or
        croak(sprintf( "Can not open '%s' for writing %s", "$opt_dirRes/$ouputName"."_stat.txt", $! ));
}
else {
    $outstream->fdopen( fileno(STDOUT), 'w' ) or
	        croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

# INPUT FILE
my $inseq = Bio::SeqIO->new(-file   => "<$opt_infile", -format => 'fasta');


#Calculate the statistics from the entries in the hash
my $cp1kb=0;
my $cp10kb=0;
my $size_internal_N = 10000;
my $nb_long_internal_N = 0;
while( my $seqObj = $inseq->next_seq() ) {

  my $sequence = $seqObj->seq();
  $nb_seq++;

  #save sequence length for N50 calculation
  push @sequencelength, length $sequence;

  #Check for Ns at the beginning or end of sequence
  if ($sequence =~ /^N/){
    $seq_with_N_at_start_or_end{$seqObj->display_id}++;
  }
  if ($sequence =~ /N$/){
    $seq_with_N_at_start_or_end{$seqObj->display_id}++;
  }
  if (length($sequence) > 1000){
   	$cp1kb++;
	 #Count total size over 1000
	 $totalcountOver1000 += length $sequence;
 	  #save sequence length for N50 calculation
   	push @sequencelengthOver1000, length $sequence;

   	if (length($sequence) > 10000){
   		$cp10kb++;
	 	  #Count total size over 10000
   		  $totalcountOver10000 += length $sequence;
	 	  #save sequence length for N50 calculation
   		  push @sequencelengthOver10000, length $sequence;
   	}
  }

	# Count long internal NNNNN
  my @internal_ns = $sequence =~ /[ATGCURYSWKMBDHVatgcuryswkmbdhv]([Nn]+)[ATGCURYSWKMBDHVatgcuryswkmbdhv]/g;
	foreach my $internal_n ( @internal_ns ){
		$Ncount++;
		if(length($internal_n) > $size_internal_N){
			$nb_long_internal_N++;
		}
	}

  #Count GC
  $gccount += ($sequence =~ tr/gGcC/gGcC/);
  #Count size total
  $totalcount += length $sequence;
  #Count size with No Ns
  my $noNs=$sequence;
  $noNs =~ s/[Nn]//g; # remove Ns and ns
  $total_noNs += length $noNs;
  if(length $noNs == 0){
  	$pureNseq++;
  }
  #Count lowercase outside Ns
  my $lowerCaseCount += ($noNs =~ tr/atgcunryswkmbdhv/atgcunryswkmbdhv/);
  $total_lowerCaseCount += $lowerCaseCount;
  if($lowerCaseCount > 0){
	$nbseq_withLowerCase++;
  }
	#Count iupac
	my $nb_IUPAC_here += ($noNs =~ tr/RYSWKMBDHVryswkmbdhv/RYSWKMBDHVryswkmbdhv/);
	$nb_IUPAC += $nb_IUPAC_here;
	if($nb_IUPAC_here > 0){
		$nb_seq_with_IUPAC++;
	}
	#Count Uracile
	my $nb_U_here += ($noNs =~ tr/Uu/Uu/);
	$nb_U += $nb_U_here;
	if($nb_U_here > 0){
		$nb_seq_with_U++;
	}
}

# ------------

#Calculate some statistics
my $GCpercentage = ($gccount/$totalcount*100);
my $GCnoNs = ($gccount/$total_noNs*100);
my $totalNs =$totalcount-$total_noNs;

#################
# Calculate N50 #
@sequencelength = reverse sort { $a <=> $b } @sequencelength;
my $N50=$totalcount/2;
my $sum=0;
my $entry;
# copy of sequencelength to keep it intactfor R calculation purpose later
my @sequencelengthForN50Calcul=@sequencelength;

my $L50=0;
while ($sum < $N50){
  $entry = shift @sequencelengthForN50Calcul;
  $sum += $entry;
  $L50++;
}

#################
# Calculate N90 #
# copy of sequence length to keep it intact for R calculation purpose later
my @sequencelengthForN90Calcul=@sequencelength;
@sequencelengthForN90Calcul = reverse sort { $a <=> $b } @sequencelengthForN90Calcul;
my $NinetyPercGenomeSize=( 90*$totalcount / 100);
$sum=0;
my $N90;
my $L90=0;
while ($sum < $NinetyPercGenomeSize){
  $N90 = shift @sequencelengthForN90Calcul;
	$L90++;
  $sum += $N90;
}

###########################
# Calculate N50 over 1000 #
@sequencelengthOver1000 = reverse sort { $a <=> $b } @sequencelengthOver1000;
my $HalfGenomeSizeOver1000=$totalcountOver1000/2;
$sum=0;
my $N50over1000=0;
my $L50over1000=0;
while ($sum < $HalfGenomeSizeOver1000){
  $N50over1000 = shift @sequencelengthOver1000;
  $sum += $N50over1000;
	$L50over1000++;
}

############################
# Calculate N50 over 10000 #
@sequencelengthOver10000 = reverse sort { $a <=> $b } @sequencelengthOver10000;
my $HalfGenomeSizeOver10000=$totalcountOver10000/2;
$sum=0;
my $N50over10000=0;
my $L50over10000=0;
while ($sum < $HalfGenomeSizeOver10000){
  $N50over10000 = shift @sequencelengthOver10000;
	$L50over10000++;
  $sum += $N50over10000;
}


###########################
#print out the statistics #
my $stingToPrint=undef;
my $date = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
my $lineA= "-" x 80;
my $lineB= "-" x 78;
$lineB = "|".$lineB."|\n";
$stingToPrint .= $lineA."\n";
$stingToPrint .= "|".sizedPrint(basename($opt_infile),78)."|\n";
$stingToPrint .= "|".sizedPrint("Analysis launched the $date",78)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences",57,"L")."|".sizedPrint("$nb_seq",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences >1kb",57,"L")."|".sizedPrint("$cp1kb",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences >10kb",57,"L")."|".sizedPrint("$cp10kb",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of nucleotides (counting Ns)",57,"L")."|".sizedPrint("$totalcount",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of nucleotides U",57,"L")."|".sizedPrint("$nb_U",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with U nucleotides",57,"L")."|".sizedPrint("$nb_seq_with_U",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of IUPAC nucleotides",57,"L")."|".sizedPrint("$nb_IUPAC",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with IUPAC nucleotides",57,"L")."|".sizedPrint("$nb_seq_with_IUPAC",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of Ns",57,"L")."|".sizedPrint("$totalNs",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of internal N-regions (possibly links between contigs)",58,"L")."|".sizedPrint("$Ncount",19)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of long internal N-regions >$size_internal_N",57,"L")."|".sizedPrint("",20)."|\n";
$stingToPrint .= "|".sizedPrint("/!\\ This is problematic for Genemark",57,"L")."|".sizedPrint("$nb_long_internal_N",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of pure (only) N sequences",57,"L")."|".sizedPrint("$pureNseq",20)."|\n$lineB";
my $nb_seq_with_N_at_start_or_end= scalar keys %seq_with_N_at_start_or_end;
$stingToPrint .= "|".sizedPrint("Nb of sequences that begin or end with Ns",57,"L")."|".sizedPrint("$nb_seq_with_N_at_start_or_end",20)."|\n$lineB";
my $GC_content = sprintf("%.1f",$GCpercentage);
$stingToPrint .= "|".sizedPrint("GC-content (%)",57,"L")."|".sizedPrint("$GC_content",20)."|\n$lineB";
my $GC_content2 = sprintf("%.1f",$GCnoNs);
$stingToPrint .= "|".sizedPrint("GC-content not counting Ns(%)",57,"L")."|".sizedPrint("$GC_content2",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of sequences with lowercase nucleotides",57,"L")."|".sizedPrint("$nbseq_withLowerCase",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("Nb of lowercase nucleotides",57,"L")."|".sizedPrint("$total_lowerCaseCount",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("N50",57,"L")."|".sizedPrint("$entry",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("L50",57,"L")."|".sizedPrint("$L50",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("N90",57,"L")."|".sizedPrint("$N90",20)."|\n$lineB";
$stingToPrint .= "|".sizedPrint("L90",57,"L")."|".sizedPrint("$L90",20)."|\n$lineB";
#$stingToPrint .= "|".sizedPrint("N50 for sequences over 1000bp",57)."|".sizedPrint("$N50over1000",20)."|\n$lineB";
#$stingToPrint .= "|".sizedPrint("L50 for sequences over 1000bp",57)."|".sizedPrint("$L50over1000",20)."|\n$lineB";
#$stingToPrint .= "|".sizedPrint("N50 for sequeces over 10000bp",57)."|".sizedPrint("$N50over10000",20)."|\n$lineB";
#$stingToPrint .= "|".sizedPrint("L50 for sequeces over 10000bp",57)."|".sizedPrint("$L50over10000",20)."|\n$lineB";

print $outstream "$stingToPrint";

#######
#
# Plot
#
######
if($opt_dirRes){
	print $stingToPrint;
	if($nb_seq > 1){ #If only 1 seq we get error like: density.default(listValues) : need at least 2 points to select a bandwidth automatically
        # temporary file name
        my $tempFile1="dump.tmp";

		try {
			print "This result is saved in the <$opt_dirRes> directory along with plots in <pdf> format.\n";

			# write the data in temporary file
			open(FILE, ">$tempFile1") || die "Erreur E/S:$!\n";
			foreach my $size ( @sequencelength ) {
			  print FILE "$size\n";
			}
			close(FILE);

			my $ouputPlot=$opt_dirRes."/".$ouputName;

			# Calcul percentage of contig right and left to N50
			my $percentContigRightN50=($L50*100)/($#sequencelength+1);
			my $percentContigLeftN50=100-$percentContigRightN50;
			$percentContigRightN50=sprintf ("%0.2f",$percentContigRightN50)."%";
			$percentContigLeftN50=sprintf ("%0.2f",$percentContigLeftN50)."%";
			# Name of different outputs
			my $outputPlotLog=$ouputPlot."_PlotLog.pdf";
			my $outputPlotDensity=$ouputPlot."_PlotDensity.pdf";
			my $outputPlotHist=$ouputPlot."_PlotHisto.pdf";
			# calcul right and left position to write percentContigRightN50 and percentContigLeftN50
			my $biggestValue=shift @sequencelength;
			my $positionright=(5*$biggestValue)/100+$entry;
			my $positionleft=$entry-(5*$biggestValue)/100;
			# Tab=as.matrix(read.table("$tempFile1", sep="\t", he=T))
			# R object Declaration
			my $R = Statistics::R->new() or die "Problem with R : $!\n";

		## info ##
		#myhist$breaks contient la valeur minimale de chaque intervalle
		#myhist$mids contient la valeur au milieu de chaque intervalle
		#myhist$counts contient le nombre de valeurs situées dans cet intervalle
		#myhist$density contient la proportion de valeurs situées dans cet intervalle (autrement dit, tec.hist$counts / length (tec)).
		#print "tempFile1:$tempFile1 outputPlotLog:$outputPlotLog entry:$entry percentContigRightN50:$percentContigRightN50 percentContigLeftN50:$percentContigLeftN50\n";
		#print "outputPlotDensity:$outputPlotDensity outputPlotHist:$outputPlotHist positionright:$positionright positionleft:$positionleft\n";

			# R command
			 $R->send(
			     qq`
			  listValues=as.matrix(read.table("$tempFile1", sep="\t", he=F))
			  myhist<-hist(listValues)
			  legendToDisplay=paste("Number of value used : ",length(listValues))

			  pdf("$outputPlotLog")
			  plot(x = log(myhist\$mids), y = log(myhist\$counts), xlab="log(Contig size)", ylab="log(Frequency)", main="Size distribution of contigs")
			  abline(v =log($entry), col=2)
			  axisValues=par("usr")
			  ymax=(axisValues[4]*90)/100
			  ymax2=(axisValues[4]*85)/100
			  shiftFivePercent=(5*(axisValues[2]-axisValues[1]))/100
			  text(x = log($entry)+shiftFivePercent, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
			  text(x = log($entry)-shiftFivePercent, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
			  text(x = log($entry), y = ymax, paste("N50"), cex = 1, col = "red")
			  legend("topright", col=(1), lty=1, c(legendToDisplay))
			  dev.off()

			  pdf("$outputPlotDensity")
			  plot(density(listValues), xlab="Contig size", main="Size distribution of contigs")
			  abline(v =$entry, col=2)
			  axisValues=par("usr")
			  ymax=(axisValues[4]*90)/100
			  ymax2=(axisValues[4]*85)/100
			  text(x = $positionright, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
			  text(x = $positionleft, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
			  text(x = $entry, y = ymax, paste("N50"), cex = 1, col = "red")
			  legend("topright", col=(1), lty=1, c(legendToDisplay))
			  dev.off()

			  pdf("$outputPlotHist")
			  hist(listValues, xlab="Contig size",main="Size distribution of contigs")
			  abline(v =$entry, col=2)
			  axisValues=par("usr")
			  ymax=(axisValues[4]*90)/100
			  ymax2=(axisValues[4]*85)/100
			  text(x = $positionright, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
			  text(x = $positionleft, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
			  text(x = $entry, y = ymax, paste("N50"), cex = 1, col = "red")
			  legend("topright", col=(1), lty=1, c(legendToDisplay))
			  dev.off()`
			     );

			# Close the bridge
			$R->stopR();

			# Delete temporary file
			unlink $tempFile1;
			unlink "Rplots.pdf"; #created by " myhist<-hist(listValues)". I do not know how do differently...
		}
		catch{
			warn "caught error: $_";
			unlink $tempFile1;
		}
	}
	else{
		print "Not enough sequence to perform any plot.\n";
	}
}

__END__

=head1 NAME

gaas_fasta_statistics.pl

=head1 DESCRIPTION

Get some basic statistics about a nucleotide fasta file.
e.g Number of sequence, Number of nucleotide, N50, GC-content, etc.
It can also create R plots about contig size distribution.
The plots are performed only if an output is given.
This script is not designed for AA/Protein sequences.

=head1 SYNOPSIS

    ./gaas_fasta_statistics.pl --f=infile [--output=Directory]
    ./gaas_fasta_statistics.pl --help

=head1 OPTIONS

=over 8

=item B<--f>, B<--infile> or B<-f>

Input fasta file containing DNA sequences.

=item  B<--out>, B<--output> or B<-o>

[OPTIONAL] Output directory where diffrent output files will be written. If no output is specified, the result will written to STDOUT.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut

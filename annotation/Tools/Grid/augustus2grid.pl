#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use File::Basename;
use Bio::SeqIO;
use Cwd;
use Carp;
no strict qw(subs refs);
use GAAS::Grid::Bsub;
use GAAS::Grid::Sbatch;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outdir = "augustus_output";
my $genome = undef;
my $species = undef;
my $hints = undef;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my $grid="Slurm";
my $queue=undef;
my $help;

if ( ! GetOptions(
				    "hints=s" => \$hints,
						"species=s" => \$species,
						"genome|g=s" => \$genome,
						"grid=s"  => \$grid,
						"quiet|q!"  => \$quiet,
						"queue=s"  => \$queue,
				    "outdir=s" => \$outdir,
						"help|h!" => \$help ) )

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
        pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! defined($genome) or ! defined($species)){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\nInput genome fasta file (--genome) and augustus species hmm profile (--species)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# set grid option properly
my @grid_choice=('slurm','lsf','none');
$grid=lc($grid);
if (! grep( /^$grid/, @grid_choice ) ) {
  print "$grid is not a value accepted for grid parameter.";exit;
}
$grid= undef if lc($grid) eq 'none';

# .. Check that all binaries are available in $PATH

my @tools = ( "augustus" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

my $augustus_config_pathj = $ENV{'AUGUSTUS_CONFIG_PATH'} or die "AUGUSTUS_CONFIG_PATH is not set, aborting." ;

# .. Create output directory

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/augustus.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$genome", -format => 'fasta');

# .. and create chunks
msg("Creating chunks for GRID\n");

my $counter = 10000;
my $seq;

while( $seq = $inseq->next_seq() ) {
	$counter += 1;
  my $outfile = $outdir . "/seq_" . $counter . ".fa";
	my $cmd=undef;
	if($hints){
		$cmd = "augustus --species=$species --hintsfile=$hints --alternatives-from-evidence=true --gff3=on --extrinsicCfgFile=/references/software/augustus/config/extrinsic/extrinsic.E.cfg --uniqueGeneId=true $outfile > $outfile.augustus" ;
	}
	else{
		$cmd = "augustus --species=$species --alternatives-from-evidence=true --gff3=on --extrinsicCfgFile=/references/software/augustus/config/extrinsic/extrinsic.E.cfg --uniqueGeneId=true $outfile > $outfile.augustus" ;
	}
	push(@cmds,$cmd);
	my $seq_out = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
	$seq_out->write_seq($seq);
}

# Submit job chunks to grid
msg("submitting chunks\n");

if( $grid){
  msg("Sending $#cmds jobs to the grid\n");
  chomp(@cmds); # Remove empty indices
  # Submit job chunks to grid
  my $grid_runner;
  if ( $grid eq 'lsf'){
    $grid_runner = Bsub->new( cmds_list => \@cmds);
  }
  elsif( $grid eq 'slurm'){
    $grid_runner = Sbatch->new( cmds_list => \@cmds);
  }
  if($queue){$grid_runner->queue($queue)}
  $grid_runner->run();
}
else{
 	foreach my $command (@cmds){

 		system($command);

 		if ($? == -1) {
    			 print "failed to execute: $!\n";
		}
 		elsif ($? & 127) {
 		    printf "child died with signal %d, %s coredump\n",
 		    ($? & 127),  ($? & 128) ? 'with' : 'without';
 		}
 		else {
 		    printf "child exited with value %d\n", $? >> 8;
 		}
 	}
}

# ..Postprocessing here, like merging of output files

# --------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

# --------------------

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

# --------------------

sub check_bin {
	length(`which @_`) > 0 ? return 1 : return 0;
}

#----------------------------------------------------------------------

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}

__END__

=head1 NAME

gaas_augustus2grid.pl -

=head1 DESCRIPTION

Chunk input data to run multiple augustus jobs in parallel to grid

=head1 SYNOPSIS

    gaas_augustus2grid.pl -f fasta_file --species species_name
    gaas_augustus2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--genome> or B<-g>

The name of the genome fasta file.

=item B<--hints>

Augustus hints file (e.g Intron)

=item B<--species>

Species name for the hmm profile to use within Augustus

=item B<--queue>

If you want to define a particular queue to run the jobs

=item B<--grid>

Define which grid to use, Slurm, Lsf or None. Default = Slurm.

=item B<--quiet> or B<-q>

Quiet mode

=item B<--outdir> or B<-o>

The name of the output directory.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat

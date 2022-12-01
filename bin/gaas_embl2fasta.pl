#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outfile = undef;
my $embl = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "embl=s" => \$embl,
    "outfile|output|o|out|gb=s" => \$outfile))
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($embl)) ){
    pod2usage( {
           -message => "$header\nMissing the --embl argument",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $fasta_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $fasta_out= Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
}
else{
  $fasta_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta');
}

### Read gb input file.
my $embl_in = Bio::SeqIO->new(-file => $embl, -format => 'embl');


### MAIN ###

while( my $seq = $embl_in->next_seq) {
  $fasta_out->write_seq($seq)
}

__END__

=head1 NAME

gaas_embl2gb.pl

=head1 DESCRIPTION

The script take a EMBL file as input, and will translate into fasta format.

=head1 SYNOPSIS

    gaas_embl2fasta.pl --embl=infile.gff [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<-embl>

Input EMBL file that will be read

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gb>

Output fasta file. If no output file is specified, the output will be
written to STDOUT.

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

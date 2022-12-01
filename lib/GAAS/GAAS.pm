#!/usr/bin/perl -w

package GAAS::GAAS;

use strict;
use warnings;
use Exporter;

use GAAS::Utilities;

our $VERSION     = "v1.3.0";
our @ISA         = qw(Exporter);
our @EXPORT      = qw(get_gaas_header);
sub import {
  GAAS::GAAS->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
}

=head1 SYNOPSIS

  Meta package for conveniency.

=head1 DESCRIPTION



=head1 AUTHOR

    Jacques Dainat - jacques.dainat@nbis.se

=cut

# Provide meta information
sub get_gaas_header{

  my ($verbose) = @_;

  my $header = qq{
 ------------------------------------------------------------------------------
|   Genome Assembly Annotation Service (GAAS) - Version: $VERSION              |
|   https://github.com/NBISweden/AGAT                                          |
|   National Bioinformatics Infrastructure Sweden (NBIS) - www.nbis.se         |
 ------------------------------------------------------------------------------
  };

return $header;

}
1;

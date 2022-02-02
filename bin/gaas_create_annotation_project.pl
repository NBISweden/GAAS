#!/usr/bin/env perl

# file: gaas_create_annotation_project.pl
# Last modified: ons feb 02, 2022  02:24
# Sign: Johan Nylander

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path);

## Globals
my $annotation_root = "/projects/annotation";
my $default_name = "Genus_species-annotation_version-NBIS_ID";
my $logname = $ENV{LOGNAME};
my $time_stamp = localtime();
my $full_path = q{};
my @folders = (
    "abinitio",
    "customer_data",
    "genome",
    "maker",
    "organelles",
    "repeats",
    "rfam",
    "RNAseq");

## Options
my $name = undef;
my $assembly_version = undef;
my $id = undef;
my $path = undef;
my $version = undef;
my $help = undef;
my $man = undef;

GetOptions(
    "s|name=s" => \$name,
    "assembly-version=s" => \$assembly_version,
    "i|id=s" => \$id,
    "path=s" => \$path,
    "version" => \$version,
    "help" => \$help,
    "man" => \$man,
);

pod2usage(1) if ($help);
pod2usage(-exitval => 0, -verbose => 99, -sections => 'VERSION') if ($version);
pod2usage(-exitval => 0, -verbose => 2) if ($man);

## Check if path: then use this string, else use name
if ($path) {
    $full_path = $path;
}
elsif ($name) {
    $full_path = $annotation_root . "/" . $name;
}
else {
    $full_path = $annotation_root . "/" . $default_name;
}

## Check if folder exists, otherwise create it and the subfolders
if ( -e $full_path) {
    die "$0 WARNING:\nFolder $full_path already exists.\nCowardly refuses to overwrite. Exiting.\n";
}
else {
    for my $folder (@folders) {
        my $f = $full_path . "/" . $folder;
        #system("mkdir -p $f");
        make_path($f)
    }
}

## Create README.md
my $readme_file = $full_path . "/" . "README.md";

open my $FILE, ">", $readme_file or die "$0 ERROR: Could not open file \'$readme_file\' for writing: $!\n";

print $FILE "# README";
if ($name) {
    print $FILE " for $name";
}
if ($id) {
    print $FILE ", ID $id";
}
print $FILE "\n\n";
print $FILE "- Created: $time_stamp\n";
print $FILE "- Last modified: $time_stamp\n";
print $FILE "- Sign: $logname\n\n";
print $FILE "## Description:\n\n";
print $FILE "Text here\n";

close($FILE);

## Last check
if ( -e $full_path) {
    print STDERR "Created folders in project $full_path\n";
}
else {
    die "$0 ERROR: No project folder created\n";
}


__END__

=pod

=head1 NAME

gaas_create_annotation_project.pl - Create annotation project file hierarchy

=head1 VERSION

2.0

=head1 SYNOPSIS

gaas_create_annotation_project.pl [options]

 Options:
   --name               name of project
   --assembly-version   version string for genome assembly
   --id                 ID
   --help               brief help message
   --version            script version
   --man                full documentation


=head1 OPTIONS

Mandatory arguments to long options are mandatory for short options too.

Note that all arguments are optional.

=over 8

=item B<-n, --name=>I<string>

Name of the project

=item B<-s, --species=>I<string>

Same as B<--name>

=item B<-a, --assembly-version=>I<integer>

Version of the assembly used for the project

=item B<-i, --id=>I<string>

ID (e.g. NBIS redmine ID)

=item B<--help>

Print a brief help message and exits

=item B<--man>

Prints the manual page and exits

=back

=head1 DESCRIPTION

This script will create a project file hierarchy.
Without any arguments, this is the expected output

    Genus_species-annotation_version-NBIS_ID/
        |── abinitio/
        |── customer_data/
        |── genome/
        |── maker/
        |── organelles/
        |── repeats/
        |── rfam/
        |── RNAseq/
        |── README.md

The name of the parent folder can be changed using the B<--name> or the
B<--path> options (B<--path> will have precedence over B<--name>).

Default root of the annotation project is I</projects/annotation/>,
but this can be overridden by providing the full path using B<--path>.

=head1 EXAMPLES

    gaas_create_annotation_project.pl -n Apa_bpa
    gaas_create_annotation_project.pl -p /full/path/to/Apa_bpa
    gaas_create_annotation_project.pl -n Apa_bpa -a 1 -i 666

=cut


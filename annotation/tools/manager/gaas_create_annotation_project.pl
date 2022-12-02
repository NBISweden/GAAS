#!/usr/bin/env perl

# file: gaas_create_annotation_project.pl
# Last modified: tor maj 05, 2022  04:10
# Sign: Johan Nylander

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);

## Globals
my $annotation_root = "/projects/annotation"; # For testing, change this to full path to folder with write access
my $default_name = "Genus_species-annotation_version-NBIS_ID";
my $logname = $ENV{LOGNAME};
my $time_stamp = localtime();
my $full_path = q{};
my @folders = (
    "abinitio",
    "customer_data",
    "assembly",
    "maker/maker_evidence",
    "maker/maker_abinitio",
    "organelles",
    "public_data",
    "repeats",
    "rfam",
    "RNAseq");

## Options
my $copy_rnadata;
my $link_rnadata;
my $assembly_version;
my $man;
my $help;
my $id;
my $name;
my $path;
my $version;

GetOptions(
    "s|name=s" => \$name,
    "assembly-version=s" => \$assembly_version,
    "id=s" => \$id,
    "path=s" => \$path,
    "version" => \$version,
    "copy-rnadata=s" => \$copy_rnadata,
    "link-rnadata=s" => \$link_rnadata,
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
        make_path($f)
    }
}

## Copy or symlink rna data 
if ($copy_rnadata || $link_rnadata) {
    my $dest_folder = $full_path . "/" . "RNAseq";
    die "$0 ERROR:\nFolder $dest_folder does not exist.\n" unless ( -d $dest_folder);
    if ($copy_rnadata) {
        die "$0 ERROR:\nCan not find file $copy_rnadata.\n" unless ( -e $copy_rnadata);
        my $file = basename($copy_rnadata);
        my $copy = $dest_folder . "/" . $file;
        copy($file, $copy);
    }
    elsif ($link_rnadata) {
        die "$0 ERROR:\nCan not find file $link_rnadata.\n" unless ( -e $link_rnadata);
        my $file = basename($link_rnadata);
        my $symlink = $dest_folder . "/" . $file;
        my $absfile = abs_path($link_rnadata);
        symlink($absfile, $symlink);
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
   --copy-rnadata       copy rnadata file
   --link-rnadata       link rnadata file
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

=item B<-c, --copy-rnadata=>I<FILE>

Copy RNA data from I<FILE> to subfolder RNAseq

=item B<-l, --link-rnadata=>I<FILE>

Link (symbolic) RNA data from I<FILE> to subfolder RNAseq

=item B<--help>

Print a brief help message and exits

=item B<--man>

Prints the manual page and exits

=back

=head1 DESCRIPTION

This script will create a project file hierarchy.
Without any arguments, this is the default output

    Genus_species-annotation_version-NBIS_ID/
        L abinitio/
        L assembly/
        L customer_data/
        L maker/
            L maker_abinitio/
            L maker_evidence/
        L organelles/
        L public_data/
        L repeats/
        L rfam/
        L RNAseq/
        L README.md

The name of the parent folder can be changed using the B<--name> or the
B<--path> options (B<--path> will have precedence over B<--name>).

Default root of the annotation project is I</projects/annotation/>,
but this can be overridden by providing the full path using B<--path>.

=head1 EXAMPLES

    gaas_create_annotation_project.pl
    gaas_create_annotation_project.pl -n Blaps_mortisaga
    gaas_create_annotation_project.pl -p /full/path/to/Blaps_mortisaga
    gaas_create_annotation_project.pl -n Blaps_mortisaga -a 1 -i 666 -l /path/to/RNAdata.fa

=cut


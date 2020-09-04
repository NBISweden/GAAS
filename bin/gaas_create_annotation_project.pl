#!/usr/bin/env perl

#
#This script creates directories needed for a new project:
#
#NBIS 2018
#Nima Rafati
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Path;


my $logname = $ENV{ LOGNAME };
my $help;
my $fasta = undef;
my $version = undef;
my $species = undef;
my $annotation_root = "/projects/annotation";
#my $annotation_root = "~/test_annotation/";
my $usage = qq{
perl $0
	Getting help:
	[-help]
	Input:
	[-s species_name]
	[-g genome.fa]
	[-v assembly version]
	[-path path to project directory (default: /projects/annotation)]
};

GetOptions( 
"help" => \$help,
"s=s" => \$species,
"g=s" => \$fasta,
"v=i" => \$version,
"path=s" => \$annotation_root);

# Print Help and exit
if ($help) {
	print $usage;
	exit(0);
}
#check all parameters
if ( !( defined($fasta))){
	 pod2usage( { -message => "\nA fasta file for genome assembly must be provided (-g)\n$usage",
-verbose => 0,
-exitval => 2 } );
}
if ( !( defined($species))){
	 pod2usage( { -message => "\nA species name must be provided (-s)\n$usage",
-verbose => 0,
-exitval => 2 } );
}
if ( !( defined($version))){
	 pod2usage( { -message => "\nVersion of the genome assembly must be provided (-v)\n$usage",
-verbose => 0,
-exitval => 2 } );
}

my $project_path = "$annotation_root/$species/";
print "This project \" $species v.$version\" was created by \$logname = $logname on ". localtime().".\n You can find the working directory here:$project_path\n";



##Resources
make_dir($project_path, "Refseqs");
make_dir($project_path, "EST");
make_dir($project_path, "RNAseq");
make_dir($project_path, "Genome");
make_dir($project_path, "Mito");
make_dir($project_path, "Delivery");

##Ab initio
make_dir($project_path, "ab-initio");
make_dir("$project_path/ab-initio", "SNAP");
make_dir("$project_path/ab-initio", "Augustus");
make_dir("$project_path/ab-initio", "GeneMark_ET");
make_dir("$project_path/ab-initio", "GeneMark_EP");

##Maker
make_dir($project_path, "Maker");
make_dir("$project_path/Maker", "Evidence_build");
make_dir("$project_path/Maker", "gene_build");

##Transcriptome assembly
make_dir($project_path, "RNAseq_alignment");
make_dir($project_path, "Transcriptome_assembly");
make_dir("$project_path/Transcriptome_assembly", "Genome_guided");
make_dir("$project_path/Transcriptome_assembly", "Denovo");

##Prepare fasta file
prepare_fasta($fasta);

sub prepare_fasta {
	print "$fasta $project_path/Genome/";
	system("mv $fasta $project_path/Genome/");
	system("ln -s $project_path/Genome/$fasta $project_path/Genome/genome.fa");
}

sub make_dir {
#    my $directory = $_[0],"/",$_[1];
	my ($tmp_path, $tmp_dir) = @_;
	my $directory = "$tmp_path/$tmp_dir";
#	print "$directory"; <STDIN>;
	system("mkdir -p $directory");
#	mkdir $directory;
#    unless(mkdir $directory) {
#	unless(-e $directory or mkdir($directory,775)) {
#		"Unable to create $directory\n";
#	}
}

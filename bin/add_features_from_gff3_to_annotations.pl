#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin);
use lib $ENV{'CONDA_PREFIX'} . "/opt/jbrowse/src/perl5";
use JBlibs;

use File::Basename;
use lib dirname($0);

use Bio::GFF3::LowLevel::Parser;
use Getopt::Long qw(:config no_ignore_case bundling);
use IO::File;
use LWP::UserAgent;
use JSON;
use Data::Dumper;

my $gene_types_in = "gene";
my $pseudogene_types_in = "pseudogene";
my $transcript_types_in = "transcript|mRNA";
my $noncoding_rna_types = "ncRNA|rRNA|tRNA|snRNA|snoRNA|miRNA";
my $exon_types_in = "exon";
my $cds_types_in = "CDS";
my $ontology = "sequence";
my $gene_type_out = "gene";
my $pseudogene_type_out = "pseudogene";
my $mrna_type_out = "mRNA";
my $transcript_type_out = "transcript";
my $exon_type_out = "exon";
my $cds_type_out = "CDS";
my $sequence_alteration_types = "insertion|deletion|substitution";
my $annotation_track_prefix = "";
my $property_ontology = "feature_property";
my $comment_type_out = "comment";
my $property_type_out = "feature_property";
my $disable_cds_recalculation = 0;
my $organism;
my $in = \*STDIN;
my $username;
my $password;
my $url;
my $session_id;
my $name_attributes="Name";
my $test = 0;
my $use_name_for_feature = 0;

my $success_log;
my $error_log;

my $chunk_size = 100;

my $json = new JSON();
my $ua = new LWP::UserAgent();
$ua->timeout(3000); # 5 minutes

my %source_features = ();
my %skip_ids = ();

my %reserved_properties = (     symbol      => 1,
                                description => 1,
                                status      => 1
                          );
my %ignored_properties = (  owner               => 1,
                            date_creation       => 1,
                            date_last_modified  => 1
                         );

$| = 1;

parse_options();
process_gff();

sub parse_options {
    my $help;
    my $input_file;
    my $success_log_file;
    my $error_log_file;
    my $skip_file;
    GetOptions("input|i=s"          => \$input_file,
           "username|u=s"       => \$username,
           "password|p=s"       => \$password,
           "url|U=s"            => \$url,
           "gene_types_in|g=s"      => \$gene_types_in,
           "pseudogene_type_in|n=s" => \$pseudogene_types_in,
           "transcript_types_in|t=s"    => \$transcript_types_in,
           "exon_types_in|e=s"      => \$exon_types_in,
           "organism|o=s"      => \$organism,
           "cds_types_in|d=s"       => \$cds_types_in,
           "ontology|O=s"       => \$ontology,
           "gene_type_out|G=s"      => \$gene_type_out,
           "pseudogene_type_out|N=s" => \$pseudogene_type_out,
           "mrna_type_out|M=s"    => \$mrna_type_out,
           "transcript_type_out|T=s" => \$transcript_type_out,
           "exon_type_out|E=s"      => \$exon_type_out,
           "cds_type_out|D=s"       => \$cds_type_out,
           "property_ontolgy|R=s"   => \$property_ontology,
           "comment_type_out|C=s"   => \$comment_type_out,
           "property_type_out|S=s"  => \$property_type_out,
           "track_prefix|P=s"       => \$annotation_track_prefix,
           "disable_cds_recalculation|X"   => \$disable_cds_recalculation,
           "success_log|l=s"        => \$success_log_file,
           "error_log|L=s"      => \$error_log_file,
           "skip|s=s"           => \$skip_file,
           "test|x"           => \$test,
           "help|h"         => \$help,
           "name_attributes=s"   => \$name_attributes,
           "use_name_for_feature|a" => \$use_name_for_feature);

    print_usage() if $help;
    $in = new IO::File($input_file) or die "Error reading $input_file: $!\n"
        if $input_file;
    $success_log = new IO::File($success_log_file, "w") or die "Error writing success log $success_log_file: $!\n" if $success_log_file;
    $error_log = new IO::File($error_log_file, "w") or die "Error writing error log $error_log_file: $!\n" if $error_log_file;
    if ($skip_file) {
        my $fh = new IO::File($skip_file) or die "Error reading $skip_file: $!\n";
        while (my $id = <$fh>) {
            chomp $id;
            ++$skip_ids{$id};
        }
    }
	print_usage() if(!$username && !$password && !$url);
    die "Missing required parameter: username\n" if !$username;
    die "Missing required parameter: password\n" if !$password;
    die "Missing required parameter: url\n" if !$url;
    $url = "http://$url" if $url !~ m%https?://%;
}

sub print_usage {
    my $progname = basename($0);
    die << "END";

usage: $progname
    --url|-U <URL to Apollo instance>
    --username|-u <username>
    --password|-p <password>
    [--genes_types_in|-g <gene types for input>]
    [--pseudogene_type_in|-n <pseudogene type for input>]
    [--transcript_types_in|-t <transcript types for input>]
    [--exon_types_in|-e <exon types for input>]
    [--cds_types_in|-d <CDS types for input>]
    [--ontology|-O <ontology name used in server>]
    [--gene_type_out|-G <gene type used in server>]
    [--pseudogene_type_out|-N <pseudogene type used in server>]
    [--mrna_type_out|-T <mRNA type used in server>]
    [--transcript_type_out|-T <transcript type used in server>]
    [--exon_type_out|-E <exon type used in server>]
    [--cds_type_out|-D <CDS type used in server>]
    [--property_ontology|-R <property ontology name used in server>]
    [--comment_type_out|-C <comment type used in server>]
    [--property_type_out|-S <feature property type used in server>]
    [--track_prefix|-P <annotation track prefix>]
    [--disable_cds_recalculation|-X]
    [--input|-i <GFF3 file>]
    [--organism|-o <the name field used for your organism in WA2>]
    [--success_log|-l <success log file>]
    [--error_log|-L <error log file>]
    [--skip|-s <skip id file>]
    [--test|-x]
    [--help|-h]
    [--name_attributes <feature attribute to be used as name, first found used>]
    [--use_name_for_feature|-a]

    U: URL to Apollo instance
    u: username to access Apollo
    p: password to access Apollo
    g: string/regex to define the GFF3 types to treat as genes
       [default: "$gene_types_in"]
    n: string/regex to define the GFF3 types to treat as pseudogenes
       [default: "$pseudogene_types_in"]
    t: string/regex to define the GFF3 types to treat as transcripts
       [default: "$transcript_types_in"]
    e: string/regex to define the GFF3 types to treat as exons
       [default: "$exon_types_in"]
    d: string/regex to define the GFF3 types to treat as CDS
       [default: "$cds_types_in"]
    O: ontology name used in Apollo instance
       [default: "$ontology"]
    G: gene type used in Apollo instance
       [default: "$gene_type_out"]
    N: pseudogene type used in Apollo instance
       [default: "$pseudogene_type_out"]
    M: mRNA type used in Apollo instance
       [default: "$mrna_type_out"]
    T: transcript type used in Apollo instance
       [default: "$transcript_type_out"]
    E: exon type used in Apollo instance
       [default: "$exon_type_out"]
    D: CDS type used in Apollo instance
       [default: "$cds_type_out"]
    R: property ontology name used in Apollo instance
       [default: "$property_ontology"]
    C: comment type used in Apollo instance
       [default: "$comment_type_out"]
    S: feature property type used in Apollo instance
       [default: "$property_type_out"]
    P: annotation track prefix
       [default: "$annotation_track_prefix"]
    X: disable the recalculation of CDS by Apollo for protein coding features
    i: input GFF3 file
       [default: STDIN]
    o: organism common name in Apollo instance
    a: preserve feature names, as-is, in Apollo
    l: log file for ids that were successfully processed
    L: log file for ids that were erroneously processed
    s: file with ids to skip
    x: test mode
    h: this help screen
END
}

sub process_gff {
    my $gffio = Bio::GFF3::LowLevel::Parser->open($in);
    my $seq_ids_to_transcripts = {};
    my $seq_ids_to_genes = {};
    my $seq_ids_to_sequence_alterations = {};

    while (my $features = $gffio->next_item()) {
        next if ref $features ne "ARRAY";
        process_gff_entry($features, $seq_ids_to_genes, $seq_ids_to_transcripts, $seq_ids_to_sequence_alterations);
    }
    if ($test) {
        print "[";
        if (keys %{$seq_ids_to_genes} > 0) {
            print_features($seq_ids_to_genes, "addFeature");
        }
        if (keys %{$seq_ids_to_transcripts}) {
            if (keys %{$seq_ids_to_genes} > 0) {
                print ",";
            }
            print_features($seq_ids_to_transcripts, "addTranscript");
        }
        if (keys %{$seq_ids_to_sequence_alterations}) {
            if (keys %{$seq_ids_to_sequence_alterations} > 0) {
                print ",";
            }
            print_features($seq_ids_to_sequence_alterations, "addSequenceAlteration");
        }
        print "]";
    }
    else {
        write_features("addFeature", $seq_ids_to_genes);
        write_features("addTranscript", $seq_ids_to_transcripts);
        write_features("addSequenceAlteration", $seq_ids_to_sequence_alterations);
        if (!scalar(keys(%{$seq_ids_to_genes})) && !scalar(keys(%{$seq_ids_to_transcripts})) && !scalar(keys(%{$seq_ids_to_sequence_alterations}))) {
            print "No genes of type \"$gene_types_in\" or transcripts of type \"$transcript_types_in\" or sequence alterations of type \"$sequence_alteration_types\" found\n";
        }
    }

}

sub write_features {
    my $operation = shift;
    my $data = shift;
    while (my ($seq_id, $features) = each(%{$data})) {
        print "Processing $seq_id\n";
        my @to_be_written = ();
        my @ids = ();
        my @error = ();
        my $chunk = 0;
        for (my $i = 0; $i < scalar(@{$features}); ++$i) {
            my $feature = $features->[$i];
            push @to_be_written, $feature->[0];
            push @ids, $feature->[1];
            if (($i + 1) % $chunk_size == 0 || $i == scalar(@{$features}) - 1) {
                print "Processing chunk " . ++$chunk . ": ";
                my $res = send_request($seq_id, $operation, \@to_be_written);
                if ($res->{error}) {
                    print "error: $res->{error}\n";
                    for (my $j = 0; $j < scalar(@to_be_written); ++$j) {
                        push @error, [$to_be_written[$j], $ids[$j]];
                    }
                }
                else {
                    print "success\n";
                    foreach my $id (@ids) {
                        print $success_log "$id\n" if $success_log;
                    }
                }
                @to_be_written = ();
                @ids = ();
            }
        }
        if (scalar(@error)) {
            print "Processing individual features for failed chunks\n";
            foreach my $feature (@error) {
                print "Processing $feature->[1]: ";
                my $res = send_request($seq_id, $operation, [$feature->[0]]);
                if ($res->{error}) {
                    print "error: $res->{error}\n";
                    print $error_log "$feature->[1]\n" if $error_log;
                }
                else {
                    print "success\n";
                    print $success_log "$feature->[1]\n" if $success_log;
                }
            }
        }
    }
}

sub convert_feature {
    my $features = shift;
    my $type = shift;
    my $name_attributes = shift;
    my $start = get_start($features);
    my $end = get_end($features);
    my $strand = get_strand($features);
    my $name = get_name($features,$name_attributes);
    my $comments = get_comments($features);
    my $dbxrefs = get_dbxrefs($features);
    my $ontology_terms = get_ontology_terms($features);
    my $json_feature = {
        location => {
            fmin => $start - 1,
            fmax => $end,
            strand => $strand
        },
        type => {
            cv => {
                name => $ontology
            },
            name => $type
        }
    };
    if ($name) {
        $json_feature->{name} = $name;
    }
    my $json_properties = [];
    foreach my $comment (@{$comments}) {
        my $json_comment = {
            value => $comment,
            type => {
                cv => {
                    name => $property_ontology
                },
                name => $comment_type_out
            }
        };
        push @{$json_properties}, $json_comment;
    }
    foreach my $feature (@{$features}) {
        while (my ($tag, $values) = each(%{$feature->{attributes}})) {
            next if $tag =~ /^[A-Z]/;
            next if $ignored_properties{$tag};
            if ($tag =~ /symbol/) {
                $json_feature->{symbol} = $values->[0];
                next;
            }
            if ($tag =~ /description/) {
                $json_feature->{description} = $values->[0];
                next;
            }

            my $type = $tag;
            foreach my $value (@{$values}) {
                my $json_property = {
                    value => $value,
                    type => {
                        cv => {
                            name => $property_ontology
                        },
                        name => $type
                    }
                };
                push @{$json_properties}, $json_property;
            }
        }
        if (scalar(@{$json_properties})) {
            $json_feature->{properties} = $json_properties;
        }
    }
    my $json_dbxrefs = [];
    foreach my $dbxref (@{$dbxrefs}) {
        my ($db, $accession) = split ":", $dbxref;
        my $json_dbxref = {
            accession => $accession,
            db => {
                name => $db
            }
        };
        push @{$json_dbxrefs}, $json_dbxref;
    }
    foreach my $ontology_term (@{$ontology_terms}) {
        my @ontology_term_parts = split /:/, $ontology_term;
        my $json_ontology_term = {
            accession => $ontology_term_parts[1],
            db => {
                name => $ontology_term_parts[0]
            }
        };
        push @{$json_dbxrefs}, $json_ontology_term;
    }
    if (scalar(@{$dbxrefs})) {
        $json_feature->{dbxrefs} = $json_dbxrefs;
    }
    return $json_feature;
}

sub convert_sequence_alteration {
    my $features = shift;
    my $type = shift;
    my $name_attributes = shift;
    my $start = get_start($features);
    my $end = get_end($features);
    my $strand = get_strand($features);
    #my $name = get_name($features,$name_attributes);
    my $comments = get_comments($features);

    if ($type =~ /insertion/) {
        $start--;
        $end--;
    }
    else {
        $start--;
    }

    my $json_feature = {
        location => {
            fmin => $start,
            fmax => $end,
            strand => $strand
        },
        type => {
            cv => {
                name => $ontology
            },
            name => $type
        }
    };

    #if ($name) {
    #    $json_feature->{name} = $name;
    #}

    my $json_properties = [];
    foreach my $comment (@{$comments}) {
        my $json_comment = {
            value => $comment,
            type => {
                cv => {
                    name => $property_ontology
                },
                name => $comment_type_out
            }
        };
        push @{$json_properties}, $json_comment;
    }


    foreach my $feature(@{$features}) {
        while (my ($tag, $values) = each(%{$feature->{attributes}})) {
            next if $tag =~ /^[A-Z]/;
            next if $ignored_properties{$tag};
            if ($tag =~ /residues/) {
                $json_feature->{residues} = $values->[0];
                next;
            }

            my $type = $tag;
            foreach my $value (@{$values}) {
                my $json_property = {
                    value => $value,
                    type => {
                        cv => {
                            name => $property_ontology
                        },
                        name => $type
                    }
                };
                push @{$json_properties}, $json_property;
            }
        }
        if (scalar(@{$json_properties})) {
            $json_feature->{properties} = $json_properties;
        }
    }
    return $json_feature;
}

sub get_editor_url {
    my $cleanurl = $url =~ s,/*$,,;
    return "$url/annotationEditor";
}


sub send_request {
    my $seq_id = shift;
    my $operation = shift;
    my $features = shift;
    my $request = {
            track => "$annotation_track_prefix$seq_id",
            operation => $operation,
            username => $username,
            password => $password
    };
    if($organism) {
        $request->{organism}= $organism;
    }

    if ($features) {
        foreach my $feature (@{$features}) {
            push(@{$request->{features}}, $feature);
        }
    }
    my $req = new HTTP::Request();
    $req->method("POST");

    my $newurl = get_editor_url();
    $newurl = "$newurl/$operation";
    print "$newurl \n";
    $req->uri($newurl);
    $req->content($json->encode($request));
    my $res = $ua->request($req);
    if ($res->is_success()) {
        return $json->decode($res->content());
    }
    else {
        my $content;
        my $message;
        eval {
            $content = $json->decode($res->content());
        };
        if ($@) {
            $message = $res->status_line;
        }
        else {
            $message = $content->{error} ? $content->{error} : $res->status_line;
        }
        return { error => $message };
    }
}

sub get_name {
    my $features = shift;
    my $name_attributes = shift;
    my @attrs = split /\|/, $name_attributes;
    my $name;
    ATTR:
    foreach my $attr (@attrs) {
        $name = $features->[0]->{attributes}->{$attr}->[0];
        last ATTR unless !$name;
    }
    return $name;
}

sub get_comments {
    my $features = shift;
    my @comments;
    foreach my $feature (@{$features}) {
        my $values = $feature->{attributes}->{Note};
        push @comments, @{$values} if $values;
    }
    return \@comments;
}

sub get_dbxrefs {
    my $features = shift;
    my @dbxrefs;
    foreach my $feature (@{$features}) {
        my $values = $feature->{attributes}->{Dbxref};
        push @dbxrefs, @{$values} if $values;
    }
    return \@dbxrefs;
}

sub get_ontology_terms {
    my $features = shift;
    my @ontology_terms;
    foreach my $feature (@{$features}) {
        my $values = $feature->{attributes}->{Ontology_term};
        push @ontology_terms, @{$values} if $values;
    }
    return \@ontology_terms;
}

sub get_seq_id {
    my $features = shift;
    return $features->[0]->{seq_id};
}

sub get_start {
    my $features = shift;
    my $start = undef;
    foreach my $feature (@{$features}) {
        if (!defined $start || $feature->{start} < $start) {
            $start = $feature->{start};
        }
    }
    return int($start);
}

sub get_end {
    my $features = shift;
    my $end = undef;
    foreach my $feature (@{$features}) {
        if (!defined $end || $feature->{end} > $end) {
            $end = $feature->{end};
        }
    }
    return int($end);
}

sub get_strand {
    my $features = shift;
    my $strand = $features->[0]->{strand};
    if (!defined $strand) {
        $strand = "."
    }
    if ($strand eq "+") {
        return 1;
    }
    elsif ($strand eq "-") {
        return -1;
    }
    return 0;
}

sub get_id {
    my $features = shift;
    return $features->[0]->{attributes}->{ID}->[0];
}

sub get_type {
    my $features = shift;
    return $features->[0]->{type};
}

sub process_gene {
    my $features = shift;
    my $gene_type = get_type($features);
    my $gene;
    if ($gene_type =~ $pseudogene_types_in) {
        $gene = convert_feature($features, $pseudogene_type_out, $name_attributes);
    }
    else {
        $gene = convert_feature($features, $gene_type_out, $name_attributes);
    }
    if ($use_name_for_feature) {
        $gene->{use_name} = "true";
    }
    my $subfeatures = get_subfeatures($features);
    foreach my $subfeature (@{$subfeatures}) {
        my $type = get_type($subfeature);
        if ($type =~ /$transcript_types_in/) {
            my $transcript = process_transcript($subfeature);
            push(@{$gene->{children}}, $transcript) if $transcript;
        }
        elsif ($type =~ /$noncoding_rna_types/) {
            my $rna = process_rna($subfeature);
            push(@{$gene->{children}}, $rna) if $rna;
        }
    }
    return $gene;
}

sub process_transcript {
    my $features = shift;
    my $original_type = get_type($features);
    my $transcript;
    if ($original_type =~ /transcript/) {
        $transcript = convert_feature($features, $transcript_type_out, $name_attributes);
    }
    else {
        $transcript = convert_feature($features, $mrna_type_out, $name_attributes);
    }

    if ($use_name_for_feature) {
        $transcript->{use_name} = "true";
    }
    my $cds_feature = undef;
    my $subfeatures = get_subfeatures($features);
    foreach my $subfeature (@{$subfeatures}) {
        my $type = get_type($subfeature);
        if ($type =~ /$exon_types_in/) {
            my $exon = convert_feature($subfeature, $exon_type_out,$name_attributes);
            push(@{$transcript->{children}}, $exon);
        }
        elsif ($type =~ /$cds_types_in/) {
            my $start = get_start($subfeature);
            my $end = get_end($subfeature);
            if (!defined $cds_feature) {
                $cds_feature = $subfeature;
            }
            else {
                my $cds_start = get_start($cds_feature);
                my $cds_end = get_end($cds_feature);
                if ($start < $cds_start) {
                    $cds_feature->[0]->{start} = $start;
                }
                if ($end > $cds_end) {
                    $cds_feature->[0]->{end} = $end;
                }
            }
        }
    }
    if ($cds_feature) {
        my $cds = convert_feature($cds_feature, $cds_type_out, $name_attributes);
        push(@{$transcript->{children}}, $cds);
    }
    return $transcript;
}

sub process_rna {
    my $features = shift;
    my $type = get_type($features);
    my $rna = convert_feature($features, $type, $name_attributes);
    if ($use_name_for_feature) {
        $rna->{use_name} = "true";
    }
    my $subfeatures = get_subfeatures($features);
    foreach my $subfeature (@{$subfeatures}) {
        my $type = get_type($subfeature);
        if ($type =~ /$exon_types_in/) {
            my $exon = convert_feature($subfeature, $exon_type_out, $name_attributes);
            push(@{$rna->{children}}, $exon);
        }
    }
    return $rna;
}

sub get_subfeatures {
    my @subfeatures;
    my $features = shift;
    foreach my $feature (@{$features}) {
        push @subfeatures, @{$feature->{child_features}};
    }
    return \@subfeatures;
}

sub process_gff_entry {
    my $features = shift;
    my $seq_ids_to_genes = shift;
    my $seq_ids_to_transcripts = shift;
    my $seq_ids_to_sequence_alterations = shift;
    my $type = get_type($features);
    if ($type !~ /$gene_types_in/ && $type !~ /$transcript_types_in/) {
        my @subfeatures_array = get_subfeatures($features);
        if ($#subfeatures_array) {
            foreach my $subfeature (@subfeatures_array) {
                process_gff_entry($subfeature, $seq_ids_to_genes, $seq_ids_to_transcripts);
            }
        }
        else {
            my @json_features = process_singleton_feature($features);
            if (scalar(@json_features) != 0) {
                foreach my $json_feature (@json_features) {
                    my $id = get_id($features);
                    if ($skip_ids{$id}) {
                        print "Skipping $id\n";
                        next;
                    }
                    my $seq_id = get_seq_id($features);
                    if ($json_feature->{type}->{name} =~ /$sequence_alteration_types/) {
                        push(@{$seq_ids_to_sequence_alterations->{$seq_id}}, [$json_feature, $id]);
                    }
                    else {
                        push(@{$seq_ids_to_genes->{$seq_id}}, [$json_feature, $id]);
                    }
                }
            }
        }

    }
    else {
        my @json_features = @{process_feature($features)};
        if (scalar(@json_features) != 0) {
            foreach my $json_feature (@json_features) {
                my $id = get_id($features);
                if ($skip_ids{$id}) {
                    print "Skipping $id\n";
                    next;
                }
                my $seq_id = get_seq_id($features);
                if ($json_feature->{type}->{name} =~ /$gene_types_in/) {
                    push(@{$seq_ids_to_genes->{$seq_id}}, [$json_feature, $id]);
                }
                else {
                    push(@{$seq_ids_to_transcripts->{$seq_id}}, [$json_feature, $id]);
                }
            }
        }
    }
}

sub process_feature {
    my $features = shift;
    my $type = get_type($features);
    my $transcript_type = get_subfeature_type($features);

    if ($type =~ /$gene_types_in/) {
        if ($transcript_type =~ /mRNA/) {
            # process with mRNA at top-level and gene information as 'parent'
            return process_mrna($features);
        }
        elsif ($transcript_type =~ /$noncoding_rna_types/) {
            # process noncoding RNAs
            return [process_gene($features)];
        }
        else {
            # process top-level
            return [process_gene($features)];
        }
    }
    elsif ($type =~ /$transcript_types_in/) {
        return [process_transcript($features)];
    }
    return undef;
}

sub process_singleton_feature {
    my $features = shift;
    my $type = get_type($features);
    my $feature;
    if ($type =~ /$sequence_alteration_types/) {
        # process sequence alterations
        $feature = convert_sequence_alteration($features, $type, $name_attributes);
    }
    else {
        # process other singletons
        $feature = convert_feature($features, $type, $name_attributes);
        if ($use_name_for_feature) {
            $feature->{use_name} = "true";
        }
    }
    return $feature;
}

sub get_subfeature_type {
    my $features = shift;
    return ${features}->[0]->{child_features}->[0]->[0]->{type};
}

sub process_mrna {
    my $features = shift;
    my $gene = convert_feature($features, $gene_type_out, $name_attributes);
    my $subfeatures = get_subfeatures($features);
    my @processed_mrnas;

    foreach my $subfeature (@{$subfeatures}) {
        my $type = get_type($subfeature);
        if ($type =~ /mRNA/) {
            my $mrna_json_feature = convert_mrna_feature($subfeature, $gene);
            push(@processed_mrnas, $mrna_json_feature);
        }
    }
    return \@processed_mrnas;
}

sub convert_mrna_feature {
    my $features = shift;
    my $gene_json_feature = shift;
    my $cds_feature = undef;
    my $rtsc_json_feature = undef;
    my $mrna_json_feature = convert_feature($features, $mrna_type_out, $name_attributes);
    if ($disable_cds_recalculation) {
        $mrna_json_feature->{use_cds} = "true";
    }
    if ($use_name_for_feature) {
        $mrna_json_feature->{use_name} = "true";
    }
    $mrna_json_feature->{parent} = $gene_json_feature;
    my $subfeatures = get_subfeatures($features);

    foreach my $subfeature (@{$subfeatures}) {
        my $type = get_type($subfeature);
        if ($type =~ /$exon_types_in/) {
            my $exon = convert_feature($subfeature, $exon_type_out, $name_attributes);
            push(@{$mrna_json_feature->{children}}, $exon);
        }
        elsif ($type =~ /$cds_types_in/) {
            my $start = get_start($subfeature);
            my $end = get_end($subfeature);
            if (!defined $cds_feature) {
                $cds_feature = $subfeature;
            }
            else {
                my $cds_start = get_start($cds_feature);
                my $cds_end = get_end($cds_feature);
                if ($start < $cds_start) {
                    $cds_feature->[0]->{start} = $start;
                }
                if ($end > $cds_end) {
                    $cds_feature->[0]->{end} = $end;
                }
            }
            my $cds_subfeatures = get_subfeatures($cds_feature);
            if (scalar(@{$cds_subfeatures}) > 0) {
                $rtsc_json_feature = convert_feature(${$cds_subfeatures}[0], "stop_codon_read_through", $name_attributes);
            }
        }
        else {
            print "Ignoring unsupported sub-feature type: $type\n";
        }
    }

    if ($cds_feature) {
        my $cds_json_feature = convert_feature($cds_feature, $cds_type_out, $name_attributes);
        push(@{$cds_json_feature->{children}}, $rtsc_json_feature) if ($rtsc_json_feature);
        push(@{$mrna_json_feature->{children}}, $cds_json_feature);
    }

    return $mrna_json_feature;
}

sub print_features {
    my $data = shift;
    my $operation = shift;
    while (my ($seq_id, $features) = each(%{$data})) {
        my @to_be_written = ();
        my @ids = ();
        print "{\"$operation\": [";
        for (my $i = 0; $i < scalar(@${features}); ++$i) {
            my $feature = $features->[$i];
            print to_json(@$feature[0]);
            print "," if ($i != scalar(@$features) - 1);
        }
        print "]}";
    }
}

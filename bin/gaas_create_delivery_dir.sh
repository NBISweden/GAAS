#!/bin/bash

# JD 2015
#
#This script will soft-mask (lowerCase) a genome using a gff file as input
####



# Arguments and Paths
############################################################################
dirname=""

if (( $# != 1 )) ; then
        echo -e "No directory name given. By default we will call it <NBIS_delivery> "
	dirname="NBIS_delivery"
else
	if [[ $1 =~ \-?[hH]{1}elp$ || $1 =~ \-?[hH]{1}$ ]];then
		echo -e "This script prepare the five directories for genome annotation results as well as a Readme file.\nCommand:\n========\nscript.sh [dirName]\nBy default the directory name is called <delivery>"
		exit
	else
		dirname=$1
	fi
fi

if [[ -d $dirname ]];then
	echo "file <$dirname> already exists"; exit;
else
	mkdir -p "$dirname/metadata"
	mkdir -p "$dirname/summary"
	mkdir -p "$dirname/gff/repeats"
	mkdir -p "$dirname/gff/lift-over"
	mkdir -p "$dirname/gff/ncRNA"
	mkdir -p "$dirname/gff/gene_build"
	mkdir -p "$dirname/gff/transcript"
mkdir -p "$dirname/gff/proteins"
mkdir -p "$dirname/fasta/transcript"
  mkdir -p "$dirname/fasta/genome"
  mkdir -p "$dirname/fasta/repeat_library"
  mkdir -p "$dirname/miscellaneous/maker_evidence/parameters"
  mkdir -p "$dirname/miscellaneous/maker_evidence/tracks"
  mkdir -p "$dirname/miscellaneous/maker_evidence/busco"
  mkdir -p "$dirname/miscellaneous/maker_abinitio/parameters"
  mkdir -p "$dirname/miscellaneous/maker_abinitio/tracks"
  mkdir -p "$dirname/miscellaneous/maker_abinitio/busco"
  mkdir -p "$dirname/miscellaneous/abinitio_profiles"
  mkdir -p "$dirname/miscellaneous/completness/busco_assembly"
  mkdir -p "$dirname/miscellaneous/completness/busco_annotation"


	#write Readme
	echo -e "* NBIS Genome Annotation Platform team *" >> $dirname/ReadMe.txt
	mydate=$(date)
  	echo -e "<Species_Annotation_Report.docx> file => A complete report of work/process performed.
 ================================

 <fasta> directory => It contains data in FASTA format.
 =================
                    <cds.fasta> file => It contains cds of annotated genes.
                    <protein.fasta> file => It contains the proteins produced by the annotated genes.

                    <transcripts> sub directory
                    ===========
                    <Trinity_*_.fasta> is the output of our Trinity pipeline used as inout for Maker.
                    Trinity commando line in a txt file

                    <genome> sub directory
                    ===========
                    <genome_raw.fa> Raw genome/assembly received
                    <genome.fa> genome/assembly used for the annotation

                    <repeat_library> sub directory
                    ============
                    <repeat.fa> The library we created specifically for the species.

 <gff> directory => It contains data in gff3 format.
 ===============

                    <gene_builds> sub directory
                    ===========
                    <species_rcX.gff> Are the “release candidate” annotation for the species. The higher release candidate corresponds to the most successful achievement.

                    <lift-overs> sub directory (depending on customer demand)
                    ============
                    <referenceSpecies2studiedSpecies.gff> It is the lift-over of a reference species genome annotation (from Ensembl) on the genome under investigation.

                    <ncRNA> sub directory
                    =======
                    <ncRNA_rfam.gff> It contains annotation of non-coding element annotated using the Eucaryotes data from Rfam database.
			              <*_tRNA.gff> It contains transfert RNA annotated using tRNAscan.

                    <repeats> sub directory
	                  =======
                    <repeatmasker.gff> It contains repeats annotated thanks to repeatmasker.
		                <repeatrunner.gff> It contains repeats annotated thanks to repeatrunner

                    <proteins> sub directory
	                  =======
                    <protein2genome.gff> It contains proteins mapping the genome (using blast).

                    <transcripts> sub directory
                    ===========
                    <est2genome.gff> It contains transcripts obtained by mapping the trinity (fasta) file.
                    <est_gff_stringtie.gff> It contains transcripts obtained by mapping the RNAseq (Stringtie gff) file.
                    <*> Depending of RNAseq data types, and the most appropriate Transcriptome assembly method chosen.

 <metadata> directory
 ==========
                    <go.txt> Tab-delimited format file (2 columns) containing GO terms of functions retrieved for each transcript.
                    <pfam.txt> Tab-delimited format file (2 columns) containing pram terms of functions retrieved for each transcript.
                    <interpro.txt> Tab-delimited format file (2 columns) containing interpro terms of functions retrieved for each transcript.

 <summary> directory
 =========
                    <coding_gene.txt> General information/statistics about annotated coding genes
                    <tRNA.txt> General information about annotated tRNA
                    <repeats.txt> General statistics about the repeats

 <miscellaneous> directory
 =============
                    <maker_evidence> sub directory
                    ==============
                                    <maker_annotation_stat.txt> (Annotation statistics)
                                    <parameters> sub sub directory
                                     ==========
                                     <maker.opt>
                                     <maker.ctl>
                                     <maker.exe>
                                     <maker.evm>

                                    <tracks> sub sub directory
                                     ======
                                        <protein2genome.gff>
                                        <est2genome.gff>
                                    <busco> sub sub directory
                                     ======

                    <maker_abinitio> sub directory
                    ==============
                                    <maker_annotation_stat.txt> (Annotation statistics)
                                    <parameters> sub sub directory
                                     ==========
                                     <maker.opt>
                                     <maker.ctl>
                                     <maker.exe>
                                     <maker.evm>

                                    <tracks> sub sub directory
                                     ======
                                        <protein2genome.gff>
                                        <est2genome.gff>
                                        <augustus_masked.gff>
                                    <busco> sub sub directory
                                     ======

                    <abinitio_profiles> sub directory
                    ============
                                      <genemark.hmm>
                                      <augustus.hmm>
                                      <gene.gb/gene.gff> Gene set used to train the abinitio tools


                    <completness> sub directory
                    ============
                                      <busco_assembly> The busco result from assembly.
                                      <busco_annotation> The busco result from annotation.


 <bam> /!\ To give to the customer but its's no intended to be kept by NBIS except into Webapollo if set up. /!\
 =====
">> $dirname/ReadMe.txt


	echo -e "Well done ! Directories and ReadMe have been created. Think to check The ReadMe..."

fi

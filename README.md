[![Build Status](https://travis-ci.org/NBISweden/GAAS.svg?branch=master)](https://travis-ci.org/NBISweden/GAAS)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gaas/README.html)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/gaas.svg?style=flat)](https://anaconda.org/bioconda/gaas)
[![DOI](https://zenodo.org/badge/53501201.svg)](https://zenodo.org/badge/latestdoi/53501201)
[<img alt="docker_gaas" src="https://img.shields.io/badge/container-Docker-blue">](https://quay.io/repository/biocontainers/gaas)
[<img alt="singularity_gaas" src="https://img.shields.io/badge/container-Singularity-orange">](https://quay.io/repository/biocontainers/gaas)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gaas/badges/license.svg)](https://anaconda.org/bioconda/gaas)   
GAAS 
=========================================
<h2><em>G</em>enome <em>A</em>ssembly <em>A</em>nnotation <i>S</i>ervice (GAAS)</h2>  
Suite of tools related to Genome Assembly Annotation Service tasks.

[<img align="right" src="NBIS.png" width="200" height="100" />](https://nbis.se)

---------------------------

## Table of Contents

   * [What can GAAS do for you?](#what-can-gaas-do-for-you)
   * [Installation](#installation)  
       * [Using bioconda](using-bioconda)
          * [Install](#install)
          * [Update](#update)
          * [Uninstall](#uninstall)
       * [Old school](#old-school)
          * [Prerequisites](#prerequisites)
          * [Install](#install-1)
          * [Update](#update-1)
          * [Uninstall](#uninstall-1)
          * [Change to a specific version](#change-to-a-specific-version)
   * [Usage](#usage)
   * [Repository structure](#repository-structure)

---------------------------

## What can GAAS do for you?  

The repository contains mainly tools and knowledge related to bioinformatics and annotation the most often. To access and install the tools please follow the installation procedures below. For the knowledge you are invited to visit the [knowledge](annotation/knowledge) part of the repo or if you are looking specifically for genome assembly knowledge [The Genome Assembly Workshop Knowledge Base](https://github.com/NBISweden/workshop-genome_assembly/wiki).

## Installation

### Using conda

#### Install

  ```
  conda install -c bioconda gaas
  ```

#### Update

  ```
  conda update gaas
  ```

#### Uninstall
  ```
  conda uninstall gaas  
  ```

### Old school

#### Prerequisites
  * R
  * Perl
    Perl >= 5.8, and a list of perl modules that can be installed using cpan, cpanm or conda:

    * Install perl modules with cpanm
    ```
    cpanm install bioperl
    cpanm install Clone
    cpanm install Graph::Directed
    cpanm install LWP::UserAgent
    cpanm install Statistics::R
    cpanm install Sort::Naturally
    cpanm install File::Share
    cpanm install Moose
    cpanm install File::ShareDir::Install
    cpanm install Bio::DB::EUtilities
    ```
    * Install perl modules with conda

    ```
    conda env create -f conda_environment_GAAS.yml
    conda activate gaas
    ```

#### Install

  ```
  git clone https://github.com/NBISweden/GAAS.git # Clone GAAS
  cd GAAS                                         # move into GAAS folder
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```

<sup>*</sup>If dependencies are missing you can install them using cpan/cpanm or use conda and load the environment conda_environment_GAAS.yml

**Remark**: On MS Windows, instead of make you'd probably have to use dmake or nmake depending the toolchain you have.

#### Update  
From the folder where the repository is located.

  ```
  git pull                                        # Update to last GAAS
  perl Makefile.PL                                # Check all the dependencies<sup>1</sup>
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```

#### Change to a specific version
From the folder where the repository is located.  

  ```
  git pull                                        # Update the code
  git checkout v0.1.1                             # use version v0.1 (See releases tab for a list of available versions)
  perl Makefile.PL                                # Check all the dependencies<sup>1</sup>
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```

#### Uninstall

  ```
  perl uninstall_GAAS
  ```

## Usage

  ```
  script_name.pl -h
  ```    
  
---------------------------

## Repository structure

## [__annotation__](annotation)  
Annotation directory contains everything related to annotation side of the service.  

#### Shorcuts:  
   - [knowledge](annotation/knowledge)

   - [Genome annotation workshop](https://nbisweden.github.io/workshop-genome_annotation/)

   - [Tools](annotation/tools)  
     => All gff related work have been transplanted into [AGAT](https://github.com/NBISweden/AGAT) (11/2019)

   - [Pipelines](https://github.com/NBISweden/pipelines-nextflow)

## [__assembly__](assembly)  
Assembly directory contains development related to assembly side of the service.  

#### Shorcuts:  
   - [Genome assembly workshop](https://nbisweden.github.io/workshop-genome_assembly/)
   - [The Genome Assembly Workshop Knowledge Base](https://github.com/NBISweden/workshop-genome_assembly/wiki)


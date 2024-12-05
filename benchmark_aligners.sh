#!/bin/bash
set -e
#   Copyright 2020 Matthew Ralston
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# Author: Matt Ralston
# Date: 1/22/24
# Description:


source /home/matt/home/etc/optparse/optparse.bash

VERSION=0.0.1

# Miniconda3 3.12-24.9.2-0 and conda channels

# Use -c bioconda -c conda-forge -c biobuilds to install packages
# Optionally use environment.yml
export PATH=$PATH:$HOME/.pyenv/versions/miniconda3-3.12-24.9.2-0/bin


ART_ILLUMINA=/ffast/source_files/art_bin_MountRainier/art_illumina

TMPDIR=
ORIGINAL_WORKDIR=$(pwd)

# gmapper = SHRiMP
#ALIGNMENT_ROUTINES=(run_bowtie run_bowtie2 run_bbmap run_gmapper run_gsnap run_bwa_mem2 run_bwa_aln run_bwa_mem)


################################

#  OptParse | Options

################################
optparse.define short=n long=n_runs desc="Number of iterations to time each aligner run" variable=n
optparse.define short=f long=fasta desc="Fasta genome" variable=fasta
#optparse.define short=i long=fastq desc="Fastq file" variable=fastq
optparse.define short=p long=parallel desc="Core count for GNU parallel" variable=cores
#optparse.define short=o long=optional desc="Optional argument" variable=optional


source $( optparse.build )
################################
#  Required args and validations
################################
# if [ "$required" == "" ]; then
#     echo "Error: no input provided with -r|--required" 1>&2
#     exit 1
# fi

if [[ $# == 0 ]];
then
    echo "No arguments provided. See -?|--help for details."
    exit 1
fi


# Does file exist?
if [ ! -f $fasta ];
then
    echo "File doesn't exist at provided path: '$fasta'" 1>&2
    exit 1


elif ! [ $(echo $fasta | grep -E '^*.fa$') ] && ! [ $(echo $fasta | grep -E '^*.fasta$') ];
then
    echo $fasta
    echo $(echo $fasta | grep -E '^*.fasta$')
    echo "echo $fasta | grep -E '^*.fasta$'"
    ecode=$?
    echo "Exit code: $ecode"
    echo "Fasta sequence file does not have a .fasta or .fa file suffix. Possibly gzip/zlib compressed" 1>&2;
    exit 1
fi


# if [ ! -f $fastq ];
# then
#     echo "File doesn't exist at provided path: '$fastq'" 1>&2
#     exit 1
# elif ! [ $(echo $fastq | grep -E '^*.fq$') ] && ! [ $(echo $fastq | grep -E '^*.fastq$') ];
# then
#     echo "Fastq read file does not have a .fastq or .fq file suffix. Possibly gzip/zlib compressed" 1>&2;
#     exit 1
# fi


re='^[0-9]+$'
if ! [[ $n =~ $re ]];
then
    echo "error: argument '-n|--n_runs' was not a number" 1>&2
    exit 1
elif ! [[ $cores =~ $re ]];
then
    echo "error: argument '-p|--parallel' was not a number" 1>&2
    exit 1
fi


suffix="${fasta#*.}"
base_filename="${fasta%\.$suffix}" 




#####################################

#      F u n c t i o n s

#####################################

make_project_directory () {
    TMPDIR=$(mktemp -d aligner_benchmarking_temporary_directory.XXXXXX)
    full_fasta_filepath="$(readlink -f $fasta)"
    fasta_abs_dirpath="${full_fasta_filepath%/*}"
    fasta_filename="$(basename $full_fasta_filepath)"
    ln -s $fasta_abs_dirpath/$fasta_filename ./$TMPDIR/
    cd ./$TMPDIR
    echo "Benchmarking aligners in temporary project directory '$(pwd)'" 1>&2;
}

print_project_directory_and_exit() {
    echo 1>&2
    echo 1>&2
    echo 1>&2
    echo "================================" 1>&2
    echo "Completed aligner benchmarking. " 1>&2
    echo "Benchmark directory is '$TMPDIR'" 1>&2
    echo "================================" 1>&2
    echo 1>&2
    cd $ORIGINAL_WORKDIR
    echo 1>&2
    echo "DONE." 1>&2;
    echo 1>&2
}

generate_reads() {
    $ART_ILLUMINA -ss HS25 -i $fasta -p -l 150 -f 50 -m 200 -s 10 -o $base_filename
    # Additionally combine reads into single fastq file.
    cat ${base_filename}1.fq ${base_filename}2.fq > ${base_filename}_reads.fq

}


#####################################

#      I n d e x i n g

#####################################

run_bowtie_build() {
    # Generate a bowtie index
    bowtie-build $fasta $base_filename
}

run_bowtie2_build() {
    # Generate a bowtie2 index
    bowtie2-build $fasta $base_filename
}

run_bbmap_index() {
    # Generate a bbmap index
    bbmap.sh ref=$fasta
}

#####################################

#      B e n c h m a r k i n g

#####################################
run_bowtie() {
    # Generate alignment
    /usr/bin/time -a -o bowtie_timed.tsv -f "%x\t%e" bowtie -p $cores ${base_filename} -1 ${base_filename}1.fq -2 ${base_filename}2.fq -S ${base_filename}.sam
}

run_bowtie2() {
    # Run bowtie2
    /usr/bin/time -a -o bowtie2_timed.tsv -f "%x\t%e" bowtie2 -p $cores -x $base_filename -1 ${base_filename}1.fq -2 ${base_filename}2.fq -S ${base_filename}.sam
}

run_bbmap() {
    # Generate alignment
    /usr/bin/time -a -o bbmap_timed.tsv -f "%x\t%e" bbmap.sh in=${base_filename}_reads.fq threads=$cores out=${base_filename}.sam 
}

run_shrimp() {
    /usr/bin/time -a -o shrimp_timed.tsv -f "%x\t%e" gmapper -1 ${base_filename}1.fq -2 ${base_filename}2.fq $fasta
}

run_gsnap() {
    /usr/bin/time -a -o gsnap_timed.tsv -f "%x\t%e" gsnap -d $fasta -t $cores --novelsplicing=0 --format=sam ${base_filename}_reads.fq > ${base_filename}.sam
}





#####################################

#      M a i n

#####################################

make_project_directory


# Generate reads at 50x fold coverage, 150bp read length, 200bp insert size
generate_reads
# Index 
run_bowtie_build
run_bowtie2_build
run_bbmap_index

# Programs
run_bowtie
run_bowtie2
run_bbmap
#run_shrimp
#run_gsnap





print_project_directory_and_exit



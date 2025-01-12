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


export ART_ILLUMINA=/ffast/source_files/art_bin_MountRainier/art_illumina


FASTA=
GTF=

export ORIGINAL_WORKDIR=$(pwd)

# gmapper = SHRiMP
#ALIGNMENT_ROUTINES=(run_bowtie run_bowtie2 run_bbmap run_gmapper run_gsnap run_bwa_mem2 run_bwa_aln run_bwa_mem)


################################

#  OptParse | Options

################################
optparse.define short=n long=n_runs desc="Number of iterations to time each aligner run" variable=n
optparse.define short=f long=fasta desc="Fasta genome" variable=fasta
#optparse.define short=i long=fastq desc="Fastq file" variable=fastq
optparse.define short=p long=parallel desc="Core count for GNU parallel" variable=cores
optparse.define short=g long=gtf desc="GTF file for alignment benchmarking/assessment" variable=gtf
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


elif ! [ $(echo $fasta | grep -E '^*.fa$') ] && ! [ $(echo $fasta | grep -E '^*.fasta$') ] && ! [ $(echo $fasta | grep -E '^*.fna$') ];
then
    echo $fasta
    echo $(echo $fasta | grep -E '^*.fasta$')
    echo "echo $fasta | grep -E '^*.fasta$'"
    ecode=$?
    echo "Exit code: $ecode"
    echo "Fasta sequence file does not have a .fasta or .fa file suffix. Possibly gzip/zlib compressed. Please provide a decompressed fasta file to continue." 1>&2;
    exit 1
else
    export fasta
fi


if [ ! -f $gtf ];
then
    echo "GTF File doesn't exist at provided path: '$gtf'" 1>&2
    exit 1
else
    export gtf
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
    local -n project_dir=$1
    fasta=$2
    gtf=$3

    
    TMPDIR=$(mktemp --tmpdir=$TMPDIR -d aligner_benchmarking_temporary_directory.XXXXXX)
    full_fasta_filepath="$(readlink -f $fasta)"
    fasta_abs_dirpath="${full_fasta_filepath%/*}"
    fasta_filename="$(basename $full_fasta_filepath)"

    full_gtf_filepath="$(readlink -f $gtf)"
    gtf_abs_dirpath="${full_gtf_filepath%/*}"
    gtf_filename="$(basename $full_gtf_filepath)"
    
    cp $fasta_abs_dirpath/$fasta_filename $TMPDIR/
    FASTA=$TMPDIR/$fasta_filename

    cp $gtf_abs_dirpath/$gtf_filename $TMPDIR/
    GTF=$TMPDIR/$gtf_filename

    
    cd $TMPDIR
    echo "Benchmarking aligners in temporary project directory '$(pwd)'" 1>&2;

    project_dir=($TMPDIR $FASTA $GTF)
}

clean_project_directory() {
    echo 1>&2
    echo 1>&2
    echo 1>&2
    echo "Compressing fastq files from $TMPDIR..." 1>&2
    for f in $(/bin/ls $TMPDIR/*.fq);
    do
	gzip $f
    done
    echo "Removing old samfiles from $TMPDIR..." 1>&2
    #rm $TMPDIR/*.sam
    rm $TMPDIR/*.cleaned.bam
    echo 1>&2
    echo 1>&2
    echo 1>&2

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
    fasta=$1
    echo 1>&2
    echo 1>&2
    echo 1>&2
    echo "Generating reads with ART Illumina from fasta file '$fasta'..." 1>&2
    echo 1>&2
    echo 1>&2
    echo 1>&2

    $ART_ILLUMINA -ss HS25 -i $fasta -p -na -qL 20 -l 150 -c 10000000 -ir 0.01 -ir2 0.01 -dr 0.01 -dr2 0.01 -m 200 -s 10 -o "${fasta%.*}_"
    # Additionally combine reads into single fastq file.
    #wgsim -N 10000000 $FASTA ${FASTA%.*}_1.fq ${FASTA%.*}_2.fq
    cat ${fasta%.*}_1.fq ${fasta%.*}_2.fq > ${fasta%.*}_reads.fq
    rm *.aln
}


sort_alignment() {
    samfile=$1

    echo 1>&2
    echo 1>&2
    echo "Input samfile to sort: '$samfile'" 1>&2
    echo 1>&2
    echo 1>&2
    echo 1>&2
    
    output_bam="${samfile%.*}.sorted.bam"



    
    picard CleanSam -I $samfile -O "${samfile%.*}.cleaned.bam"
    #picard MarkDuplicates -I ${samfile%.*}.cleaned.bam -M ${samfile%.*}.markduplicates_metrics.txt -O ${samfile%.*}.marked.bam
    picard SortSam -I "${samfile%.*}.cleaned.bam" -O $output_bam --SORT_ORDER coordinate


    echo 1>&2
    echo 1>&2
    echo "Output sorted bamfile: '$output_bam'" 1>&2
    echo 1>&2
    echo 1>&2


    echo $output_bam
}


evaluate_alignment() {
    bamfile=$1
    #gtf=$2
    echo "==============================================" 1>&2
    echo 1>&2
    echo 1>&2
    echo "Evaluating alignment file ${bamfile}..." 1>&2
    echo 1>&2
    echo 1>&2
    echo "==============================================" 1>&2
    picard CollectAlignmentSummaryMetrics -I $bamfile -O ${bamfile%.*}.collectalignmentsummarymetrics.txt
    picard CollectInsertSizeMetrics -I $bamfile -H ${bamfile%.*}.histogram.txt -O ${bamfile%.*}.collectinsertsizemetrics.txt
    echo 1>&2
    echo 1>&2
    echo "Finished running Picard metrics" 1>&2
    echo 1>&2
    echo 1>&2
}




#####################################

#      I n d e x i n g

#####################################

run_bwa_mem_index(){
    fasta=$1
    base_filename="${fasta%.*}"
    # Generate a bwa-mem index
    bwa index $fasta -p $base_filename
}

run_bwa_mem2_index(){
    fasta=$1
    base_filename="${fasta%.*}"

    # Generate a bwa-mem2 index
    bwa-mem2 index -p $base_filename $fasta
}

run_bowtie_build() {
    fasta=$1
    base_filename="${fasta%.*}"
    
    # Generate a bowtie index
    bowtie-build $fasta $base_filename
}

run_bowtie2_build() {
    fasta=$1
    base_filename="${fasta%.*}"

    # Generate a bowtie2 index
    bowtie2-build $fasta $base_filename
}

run_bbmap_index() {
    fasta=$1

    # Generate a bbmap index
    bbmap.sh ref=$fasta
}

## Literally unusable
# run_gmap_index() {
#     Generate a GMap/GSnap index
#     echo 1&>2
#     echo 1&>2
#     echo "Genome is '${FASTA%.*}'" 1&>2
#     echo "FASTA is '$FASTA'" 1&>2
#     echo "GTF is '$GTF'" 1&>2
#     echo 1&>2
#     echo 1&>2
#     gmap_make_ref.sh ${FASTA%.*} $FASTA $GTF
#     gmap_build --dir=$(pwd) --db=${FASTA%.*} ${FASTA%.*} $FASTA
#     make gmapdb
# }




#####################################

#      B e n c h m a r k i n g

#####################################
run_bwa_mem() {
    fasta=$1
    base_filename="${fasta%.*}"
    output_bam="${base_filename}.bwa.sam"

    
    /usr/bin/time -a -o bwa_timed.tsv -f "%x\t%e" bwa mem -u -o $output_bam $base_filename ${base_filename}_1.fq ${base_filename}_2.fq
    echo 1>&2
    echo 1>&2
    echo "Completed running bwa..." 1>&2
    echo 1>&2
    echo 1>&2
    
    echo $output_bam
    
}

run_bwa_mem2() {
    fasta=$1
    base_filename="${fasta%.*}"
    output_bam="${base_filename}.bwa2.sam"
    
    /usr/bin/time -a -o bwa_timed.tsv -f "%x\t%e" bwa-mem2 mem -o $output_bam $base_filename ${base_filename}_1.fq ${base_filename}_2.fq
    echo 1>&2
    echo 1>&2
    echo "Completed running bwa2..." 1>&2
    echo 1>&2
    echo 1>&2
    
    echo $output_bam
    
}


run_bowtie() {
    fasta=$1
    base_filename="${fasta%.*}"
    output_bam="${base_filename}.bowtie.sam"
    
    # Generate alignment
    /usr/bin/time -a -o bowtie_timed.tsv -f "%x\t%e" bowtie $base_filename -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq -S $output_bam
    echo 1>&2
    echo 1>&2
    echo "Completed running bowtie..." 1>&2
    echo 1>&2
    echo 1>&2
    
    echo $output_bam
}

run_bowtie2() {
    fasta=$1
    base_filename="${fasta%.*}"
    output_bam="${base_filename}.bowtie2.sam"

    # Run bowtie2
    /usr/bin/time -a -o bowtie2_timed.tsv -f "%x\t%e" bowtie2 $cores -x $base_filename -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq -S $output_bam
    echo 1>&2
    echo 1>&2
    echo "Completed running bowtie2..." 1>&2
    echo 1>&2
    echo 1>&2

    echo $output_bam
}

run_bbmap() {
    fasta=$1
    base_filename="${fasta%.*}"


    output_bam="${base_filename}.bbmap.sam"
    
    # Generate alignment
    /usr/bin/time -a -o bbmap_timed.tsv -f "%x\t%e" bbmap.sh ref=$fasta in=${base_filename}_1.fq in2=${base_filename}_2.fq out=$output_bam
    echo 1>&2
    echo 1>&2
    echo "Completed running bbmap..." 1>&2
    echo 1>&2
    echo 1>&2


    
    echo $output_bam
}

run_shrimp() {
    fasta=$1
    base_filename="${fasta%.*}"

    output_bam="${base_filename}.shrimp.sam"
    
    # Generate alignment
    /usr/bin/time -a -o shrimp_timed.tsv -f "%x\t%e" gmapper -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq $fasta > $output_bam
    echo 1>&2
    echo 1>&2
    echo "Completed running SHRiMP..." 1>&2
    echo 1>&2
    echo 1>&2


    
    echo $output_bam
}






#####################################

#      M a i n

#####################################




main_routine() {
    local arr
    export TMPDIR=/ffast/scratch

    # Make a project directory
    make_project_directory arr $fasta $gtf

    read TMPDIR FASTA GTF <<< "${arr[*]}"


    
    # echo 1>&2
    # echo 1>&2
    # echo "Project directory is '$TMPDIR'" 1>&2
    # echo "Fasta file is '$FASTA'" 1>&2
    # echo "GTF file is '$GTF'" 1>&2
    # echo 1>&2
    # echo 1>&2
    
    # exit 1
    
    # # Generate reads at 50x fold coverage, 150bp read length, 200bp insert size
    generate_reads $FASTA
    # # Index
    run_bwa_index $FASTA # Not working properly
    run_bwa_mem2_index $FASTA
    run_bowtie_build $FASTA
    run_bowtie2_build $FASTA
    run_bbmap_index $FASTA
    # #run_gmap_index #Unusable. Repeated errors of 'unable to find genome directory x.'


    # # Programs
    bwa_mem_sam=$(run_bwa_mem $FASTA)
    bwa_mem2_sam=$(run_bwa_mem2 $FASTA)
    bowtie_sam=$(run_bowtie $FASTA)
    bowtie2_sam=$(run_bowtie2 $FASTA)
    bbmap_sam=$(run_bbmap $FASTA)
    # shrimp_sam=$(run_shrimp $FASTA)


    # # # Sort alignments
    bwa_mem_sorted_bam=$(sort_alignment $bwa_mem_sam)
    #bwa_mem2_sorted_bam=$(sort_alignment $bwa_mem2_sam)
    bowtie_sorted_bam=$(sort_alignment $bowtie_sam)
    #bowtie2_sorted_bam=$(sort_alignment $bowtie2_sam)
    #bbmap_sorted_bam=$(sort_alignment $bbmap_sam)
    # shrimp_sorted_bam=$(sort_alignment $shrimp_sam)


    # Evaluate alignments
    evaluate_alignment $bwa_mem_sorted_bam
    evaluate_alignment $bwa_mem2_sorted_bam
    evaluate_alignment $bowtie_sorted_bam
    evaluate_alignment $bowtie2_sorted_bam
    evaluate_alignment $bbmap_sorted_bam
    #evaluate_alignment $shrimp_sorted_bam

    # Cleanup
    clean_project_directory
    # Print project directory
    print_project_directory_and_exit
}

###################
# Export functions
###################

export -f main_routine
# Prerequisites
export -f make_project_directory
export -f generate_reads
# Build index
export -f run_bwa_mem_index
export -f run_bwa_mem2_index
export -f run_bowtie_build
export -f run_bowtie2_build
export -f run_bbmap_index
# Aligners
export -f run_bwa_mem
export -f run_bwa_mem2
export -f run_bowtie
export -f run_bowtie2
export -f run_bbmap
export -f run_shrimp
# Sort alignments
export -f sort_alignment
# Alignment evaluation
export -f evaluate_alignment
# Clean directory
export -f clean_project_directory
export -f print_project_directory_and_exit
# for i in $(seq 1 $n);
# do
#     main_routine
# done

parallel -j $cores 'main_routine' ::: $(seq 1 $n)

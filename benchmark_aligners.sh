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


FASTA=
GTF=

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


elif ! [ $(echo $fasta | grep -E '^*.fa$') ] && ! [ $(echo $fasta | grep -E '^*.fasta$') ];
then
    echo $fasta
    echo $(echo $fasta | grep -E '^*.fasta$')
    echo "echo $fasta | grep -E '^*.fasta$'"
    ecode=$?
    echo "Exit code: $ecode"
    echo "Fasta sequence file does not have a .fasta or .fa file suffix. Possibly gzip/zlib compressed. Please provide a decompressed fasta file to continue." 1>&2;
    exit 1
fi


if [ ! -f $gtf ];
then
    echo "GTF File doesn't exist at provided path: '$gtf'" 1>&2
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
    rm $TMPDIR/*.sam
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

    $ART_ILLUMINA -ss HS25 -i $fasta -p -qL 30 -l 150 -f 50 -m 200 -s 10 -o "${fasta%.*}_"
    # Additionally combine reads into single fastq file.
    #wgsim -N 10000000 $FASTA ${FASTA%.*}_1.fq ${FASTA%.*}_2.fq
    cat ${fasta%.*}_1.fq ${fasta%.*}_2.fq > ${fasta%.*}_reads.fq

}


sort_alignment() {
    samfile=$1
    output_bam=${samfile%.*}.sorted.bam
    picard CleanSam -I $samfile -O ${samfile%.*}.cleaned.bam
    #picard MarkDuplicates -I ${samfile%.*}.cleaned.bam -M ${samfile%.*}.markduplicates_metrics.txt -O ${samfile%.*}.marked.bam
    picard SortSam -I ${samfile%.*}.cleaned.bam -O $output_bam --SORT_ORDER coordinate
    echo $output_bam
}


evaluate_alignment() {
    samfile=$1
    gtf=$2
    picard CollectAlignmentSummaryMetrics -I $samfile -O ${samfile%.*}.collectalignmentsummarymetrics.txt
    picard CollectInsertSizeMetrics -I $samfile -H ${samfile%.*}.histogram.txt -O ${samfile%.*}.collectinsertsizemetrics.txt
    echo 1>&2
    echo 1>&2
    echo "Finished running Picard metrics" 1>&2
    echo 1>&2
    echo 1>&2
}




#####################################

#      I n d e x i n g

#####################################

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
run_bowtie() {
    fasta=$1
    base_filename="${fasta%.*}"

    
    # Generate alignment
    /usr/bin/time -a -o bowtie_timed.tsv -f "%x\t%e" bowtie -p $cores $base_filename -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq -S ${base_filename}.bowtie.sam
    echo 1>&2
    echo 1>&2
    echo "Completed running bowtie..." 1>&2
    echo 1>&2
    echo 1>&2
    
    echo "${FASTA%.*}.bowtie.sam"
}

run_bowtie2() {
    fasta=$1
    base_filename="${fasta%.*}"


    # Run bowtie2
    /usr/bin/time -a -o bowtie2_timed.tsv -f "%x\t%e" bowtie2 -p $cores -x $base_filename -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq -S ${base_filename}.bowtie2.sam
    echo 1>&2
    echo 1>&2
    echo "Completed running bowtie2..." 1>&2
    echo 1>&2
    echo 1>&2

    echo "${FASTA%.*}.bowtie2.sam"
}

run_bbmap() {
    fasta=$1
    base_filename="${fasta%.*}"

    # Generate alignment
    /usr/bin/time -a -o bbmap_timed.tsv -f "%x\t%e" bbmap.sh in=${base_filename}_reads.fq threads=$cores out=${base_filename}.bbmap.sam
    echo 1>&2
    echo 1>&2
    echo "Completed running bbmap..." 1>&2
    echo 1>&2
    echo 1>&2
    echo "${FASTA%.*}.bbmap.sam"
}

run_shrimp() {
    fasta=$1
    base_filename="${fasta%.*}"

    # Generate alignment
    /usr/bin/time -a -o shrimp_timed.tsv -f "%x\t%e" gmapper -1 ${base_filename}_1.fq -2 ${base_filename}_2.fq $fasta > ${base_filename}.shrimp.sam
    echo 1>&2
    echo 1>&2
    echo "Completed running SHRiMP..." 1>&2
    echo 1>&2
    echo 1>&2
    echo "${FASTA%.*}.shrimp.sam"
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
    
    # Generate reads at 50x fold coverage, 150bp read length, 200bp insert size
    generate_reads $FASTA
    # Index 
    run_bowtie_build $FASTA
    run_bowtie2_build $FASTA
    run_bbmap_index $FASTA
    # #run_gmap_index #Unusable. Repeated errors of 'unable to find genome directory x.'


    # # Programs
    bowtie_sam=$(run_bowtie $FASTA)
    bowtie2_sam=$(run_bowtie2 $FASTA)
    bbmap_sam=$(run_bbmap $FASTA)
    shrimp_sam=$(run_shrimp $FASTA)


    # # Sort alignments
    bowtie_sorted_bam=$(sort_alignment $bowtie_sam)
    bowtie2_sorted_bam=$(sort_alignment $bowtie2_sam)
    bbmap_sorted_bam=$(sort_alignment $bbmap_sam)
    shrimp_sorted_bam=$(sort_alignment $shrimp_sam)

    # Evaluate alignments
    # evaluate_alignment $bowtie_sorted_bam
    # evaluate_alignment $bowtie2_sorted_bam
    # evaluate_alignment $bbmap_sorted_bam
    # evaluate_alignment $shrimp_sorted_bam

    # Cleanup
    clean_project_directory
    # Print project directory
    print_project_directory_and_exit
}

for i in $(seq 1 $n);
do
    main_routine
done

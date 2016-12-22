#!/bin/bash

. /broad/tools/scripts/useuse
use .hdfview-2.9
use Java-1.8
use .r-3.1.3-gatk-only

gdir=$1
echo $gdir

sample=$2
echo $sample

bamfile=`cat $gdir/bams.tsv | grep "${sample}" | cut -f 2`
echo $bamfile

hetsfile=$gdir/out/$sample.hets
echo $hetsfile

# intervals=/xchip/cga/reference/hg19/whole_exome_illumina_coding_v1_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list
intervals=/xchip/cga/reference/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list
echo $intervals

mem=4
jlp=/broad/software/free/Linux/redhat_6_x86_64/pkgs/hdfview_2.9/HDFView/lib/linux/
gatk_jar=/xchip/scarter/dmccabe/gatk-protected-0cdd5b1/build/libs/gatk-protected.jar

java -Xmx${mem}g -Djava.library.path=${jlp} -jar ${gatk_jar} GetBayesianHetCoverage \
  --normal $bamfile --normalHets $hetsfile \
  --reference /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
  --snpIntervals $intervals \
  --help false --version false --verbosity DEBUG --QUIET false

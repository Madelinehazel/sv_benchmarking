#!/bin/bash
####################################################################################################
#   submit_gnomad_comparison.sh:
#   creates script for merging reference and test SV vcfs by reciprocal overlap for all gnomAD samples
#   produces output bed for each sample indicating which tools called each SV
#   Madeline Couse, 2020-05-07

#   input:
#	    ref_vcf_path = path to reference vcfs, e.g. /hpf/largeprojects/ccmbio_ephemeral/SV.Validation/data/vcf/
#	    tools = comma seperated list of tool names, e.g. wham,delly,lumpy,manta,cnvkit,cnvnator,erds
#     tool_vcf_paths = comma seperated list of paths to tool vcfs, in same order as tool argument
#     rec_overlap = reciprocal overlap fraction to compare SVs, e.g. 0.5
#     output_path = path to output directory
#   output:
#     produces tsv for each sample indicating which tools called each SV

####################################################################################################

samples=(HG00512  HG00513)
ref_vcf_path=$1
tools=$2
tool_vcf_paths="$(echo $3 | tr ',' ' ')"
rec_overlap=$4
output_path=$5


if [ ! -d $output_path ]
then
  mkdir $output_path
fi


for i in ${samples[@]}; do
    ref_vcf=$ref_vcf_path/nstd152.GRCh37.variant_call.$i.reformat.vcf

    #make array of paths to output vcfs for this particular sample
    declare -a sample_vcfs
    for tool in ${tool_vcf_paths[@]}; do
      vcf="$(echo ${tool}${i}*)"
      sample_vcfs+=( $vcf )
    done

    sample_vcfs="$(echo "${sample_vcfs[@]}")"
    #make overlap script for each sample
    echo -e "source activate /hpf/largeprojects/ccmbio/mcouse/SV_comparison/envs \
    \npython /hpf/largeprojects/ccmbio/mcouse/SV_comparison/scripts/wrapper_script/sv_benchmarking/benchmark/compare_sv_vcfs.py \
    -i $ref_vcf $sample_vcfs \
    -r_overlap 0.5 -o ${output_path}/${i}_gnomad_overlap" > ${output_path}/gnomad_overlap_analysis.$i.sh

    #qsub -l vmem=16g,mem=16g,walltime=12:00:00 \
    #-N $i.gnomad.overlap \
    #-o $i.gnomad.overlap.o \
    #-e $i.gnomad.overlap.e \
    #${output_path}/gnomad_overlap_analysis.$i.sh

done

#Next, call benchmark/aggregate_sample_sv_comparisons.py; dependent on the above finishing
#Then, make summary statistics and plots

# sv_benchmarking
Scripts for benchmarking SV calls against gnomAD reference set (HG00512  HG00513  HG00514  HG00731  HG00732  HG00733  NA19238  NA19239 NA19240)

## 1. Intersect test SVs with reference SVs

submit_gnomad_comparison.sh: this creates scripts for merging reference and test SV vcfs by reciprocal overlap for all gnomAD samples and submits these as jobs to the hpf. It produces an  output bed for each sample indicating which tools called each SV.

This script takes as input:
* ref_vcf_path ($1) = path to reference vcfs
* tools ($2) = comma seperated list of tool names, e.g. wham,delly,lumpy,manta,cnvkit,cnvnator,erds
* tool_vcf_paths ($3) = comma seperated list of paths to tool vcfs, in same order as tool argument
* rec_overlap ($4) = reciprocal overlap fraction to compare SVs, e.g. 0.5
* output_path ($5) = path to output directory

An example script to run submit_gnomad_comparison might look like this: 

```

vcf_dir= ~/test_data/vcf

sh sv_benchmarking/submit_gnomad_comparison.sh ~/ref_data/vcf/ \
  manta,lumpy,wham,cnvkit,delly,erds,cnvnator,melt_alu,melt_line1,melt_SVA \
  ${vcf_dir}/manta/,${vcf_dir}/lumpy/,${vcf_dir}/wham/,/${vcf_dir}/cnvkit/,${vcf_dir}/delly/,${vcf_dir}/erds/,${vcf_dir}/cnvnator/,${vcf_dir}/melt/ALU/,${vcf_dir}/melt/LINE1/,${vcf_dir}/melt/SVA/ \
  0.5 \
  ~/overlap_output
```

## 2. Generate metrics 

sv_benchmarking/aggregate_sample_sv_comparisons.py: this script aggregates reports from overlap between SV truth sets and test sets for all nine gnomAD samples and produces metrics (recall, precision, F1) in csv and plot form. 

This script takes as input:
* overlap_dir =  path to directory containing output bed for each sample generated in step 1.
* tools  = list of tool names, e.g. wham delly lumpy manta cnvkit cnvnator erds
* output = path to output directory for summary statistics and plots


An example script to run submit_gnomad_comparison might look like this:


```
python sv_benchmarking/benchmark/aggregate_sample_sv_comparisons.py -overlap_dir ~/overlap_output \
  -tools manta lumpy wham cnvkit delly erds cnvnator melt_alu melt_line1 melt_SVA \
  -output ~/overlap_stats_output
```
   

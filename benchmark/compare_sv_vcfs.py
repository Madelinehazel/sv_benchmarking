#!/usr/bin/env python3
#copied from Arun's directory, /hpf/largeprojects/ccmbio_ephemeral/SV.Validation/scripts/compare.sv.calls2.py

import argparse
import allel
import numpy as np
import pandas as pd
import os
from pybedtools import BedTool



def list2string(col):
    return col.apply(lambda x: ', '.join(x) if isinstance(x, list) else x)

def parse_vcfs(vcfs, ins=False):

    ann_dfs = []

    for vcf in vcfs:
        print("Parsing %s" % vcf)

        vcf_dict = allel.read_vcf(vcf, ['*'])
        name = os.path.basename(vcf).replace(".vcf", "")

        for key in list(vcf_dict.keys()):
            if key not in ['variants/CHROM','variants/POS', 'variants/ALT','variants/END', 'variants/SVTYPE', 'variants/SVLEN']:
                vcf_dict.pop(key)


        vcf_dict['variants/ALT'] = [alt[0] for alt in vcf_dict['variants/ALT']]
        df = pd.DataFrame(vcf_dict)
        df['SAMPLE'] = name
        sample_list.append(name)

        # throw out BND rows
        df = df[df['variants/SVTYPE'] != "BND"]

        if 'variants/END' not in df.columns:
            df['variants/END'] = np.nan
        elif 'variants/SVLEN' not in df.columns:
            # Delly VCF
            df['variants/SVLEN'] = df['variants/END'] - df['variants/POS']
        elif 'variants/END' not in df.columns and 'variants/SVLEN' not in df.columns:
            raise Exception("VCF does not have a SVLEN or END field. Are you sure this is a structural variant VCF?")

        # compute END for ALU-insertion,LINE1-insertion,SVA-insertion SVTYPE rows
        no_end_svtypes = (df['variants/SVTYPE'] == "ALU") | (df['variants/SVTYPE'] == "LINE1") | (df['variants/SVTYPE'] == "SVA") | (df['variants/SVTYPE'] == "INS")
        if len(df.loc[no_end_svtypes, :]) != 0:
            df.loc[no_end_svtypes, ['variants/END']] = df['variants/POS']

        # workaround for START > END so BedTool doesn't freak out with MalformedBedError
        # START > END is the case for TRV, INV, DEL
        start_gt_end = df['variants/END'] < df['variants/POS']
        if len(df.loc[start_gt_end, :]) > 0:
            df.loc[start_gt_end, ['variants/END','variants/POS']] = df.loc[start_gt_end, ['variants/POS','variants/END']].values
            df['variants/POS'] = df['variants/POS'].astype(int)
            df['variants/END'] = df['variants/END'].astype(int)
            df = df.drop_duplicates()

        # compute SVLEN for empty records (SVTYPE=CNV)
        df.loc[(df['variants/SVLEN'] == '') | (df['variants/SVLEN'] == '.'), 'variants/SVLEN'] = df['variants/END'] - df['variants/POS']

        # typecast potential floats to ints (due to implicit typecasting)
        df[["variants/POS", "variants/END", 'variants/SVLEN']] = df[["variants/POS", "variants/END", 'variants/SVLEN']].astype(int)

        #print(df)

        ann_dfs.append(df)

    ann_df = pd.concat(ann_dfs)
    ann_df.columns = ann_df.columns.str.replace('variants/', '')
    ann_df = ann_df[["CHROM", "POS", "END", "ALT", "SVTYPE", "SAMPLE"]]

    if ins:
        ann_df = ann_df[ann_df['SVTYPE'].str.contains('INS')]
    else:
        ann_df = ann_df[ann_df['SVTYPE'].str.contains('DEL|DUP|INV')]
        ann_df = ann_df.drop(columns=['ALT'])

    return BedTool(list(ann_df.itertuples(index=False, name=None)))

def group_sv(bedtool, reciprocal_overlap, ins=False):

    print('Identifying equivalent structural variant calls using a reciprocal overlap of %f' % reciprocal_overlap)
    already_grouped_intervals = set()
    columns = ['CHROM', 'POS', 'END', 'SVTYPE', 'N_SAMPLES']
    columns.extend(sample_list)
    columns.extend(["%s_SV_DETAILS" % s for s in sample_list])

    master_df = pd.DataFrame(columns=columns).set_index(['CHROM', 'POS', 'END', 'SVTYPE'])

    if not ins:
        #for DEL, DUP, INV, compare with bedtools intersect reciprocal overlap
        for l in bedtool.intersect(bedtool, wa=True, wb=True, F=reciprocal_overlap, f=reciprocal_overlap):

            ref_chr, ref_start, ref_end, ref_svtype, ref_name, \
            samp_chr, samp_start, samp_end, samp_svtype, samp_name = l

            ref_interval = (ref_chr, ref_start, ref_end, ref_svtype)
            samp_interval = (samp_chr, samp_start, samp_end, samp_svtype, samp_name)
            if (samp_interval not in already_grouped_intervals) and (ref_svtype == samp_svtype):
                add_interval(master_df, ref_interval, samp_interval)
                already_grouped_intervals.add(samp_interval)

    else:
        for l in bedtool.window(bedtool, w=200):

            ref_chr, ref_start, ref_end, ref_alt, ref_svtype, ref_name, \
            samp_chr, samp_start, samp_end, samp_alt, samp_svtype, samp_name = l

            ref_interval = (ref_chr, ref_start, ref_end, ref_alt)
            samp_interval = (samp_chr, samp_start, samp_end, samp_alt, samp_name)
            if (samp_interval not in already_grouped_intervals):
                add_interval(master_df, ref_interval, samp_interval)
                already_grouped_intervals.add(samp_interval)

    return master_df.sort_index()

def add_interval(master_df, ref_interval, samp_interval):
    samp_chr, samp_start, samp_end, samp_svtype, samp_name = samp_interval

    #Get reference to row
    if ref_interval not in master_df.index:
        #make new row
        master_df.loc[ref_interval, :] = np.nan
        row = master_df.loc[ref_interval, :]
        row['N_SAMPLES'] = 0

        for name in sample_list:
            row[name] = 0
            row["%s_SV_DETAILS" % name] = []
    else:
        row = master_df.loc[ref_interval, :]

    #Set values for row
    if row[samp_name] == 0:
        row['N_SAMPLES'] += 1
        row[samp_name] = 1

    row["%s_SV_DETAILS" % samp_name].append('{}:{}-{}:{}'.format(samp_chr, samp_start, samp_end, samp_svtype))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates a report to compare structural variant calls between multiple callers (and a truth dataset) through reciprocal overlap')
    parser.add_argument('-i', nargs='+', help='VCF files containing structural variants, e.g. -i 180_230.vcf 180_231.vcf', required=True)
    parser.add_argument('-r_overlap', type=float, help='Percent reciprocal overlap required to consider 2 calls to be equivalent in decimal', required=True)
    parser.add_argument('-o', type=str, help='TSV output file prefix')
    args = parser.parse_args()

    vcfs = args.i
    sample_list = []
    del_dup_inv = parse_vcfs(vcfs)
    del_dup_inv_df = group_sv(del_dup_inv, args.r_overlap)
    del_dup_inv_df.to_csv('%s_DEL_DUP_INV.tsv'%args.o, sep="\t")

    sample_list = []
    ins = parse_vcfs(vcfs,ins=True)
    ins_df = group_sv(ins, args.r_overlap, ins=True)
    ins_df.to_csv('%s_INS.tsv'%args.o, sep="\t")

import argparse
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from sigfig import round
from sklearn.metrics import confusion_matrix, recall_score, precision_score, f1_score


def make_colnames(tools):
    col_names = ['CHROM', 'POS', 'END', 'SVTYPE', 'N_SAMPLES', 'ref']
    details = [tool + '_SV_DETAILS' for tool in tools]
    col_names = col_names + tools + ['ref_SV_details'] + details
    return col_names


def make_overlap_df(overlaps, col_names, size_bins):
    overlap_list = []
    for f in overlaps:
        df = pd.read_csv(f, sep='\t')
        df.columns = col_names
        sample = f.split('/')[-1].split('_')[0]
        df['sample'] = [sample]*len(df)
        overlap_list.append(df)
    overlap_df = pd.concat(overlap_list)
    overlap_df['LENGTH'] = overlap_df['END'] - overlap_df['POS']
    overlap_df['SIZE'] = overlap_df.apply(
        lambda row: get_size_bin(row['LENGTH'], size_bins), axis=1)
    return overlap_df


def get_size_bin(len, size_bins):
    for key, value in size_bins.items():
        if value[0] <= len <= value[1]:
            size_bin = key
    return size_bin


def get_prediction(overlap_df, tools):
    # make array of aggregate predictions from all tools
    # i.e. if 0 in all tools, prediction is 0
    # if 1 in at least one tool, prediction is 1
    predicted_agg = []
    for index, row in overlap_df.iterrows():
        agg = max([row[tool] for tool in tools])
        predicted_agg.append(agg)
    overlap_df['aggregate_prediction'] = predicted_agg
    return overlap_df


def get_summary_stats(overlap_df, tools, bin_order):
    metric_list = []
    all_tools = tools + ['aggregate_prediction']
    for svtype in ['DEL', 'DUP', 'INV']:
        for tool in all_tools:
            for size in bin_order:
                sv_true = overlap_df[(overlap_df['SVTYPE'] == svtype)
                                     & (overlap_df['SIZE'] == size)]['ref']
                sv_pred = overlap_df[(overlap_df['SVTYPE'] == svtype)
                                     & (overlap_df['SIZE'] == size)][tool]
                precision = precision_score(sv_true, sv_pred)
                recall = recall_score(sv_true, sv_pred)
                f1 = f1_score(sv_true, sv_pred)
                metric_list.append([svtype, tool, size, precision, recall, f1])
    metric_df = pd.DataFrame(metric_list, columns=[
                             'SVTYPE', 'tool', 'size', 'precision', 'recall', 'F1'])
    return metric_df


def get_summary_stats_ins(overlap_df_ins, tools):
    metric_list = []
    all_tools = tools + ['aggregate_prediction']
    for svtype in set(overlap_df_ins['SVTYPE']):
        for tool in all_tools:
            sv_true = overlap_df_ins[(overlap_df_ins['SVTYPE'] == svtype)]['ref']
            sv_pred = overlap_df_ins[(overlap_df_ins['SVTYPE'] == svtype)][tool]
            precision = precision_score(sv_true, sv_pred)
            recall = recall_score(sv_true, sv_pred)
            f1 = f1_score(sv_true, sv_pred)
            metric_list.append([svtype, tool, precision, recall, f1])
    metric_df = pd.DataFrame(metric_list, columns=['SVTYPE', 'tool', 'precision', 'recall', 'F1'])
    return metric_df


def plot_metrics(metric_df, metric, output_path):
    # plots recall, precision, F1 for DEL, DUP, INV
    sns.set(font_scale=1)
    g = sns.catplot(x="size", y=metric, col="SVTYPE",
                    data=metric_df, kind="point",
                    order=bin_order, hue="tool", ci=None, margin_titles=True, aspect=1.5)
    g.set_xticklabels(rotation=45, horizontalalignment='left')
    plt.savefig('%s/%s_DEL_DUP_by_svtype_and_size' %
                (output_path, metric), dpi=300, bbox_inches='tight')


def plot_overall_metrics(metric_df, output_path, ins=False):
    sns.set(font_scale=0.5)
    if not ins:
        svtypes = ['DEL', 'DUP', 'INV']
    else:
        svtypes = set(metric_df['SVTYPE'])
    for svtype in svtypes:
        fig, ax = plt.subplots(1, 3, figsize=(7, 4), sharey=True)
        subset_svtype = metric_df[metric_df['SVTYPE'] == svtype]
        i = 0
        for metric in ['recall', 'precision', 'F1']:
            sns.barplot(x='tool', y=metric, data=subset_svtype, ax=ax[i], ci=None, color='c')
            plt.ylim(0, 1)
            for tick in ax[i].get_xticklabels():
                tick.set_rotation(90)
                tick.set_horizontalalignment('center')
            i += 1
        plt.tight_layout()
        # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        svtype = svtype.strip('<').strip('>').replace(':', '_')
        plt.savefig('%s/%s_overall_metrics' % (output_path, svtype), dpi=300, bbox_inches='tight')


def output_TP_FP_FN(overlap_df, output_path, ins=False):
    if not ins:
        svtypes = ['DEL', 'DUP', 'INV']
    else:
        svtypes = set(metric_df['SVTYPE'])
    for svtype in svtypes:
        svtype = svtype.strip('<').strip('>').replace(':', '_')
        TP = overlap_df[(overlap_df['ref'] == 1)
                        & (overlap_df['aggregate_prediction'] == 1)
                        & (overlap_df['SVTYPE'] == svtype)]
        TP.to_csv('%s/%s_TP_SVs.tsv' % (output_path, svtype))
        FN = overlap_df[(overlap_df['ref'] == 1)
                        & (overlap_df['aggregate_prediction'] == 0)
                        & (overlap_df['SVTYPE'] == svtype)]
        FN.to_csv('%s/%s_FN_SVs.tsv' % (output_path, svtype))
        FP = overlap_df[(overlap_df['ref'] == 0)
                        & (overlap_df['aggregate_prediction'] == 1)
                        & (overlap_df['SVTYPE'] == svtype)]
        FP.to_csv('%s/%s_FP_SVs.tsv' % (output_path, svtype))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Aggregates reports from overlap between SV truth sets and test sets for all nine gnomAD samples and produces metrics')
    parser.add_argument('-overlap_dir', type=str,
                        help='Path to directory containing overlap.tsv files', required=True)
    parser.add_argument('-tools', type=str, nargs='+',
                        help='tools used to call svs', required=True)
    parser.add_argument('-output_path', type=str, help='Path to output directory', required=True)
    args = parser.parse_args()

    print(args.tools)

    overlaps = glob.glob('%s/*overlap_DEL_DUP_INV.tsv' % args.overlap_dir)
    size_bins = {"0-50bp": [0, 49], "50-100bp": [50, 99], "100-300bp": [100, 299], "300-500bp": [300, 499], "500-1000bp": [500, 999], "1-2kb": [1000, 1999], "2-5kb": [2000, 4999],
                 "5-10kb": [5000, 9999], "10kb+": [10000, 9999999999]}
    col_names = make_colnames(args.tools)
    print(col_names)

    # make dataframe aggregating SVs from all nine gnomAD samples
    overlap_df_non_ins = make_overlap_df(overlaps, col_names, size_bins)

    # add aggregate prediction (0 or 1) across all tools for each SV
    overlap_df_non_ins = get_prediction(overlap_df_non_ins, args.tools)

    # get summary stats (recall, precision, F1) for each size bin in DEL, DUP, INV
    bin_order = ["0-50bp", "50-100bp", "100-300bp", "300-500bp",
                 "500-1000bp", "1-2kb", "2-5kb", "5-10kb", "10kb+"]
    metric_df = get_summary_stats(overlap_df_non_ins, args.tools, bin_order)
    metric_df.to_csv('%s/DEL_DUP_INV_metrics_by_svtype_and_size.tsv' % args.output_path, sep='\t')
    # collapse size bins to calculate overall metrics for each SVTYPE
    metric_df_mean = metric_df.groupby(['SVTYPE', 'tool'], as_index=False).mean()
    metric_df_mean.to_csv('%s/DEL_DUP_INV_metrics_by_svtype.tsv' % args.output_path, sep='\t')

    # plot recall, precision, F1 by SVTYPE and size
    for metric in ['recall', 'precision', 'F1']:
        plot_metrics(metric_df, metric, args.output_path)

    # plot overall recall, precision, F1
    plot_overall_metrics(metric_df, args.output_path)

    # output TP, FN, FP
    output_TP_FP_FN(overlap_df_non_ins, args.output_path)

    # analyse insertions
    ins_overlaps = glob.glob('%s/*overlap_INS.tsv' % args.overlap_dir)

    # make dataframe aggregating SVs from all nine gnomAD samples
    ins_overlap_df = make_overlap_df(ins_overlaps, col_names, size_bins)

    # add aggregate prediction (0 or 1) across all tools for each SV and get metrics
    tools = args.tools[:-1]
    ins_overlap_df = get_prediction(ins_overlap_df, tools)
    ins_metric_df = get_summary_stats_ins(ins_overlap_df, args.tools)
    ins_metric_df.to_csv('%s/INS_metrics_by_svtype.tsv' % args.output_path, sep='\t')

    # plot overall recall, precision, F1 for each INS type
    plot_overall_metrics(ins_metric_df, args.output_path, ins=True)

    # output TP, FN, FP
    output_TP_FP_FN(ins_overlap_df, args.output_path, ins=True)

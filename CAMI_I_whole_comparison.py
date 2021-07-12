"""
This script is used to reproduce the plot of the CAMI I results(compared to other binners).
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def get_number_of_genomes_per_completeness(amber_path, return_pandas=False):
    genome_path = os.path.join(amber_path, 'genome')
    table = {}
    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            method_path = os.path.join(root, name)

            metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
            metric = metric.query('`Purity (bp)` >= 0.95')
            com_90_pur_95 = metric.eval('`Completeness (bp)` > 0.9').sum()
            com_80_pur_95 = metric.eval('`Completeness (bp)` > 0.8').sum() - (
                        com_90_pur_95)
            com_70_pur_95 = metric.eval('`Completeness (bp)` > 0.7').sum() - (
                        com_90_pur_95 + com_80_pur_95)
            com_60_pur_95 = metric.eval('`Completeness (bp)` > 0.6').sum() - (
                        com_90_pur_95 + com_80_pur_95 + com_70_pur_95)
            com_50_pur_95 = metric.eval('`Completeness (bp)` > 0.5').sum() - (
                        com_90_pur_95 + com_80_pur_95 + com_70_pur_95 + com_60_pur_95)
            table[method_path.split('/')[-1]]= [com_90_pur_95, com_80_pur_95, com_70_pur_95, com_60_pur_95, com_50_pur_95]
    if return_pandas:
        return pd.DataFrame(table, index=[90,80,70,60,50]).T
    return table

def plot_bar(amber_path, add_legend=True,y_label = None,title = None, output = None):
    data = get_number_of_genomes_per_completeness(amber_path, return_pandas=True)
    subset = data.loc[[
            'COCACOLA',
            'SolidBin-naive',
            'SolidBin-CL',
            'SolidBin-SFS-CL',
            'SolidBin-coalign',
            'Maxbin2',
            'Metabat2_200',
            'Vamb',
            'S3N2Bin_200']]
    subset.rename(index={
        'S3N2Bin_200': 'SemiBin',
        'Metabat2_200': 'Metabat2',
        }, inplace=True)
    high_quality_list = subset[90].sort_values().values

    print('Improvement of best binner over second best: {:.2%}'.format(
            (high_quality_list[-1] - high_quality_list[-2])/high_quality_list[-2]))
    ax = subset.plot(kind="barh", stacked=True, legend = False, color = ['#084594','#4292c6','#9ecae1','#deebf7'])

    if add_legend:
        ax.legend(['>90', '>80','>70','>60'],
                  loc='lower right', fontsize=10, title='completeness')
    ax.set_xticks(ticks=y_label)
    ax.set_xticklabels(labels=y_label,fontsize=15,color = 'black')
    ax.set_yticklabels(labels=subset.index, fontsize=15,color = 'black')
    ax.set_xlabel('Bins(< 5% contamination)', fontsize=15,color = 'black')
    ax.set_title(title, fontsize=15, alpha=1.0)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    plt.show()

def plot_f1_score(amber_path, y_label = None,title = None,size = 5, output = None):
    genome_path = os.path.join(amber_path, 'genome')
    methold_path_list = []
    method_list = {}
    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            methold_path_list.append(os.path.join(root, name))


    for method_path in methold_path_list:
        if method_path.split('/')[-1] == 'Gold standard':
            continue
        metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
        metric = metric[(metric['Completeness (bp)'] >= 0.5) & (metric['Purity (bp)'] >= 0.5)]
        com = metric[['Completeness (bp)']]
        pur = metric[['Purity (bp)']]
        com = com.fillna(0).values
        pur = pur.fillna(0).values
        f1 = 2 * (com * pur) / (com + pur)
        f1 = pd.DataFrame(f1)
        f1.columns = ['F1']
        method_list[method_path.split('/')[-1]] = f1

    method = ['COCACOLA', 'SolidBin-naive', 'SolidBin-CL', 'SolidBin-SFS-CL', 'SolidBin-coalign', 'Maxbin2',
              'Metabat2_200', 'Vamb', 'S3N2Bin_200']
    subset = method_list['COCACOLA']
    for temp in method:
        if temp == 'COCACOLA':
            continue
        subset = pd.concat([subset, method_list[temp]], axis=1)
    fig, ax = plt.subplots(nrows=1, ncols=1)
    subset.columns = [
        ['COCACOLA', 'SolidBin-naive', 'SolidBin-CL', 'SolidBin-SFS-CL', 'SolidBin-coalign', 'Maxbin2', 'Metabat2_200',
         'Vamb', 'SemiBin_200']]
    column_new = ['SemiBin_200', 'Vamb', 'Metabat2_200', 'Maxbin2', 'SolidBin-coalign', 'SolidBin-SFS-CL',
                  'SolidBin-CL', 'SolidBin-naive', 'COCACOLA']
    subset = subset[column_new]
    sns.stripplot(data=subset, size=size,
                  order=['SemiBin_200', 'Vamb', 'Metabat2_200', 'Maxbin2', 'SolidBin-coalign', 'SolidBin-SFS-CL',
                         'SolidBin-CL', 'SolidBin-naive', 'COCACOLA'], orient='h',palette= ['#a6cee3','#66a61e','#e6ab02','#e7298a','#1b9e77','#d95f02','#7570b3','#a6761d','#666666'])

    sns.boxplot(data=subset, orient="h", color='white', width=.5, fliersize=0)
    ax.set_xticks(ticks=y_label)
    ax.set_xticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_yticklabels(
        labels=['SemiBin', 'Vamb', 'Metabat2', 'Maxbin2', 'SolidBin-coalign', 'SolidBin-SFS-CL', 'SolidBin-CL',
                'SolidBin-naive', 'COCACOLA'], fontsize=15, color='black')
    ax.set_xlabel('F1-score', fontsize=15, color='black')
    ax.set_title('{}'.format(title), fontsize=15, alpha=1.0)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_SemiBin_Metabat(amber_path,if_legend=True,y_label=None, output = None):
    method_list = get_number_of_genomes_per_completeness(amber_path)

    fig, axes = plt.subplots(nrows=1, ncols=3)

    ax_position = 0

    subset = pd.DataFrame(np.array([method_list['Metabat2_200'][0:4], method_list['S3N2Bin_200'][0:4]]),
                          index=['Metabat2', 'SemiBin'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])

    if if_legend:
        ax.legend(['>90', '>80', '>70', '>60'],
                  loc='upper left', fontsize=10)
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['Metabat2', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_ylabel('Bins', fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=200'), fontsize=13, alpha=1.0)
    ax_position += 1

    subset = pd.DataFrame(np.array([method_list['Metabat2_500'][0:4], method_list['S3N2Bin_500'][0:4]]),
                          index=['Metabat2', 'SemiBin'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['Metabat2', 'SemiBin'], rotation=50, fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=500'), fontsize=13)
    ax_position += 1

    subset = pd.DataFrame(np.array([method_list['Metabat2_1000'][0:4], method_list['S3N2Bin_1000'][0:4]]),
                          index=['Metabat2', 'SemiBin'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])

    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['Metabat2', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=1000'), fontsize=13, alpha=1.0)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_CAT_mmseqs(amber_path,if_legend=True,y_label=None,output = None):
    method_list = get_number_of_genomes_per_completeness(amber_path)

    fig, axes = plt.subplots(nrows=1, ncols=3)

    ax_position = 0

    subset = pd.DataFrame(
        np.array([method_list['S3N2Bin_CAT_200'][0:4], method_list['S3N2Bin_mmseqs_200'][0:4]]),
        index=['SemiBin_CAT', 'SemiBin_mmseqs'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])

    if if_legend:
        ax.legend(['>90', '>80', '>70', '>60'],
                  loc='lower left', fontsize=10, title='completeness')
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['CAT', 'mmseqs'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_ylabel('Bins', fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=200'), fontsize=13, alpha=1.0)
    ax_position += 1

    subset = pd.DataFrame(
        np.array([method_list['S3N2Bin_CAT_500'][0:4], method_list['S3N2Bin_mmseqs_500'][0:4]]),
        index=['SemiBin_CAT', 'SemiBin_mmseqs'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['CAT', 'mmseqs'], rotation=50, fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=500'), fontsize=13)
    ax_position += 1

    subset = pd.DataFrame(
        np.array([method_list['S3N2Bin_CAT_1000'][0:4], method_list['S3N2Bin_mmseqs_1000'][0:4]]),
        index=['SemiBin_CAT', 'SemiBin_mmseqs'], columns=[90, 80, 70, 60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position], legend=False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])

    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_xticklabels(labels=['CAT', 'mmseqs'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_title('{}'.format('max_edges=1000'), fontsize=13, alpha=1.0)

    # plt.suptitle(title,y = 1)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    plt.show()

def plot_bar_semi_no_semi(amber_path,if_legend=True,y_label = None, output = None):
    method_list = get_number_of_genomes_per_completeness(amber_path)
    fig, axes = plt.subplots(nrows=1, ncols=3)

    ax_position = 0

    subset= pd.DataFrame(np.array([method_list['NoSemi_200'][0:4],method_list['S3N2Bin_200'][0:4]]),index=['NoSemi','SemiBin'],columns=[90,80,70,60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position],legend = False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])
    print('Semi improvement:', (method_list['S3N2Bin_200'][0:4][0] - method_list['NoSemi_200'][0])/ method_list['NoSemi_200'][0])
    if if_legend:
        ax.legend(['>90', '>80','>70','>60'],
                  loc='upper left', fontsize=10, title="completeness",)
    ax.set_yticklabels(labels=y_label,fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['NoSemi','SemiBin'], rotation=50,
                       minor=False, fontsize=15,color = 'black')
    ax.set_ylabel('Bins', fontsize=15,color = 'black')
    ax.set_title('{}'.format('max_edges=200'), fontsize=15, alpha=1.0,color = 'black')
    ax_position += 1

    subset= pd.DataFrame(np.array([method_list['NoSemi_500'][0:4],method_list['S3N2Bin_500'][0:4]]),index=['NoSemi','SemiBin'],columns=[90,80,70,60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position],legend = False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])
    ax.set_yticklabels(labels=y_label,fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['NoSemi','SemiBin'], rotation=50, fontsize=15,color = 'black')
    ax.set_title('{}'.format('max_edges=500'), fontsize=15,color = 'black')
    ax_position += 1

    subset= pd.DataFrame(np.array([method_list['NoSemi_1000'][0:4],method_list['S3N2Bin_1000'][0:4]]),index=['NoSemi','SemiBin'],columns=[90,80,70,60])
    ax = subset.plot(kind="bar", stacked=True,
                     ax=axes[ax_position],legend = False, color = ['#084594','#4292c6','#9ecae1','#c6dbef'])

    ax.set_yticklabels(labels=y_label,fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['NoSemi','SemiBin'], rotation=50,
                       minor=False, fontsize=15,color = 'black')
    ax.set_title('{}'.format('max_edges=1000'), fontsize=15, alpha=1.0,color = 'black')

    #plt.suptitle(title,y = 1)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    plt.show()

def plot_f1_boxplot_semi_to_nosemi(amber_path, y_label = None,size = 5, output=None):
    genome_path = os.path.join(amber_path, 'genome')
    methold_path_list = []
    method_list = {}
    for root, dirs, files in os.walk(genome_path, topdown=False):
        for name in dirs:
            methold_path_list.append(os.path.join(root, name))

    for method_path in methold_path_list:
        if method_path.split('/')[-1] == 'Gold standard':
            continue
        metric = pd.read_csv(os.path.join(method_path, 'metrics_per_bin.tsv'), sep='\t')
        metric = metric[(metric['Completeness (bp)'] >= 0.5) & (metric['Purity (bp)'] >= 0.5)]
        com = metric[['Completeness (bp)']]
        pur = metric[['Purity (bp)']]
        com = com.fillna(0).values
        pur = pur.fillna(0).values
        f1 = 2 * (com * pur) / (com + pur)
        f1 = pd.DataFrame(f1)
        f1.columns = ['F1']
        method_list[method_path.split('/')[-1]] = f1

    fig, axes = plt.subplots(nrows=1, ncols=3)
    #
    # ax_position = 0
    #
    subset = pd.concat([method_list['NoSemi_200'],method_list['S3N2Bin_200']],axis = 1)
    subset.columns = [['NoSemi_200','SemiBin_200']]
    sns.stripplot(data=subset, ax=axes[0], size=size,palette=['#b2df8a', '#a6cee3'])
    sns.boxplot(data=subset,ax=axes[0], color='white', width=.5,fliersize = 0)
    axes[0].set_yticks(ticks=y_label)
    axes[0].set_yticklabels(labels=y_label, fontsize=12, color='black')
    axes[0].set_xticklabels(labels=['NoSemi', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15,color = 'black')
    axes[0].set_ylabel('F1-score', fontsize=15,color = 'black')
    axes[0].set_title('{}'.format('max_edges=200'), fontsize=15, alpha=1.0,color = 'black')

    subset = pd.concat([method_list['NoSemi_500'],method_list['S3N2Bin_500']],axis = 1)
    subset.columns = [['NoSemi_500','SemiBin_500']]
    sns.stripplot(data=subset, ax=axes[1], size=size,palette=['#b2df8a', '#a6cee3'])
    sns.boxplot(data=subset,ax=axes[1], color='white', width=.5,fliersize = 0)
    axes[1].set_yticks(ticks=y_label)
    axes[1].set_yticklabels(labels=y_label, fontsize=12, color='black')
    axes[1].set_xticklabels(labels=['NoSemi', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15,color = 'black')
    axes[1].set_title('{}'.format('max_edges=500'), fontsize=15, alpha=1.0,color = 'black')

    subset = pd.concat([method_list['NoSemi_1000'],method_list['S3N2Bin_1000']],axis = 1)
    subset.columns = [['NoSemi_1000','SemiBin_1000']]
    sns.stripplot(data=subset, ax=axes[2], size=size,palette=['#b2df8a', '#a6cee3'])
    sns.boxplot(data=subset,ax=axes[2], color='white', width=.5,fliersize = 0)
    axes[2].set_yticks(ticks=y_label)
    axes[2].set_yticklabels(labels=y_label, fontsize=12, color='black')
    axes[2].set_xticklabels(labels=['NoSemi', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15,color = 'black')
    axes[2].set_title('{}'.format('max_edges=1000'), fontsize=15, alpha=1.0,color = 'black')


    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    plt.show()

def plot_bar_generalization(amber_path,No_semi_amber_path, title,if_legend=True,y_label = None, output=None):
    method_list = get_number_of_genomes_per_completeness(amber_path)
    no_semi_method_list = get_number_of_genomes_per_completeness(No_semi_amber_path)
    method = ['S3N2Bin_m', 'S3N2Bin_c', 'S3N2Bin_mc', 'S3N2Bin']
    data = []
    data.append(no_semi_method_list['NoSemi_200'][0])
    for temp in method:
        data.append(method_list[temp][0])
    data = np.array(data)

    subset = pd.DataFrame(data.reshape(1,data.shape[0]), columns =['NoSemi','SemiBin_m', 'SemiBin_c', 'SemiBin_mc', 'SemiBin'])

    ax = subset.plot(kind="bar",legend=False, color = ['#fdbf6f','#fb9a99', '#b2df8a', '#a6cee3', '#1f78b4'])

    if if_legend:
        ax.legend(loc = 'lower left', fontsize=10)
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_title('{}'.format(title), fontsize=15, alpha=1.0)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    plt.show()

if __name__ == '__main__':
    base_path = 'Results/Simulated/CAMI_I/whole_comparasion/'

    ### whole comparison bar plot
    # amber_path_low = base_path + 'low_whole'
    # amber_path_medium = base_path + 'medium_whole'
    # amber_path_high = base_path + 'high_whole'
    # plot_bar(amber_path_low, y_label=[0, 5, 10, 15, 20, 25, 30], title='Low-complexity(all)', output='CAMI_I_low_whole.pdf')
    # plot_bar(amber_path_medium, y_label=[0, 20, 40, 60, 80, 100], title='Medium-complexity(all)', if_legend=False, output='CAMI_I_medium_whole.pdf')
    # plot_bar(amber_path_high, y_label=[0, 100, 200, 300, 400], title='High-complexity(all)', if_legend=False, output='CAMI_I_high_whole.pdf')

    ## whole comparison F1 box plot
    # plot_f1_score(amber_path_low,y_label=[0.6,0.7,0.8,0.9,1.0],title='Low-complexity(all)', output='CAMI_I_low_whole_F1.pdf')
    # plot_f1_score(amber_path_medium,title='Medium-complexity(all)',y_label=[0.5,0.6,0.7,0.8,0.9,1.0],size = 3, output='CAMI_I_medium_whole_F1.pdf')
    # plot_f1_score(amber_path_high,title='High-complexity(all)',y_label=[0.5,0.6,0.7,0.8,0.9,1.0],size=2, output='CAMI_I_high_whole_F1.pdf')
    #
    # ### common comparison bar plot
    # amber_path_low_common = base_path + 'low_common'
    # amber_path_medium_common = base_path + 'medium_common'
    # amber_path_high_common = base_path + 'high_common'
    # plot_bar(amber_path_low_common, y_label=[0, 2, 4, 6, 8, 10], title='Low-complexity(common strain)', output='CAMI_I_low_common.pdf')
    # plot_bar(amber_path_medium_common, y_label=[0, 10, 20, 30, 40, 50], title='Medium-complexity(common strain)',
    #        add_legend=False,output='CAMI_I_medium_common.pdf')
    # plot_bar(amber_path_high_common, y_label=[0, 30, 60, 90, 120, 150], title='High-complexity(common strain)',
    #        add_legend=False,output='CAMI_I_high_common.pdf')
    
    # ### unique comparison bar plot
    # amber_path_low_unique = base_path + 'low_unique'
    # amber_path_medium_unique = base_path + 'medium_unique'
    # amber_path_high_unique = base_path + 'high_unique'
    #
    # plot_bar(amber_path_low_unique, y_label=[0, 5, 10, 15], title='Low-complexity(unique strain)',  add_legend=False, output='CAMI_I_low_unique.pdf')
    # plot_bar(amber_path_medium_unique, y_label=[0, 15, 30, 45, 60, 75], title='Medium-complexity(unique strain)',
    #           add_legend=False, output='CAMI_I_medium_unique.pdf')
    # plot_bar(amber_path_high_unique, y_label=[0, 100, 200, 300, 400], title='High-complexity(unique strain)',
    #           add_legend=False, output='CAMI_I_high_unique.pdf')
    #
    # ### compare to Metabat2 with different parameters
    # plot_SemiBin_Metabat(amber_path_low,if_legend=True,y_label=[0,5,10,15,20,25,30], output='CAMI_I_low_SemiBin_Metabat2.pdf')
    # plot_SemiBin_Metabat(amber_path_medium,y_label=[0,20,40,60,80,100], output='CAMI_I_medium_SemiBin_Metabat2.pdf')
    # plot_SemiBin_Metabat(amber_path_high,y_label=[0,100,200,300,400,500], output='CAMI_I_high_SemiBin_Metabat2.pdf')
    #
    # base_path = 'Results/Simulated/CAMI_I/mmseqs_CAT/'
    # amber_low_mmseqs_CAT = base_path + 'amber_low'
    # amber_medium_mmseqs_CAT = base_path + 'amber_medium'
    # amber_high_mmseqs_CAT = base_path + 'amber_high'
    #
    # ### compare results with CAT or mmseqs
    # plot_CAT_mmseqs(amber_low_mmseqs_CAT, if_legend=True,y_label=[0,5,10,15,20,25,30],output='CAMI_I_CAT_mmseqs_low.pdf')
    # plot_CAT_mmseqs(amber_medium_mmseqs_CAT,y_label=[0,20,40,60,80,100],output='CAMI_I_CAT_mmseqs_medium.pdf')
    # plot_CAT_mmseqs(amber_high_mmseqs_CAT,y_label=[0,100,200,300,400,500],output='CAMI_I_CAT_mmseqs_high.pdf')
    #
    # ## bar plot with Semi VS no semi
    amber_base = 'Results/Simulated/CAMI_I/No_semi_mmseqs/'
    amber_path_low_nosemi = amber_base + 'low_no_semi'
    amber_path_medium_nosemi = amber_base + 'medium_no_semi'
    amber_path_high_nosemi = amber_base + 'high_no_semi'

    plot_bar_semi_no_semi(amber_path_low_nosemi, y_label=[0,5,10,15,20,25,30],output='CAMI_I_low_semi_nosemi.pdf')
    plot_bar_semi_no_semi(amber_path_medium_nosemi,  y_label=[0,20,40,60,80,100],output='CAMI_I_medium_semi_nosemi.pdf')
    plot_bar_semi_no_semi(amber_path_high_nosemi, y_label=[0,100,200,300,400,500],output='CAMI_I_high_semi_nosemi.pdf')
    #
    # ### F1 box plot with Semi vs no semi
    # plot_f1_boxplot_semi_to_nosemi(amber_path_low_nosemi, y_label=[0.6, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00],
    #                                size=5,output='CAMI_I_low_semi_nosemi_F1.pdf' )
    # plot_f1_boxplot_semi_to_nosemi(amber_path_medium_nosemi, y_label=[0.6, 0.7, 0.8, 0.9, 1.0], size=3,output='CAMI_I_medium_semi_nosemi_F1.pdf')
    # plot_f1_boxplot_semi_to_nosemi(amber_path_high_nosemi, y_label=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0], size=2,output='CAMI_I_high_semi_nosemi_F1.pdf')
    #
    # ### comparison of generalization(S3N2Bin_c, S3N2Bin_m, S3N2Bin_mc)
    # base_path = 'Results/Simulated/CAMI_I/generalization_comparision/'
    # amber_low_path = base_path + 'amber_low_generalization'
    # amber_medium_path = base_path + 'amber_medium_generalization'
    # amber_high_path = base_path + 'amber_high_generalization'
    #
    # plot_bar_generalization(amber_low_path,amber_path_low_nosemi, 'Low complexity',y_label=[0,5,10,15,20,25], output='CAMI_I_generalization_low.pdf')
    # plot_bar_generalization(amber_medium_path, amber_path_medium_nosemi,'Medium complexity', y_label=[0,20,40,60,80],if_legend=False, output='CAMI_I_generalization_medium.pdf')
    # plot_bar_generalization(amber_high_path, amber_path_high_nosemi, 'High complexity', y_label=[0,100,200,300,400],if_legend=False, output='CAMI_I_generalization_high.pdf')
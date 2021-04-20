"""
This script is used to reproduct the plot of the real datasets
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3_unweighted


contamination = 0.05
dog_list = ['SAMN06172456', 'SAMN06172425', 'SAMN06172487', 'SAMN06172450', 'SAMN06172459', 'SAMN06172479', 'SAMN06172435', 'SAMN06172414', 'SAMN06172409', 'SAMEA103957796', 'SAMN06172442', 'SAMN06172500', 'SAMN06172437', 'SAMN06172413', 'SAMN06172514', 'SAMN06172403', 'SAMN06172471', 'SAMN06172490', 'SAMN06172448', 'SAMN06172504', 'SAMN06172457', 'SAMN06172441', 'SAMN06172422', 'SAMN06172408', 'SAMN06172429', 'SAMN06172420', 'SAMN06172503', 'SAMN06172410', 'SAMN06172458', 'SAMN06172493', 'SAMEA103957794', 'SAMN06172402', 'SAMN06172515', 'SAMN06172462', 'SAMN06172421', 'SAMN06172411', 'SAMN06172511', 'SAMN06172516', 'SAMN06172465', 'SAMN06172419', 'SAMN06172517', 'SAMN06172510', 'SAMN06172418', 'SAMN06172424', 'SAMN06172427', 'SAMN06172453', 'SAMN06172491', 'SAMN06172496', 'SAMN06172513', 'SAMN06172461', 'SAMN06172449', 'SAMN06172426', 'SAMN06172452', 'SAMN06172522', 'SAMN06172400', 'SAMN06172405', 'SAMN06172521', 'SAMN06172407', 'SAMN06172455', 'SAMN06172446', 'SAMN06172467', 'SAMN06172499', 'SAMN06172474', 'SAMN06172412', 'SAMN06172468', 'SAMN06172478', 'SAMN06172423', 'SAMN06172447', 'SAMN06172415', 'SAMN06172523', 'SAMN06172417', 'SAMN06172497', 'SAMN06172498', 'SAMN06172489', 'SAMN06172436', 'SAMN06172432', 'SAMN06172406', 'SAMN06172488', 'SAMN06172502', 'SAMN06172401', 'SAMN06172434', 'SAMN06172416', 'SAMN06172445', 'SAMN06172431', 'SAMN06172438', 'SAMN06172473', 'SAMN06172486', 'SAMN06172472', 'SAMN06172428', 'SAMEA103957793', 'SAMEA103957795', 'SAMN06172443', 'SAMN06172475', 'SAMN06172520', 'SAMN06172495', 'SAMN06172440', 'SAMN06172430', 'SAMN06172481', 'SAMN06172524', 'SAMN06172519', 'SAMN06172454', 'SAMN06172404', 'SAMN06172460', 'SAMN06172433', 'SAMN06172469', 'SAMN06172451', 'SAMN06172476', 'SAMN06172492', 'SAMN06172484', 'SAMN06172509', 'SAMN06172506', 'SAMN06172518', 'SAMN06172477', 'SAMN06172470', 'SAMN06172482', 'SAMN06172512', 'SAMN06172494', 'SAMN06172485', 'SAMN06172508', 'SAMN06172466', 'SAMN06172507', 'SAMN06172444', 'SAMN06172505', 'SAMN06172464', 'SAMN06172439', 'SAMN06172501', 'SAMN06172483', 'SAMN06172463', 'SAMN06172480']

human_list = ['CCMD41521570ST', 'CCMD75147712ST', 'CCMD18579000ST', 'CCMD53508245ST', 'CCMD19168690ST', 'CCMD52117727ST', 'CCMD42956136ST', 'CCMD79349503ST', 'CCMD89306485ST', 'CCMD76409700ST', 'CCMD31134579ST', 'CCMD71242853ST', 'CCMD89107682ST', 'CCMD76222476ST', 'CCMD10032470ST', 'CCMD17410933ST', 'CCMD38158721ST', 'CCMD35081859ST', 'CCMD54057834ST', 'CCMD28738636ST', 'CCMD98702133ST', 'CCMD30626189ST', 'CCMD32965613ST', 'CCMD53522274ST', 'CCMD37575804ST', 'CCMD68973846ST', 'CCMD25475945ST', 'CCMD65406197ST', 'CCMD21703880ST', 'CCMD50300306ST', 'CCMD51228890ST', 'CCMD59540613ST', 'CCMD49942357ST', 'CCMD95431029ST', 'CCMD41202658ST', 'CCMD15562448ST', 'CCMD21593359ST', 'CCMD92404903ST', 'CCMD50538120ST', 'CCMD49461418ST', 'CCMD72690923ST', 'CCMD85481373ST', 'CCMD39882286ST', 'CCMD18829815ST', 'CCMD51154251ST', 'CCMD85661207ST', 'CCMD71915439ST', 'CCMD39157124ST', 'CCMD22852639ST', 'CCMD35801800ST', 'CCMD27463710ST', 'CCMD59583015ST', 'CCMD89967135ST', 'CCMD52145360ST', 'CCMD95676152ST', 'CCMD45004878ST', 'CCMD67373733ST', 'CCMD99929634ST', 'CCMD89643949ST', 'CCMD26625622ST', 'CCMD23541216ST', 'CCMD31009081ST', 'CCMD99440714ST', 'CCMD66848156ST', 'CCMD65222621ST', 'CCMD98531134ST', 'CCMD45812507ST', 'CCMD46727384ST', 'CCMD73128545ST', 'CCMD30627121ST', 'CCMD50529145ST', 'CCMD98198513ST', 'CCMD93755960ST', 'CCMD35633353ST', 'CCMD56948710ST', 'CCMD27867141ST', 'CCMD32288175ST', 'CCMD29706695ST', 'CCMD72666896ST', 'CCMD10191450ST', 'CCMD49025643ST', 'CCMD74592084ST']

tara_list = ['TARA_041_SRF_0.1-0.22', 'TARA_038_SRF_0.22-1.6', 'TARA_076_SRF_0.22-3', 'TARA_023_SRF_0.22-1.6', 'TARA_042_SRF_0.22-1.6', 'TARA_124_SRF_0.22-3', 'TARA_124_SRF_0.22-0.45', 'TARA_066_SRF_lt-0.22', 'TARA_057_SRF_0.22-3', 'TARA_124_SRF_0.45-0.8', 'TARA_004_SRF_0.22-1.6', 'TARA_018_SRF_0.22-1.6', 'TARA_070_SRF_0.22-0.45', 'TARA_034_SRF_lt-0.22', 'TARA_064_SRF_0.22-3', 'TARA_125_SRF_0.22-0.45', 'TARA_111_SRF_0.22-3', 'TARA_122_SRF_0.22-0.45', 'TARA_145_SRF_0.22-3', 'TARA_099_SRF_0.22-3', 'TARA_038_SRF_lt-0.22', 'TARA_082_SRF_0.22-3', 'TARA_041_SRF_lt-0.22', 'TARA_146_SRF_0.22-3', 'TARA_151_SRF_0.22-3', 'TARA_123_SRF_0.22-3', 'TARA_110_SRF_0.22-3', 'TARA_150_SRF_0.22-3', 'TARA_072_SRF_lt-0.22', 'TARA_085_SRF_0.22-3', 'TARA_098_SRF_0.22-3', 'TARA_078_SRF_0.22-3', 'TARA_149_SRF_0.22-3', 'TARA_094_SRF_0.22-3', 'TARA_068_SRF_0.22-3', 'TARA_148_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_018_SRF_lt-0.22', 'TARA_138_SRF_0.22-3', 'TARA_093_SRF_0.22-3', 'TARA_041_SRF_0.22-1.6', 'TARA_122_SRF_0.22-3', 'TARA_078_SRF_0.45-0.8', 'TARA_070_SRF_lt-0.22', 'TARA_065_SRF_lt-0.22', 'TARA_122_SRF_0.1-0.22', 'TARA_036_SRF_0.22-1.6', 'TARA_031_SRF_0.22-1.6', 'TARA_142_SRF_0.22-3', 'TARA_124_SRF_0.1-0.22', 'TARA_036_SRF_lt-0.22', 'TARA_065_SRF_0.22-3', 'TARA_067_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_109_SRF_0.22-3', 'TARA_068_SRF_lt-0.22', 'TARA_109_SRF_lt-0.22', 'TARA_064_SRF_lt-0.22', 'TARA_048_SRF_0.22-1.6', 'TARA_034_SRF_0.22-1.6', 'TARA_070_SRF_0.45-0.8', 'TARA_025_SRF_lt-0.22', 'TARA_133_SRF_0.22-3', 'TARA_096_SRF_0.22-3', 'TARA_038_SRF_0.1-0.22', 'TARA_007_SRF_0.22-1.6', 'TARA_048_SRF_0.1-0.22', 'TARA_140_SRF_0.22-3', 'TARA_034_SRF_0.1-0.22', 'TARA_067_SRF_0.22-3', 'TARA_125_SRF_0.45-0.8', 'TARA_030_SRF_0.22-1.6', 'TARA_031_SRF_lt-0.22', 'TARA_032_SRF_0.22-1.6', 'TARA_070_SRF_0.22-3', 'TARA_132_SRF_0.22-3', 'TARA_076_SRF_lt-0.22', 'TARA_125_SRF_0.1-0.22', 'TARA_123_SRF_0.45-0.8', 'TARA_078_SRF_lt-0.22', 'TARA_068_SRF_0.45-0.8', 'TARA_068_SRF_0.22-0.45', 'TARA_067_SRF_0.22-0.45', 'TARA_100_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_137_SRF_0.22-3', 'TARA_076_SRF_0.22-0.45', 'TARA_125_SRF_0.22-3', 'TARA_078_SRF_0.22-0.45', 'TARA_076_SRF_0.45-0.8', 'TARA_084_SRF_0.22-3', 'TARA_032_SRF_lt-0.22', 'TARA_025_SRF_0.22-1.6', 'TARA_062_SRF_0.22-3', 'TARA_066_SRF_0.22-3', 'TARA_036_SRF_0.1-0.22', 'TARA_056_SRF_0.22-3', 'TARA_072_SRF_0.22-3', 'TARA_128_SRF_0.22-3', 'TARA_052_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_123_SRF_0.22-0.45', 'TARA_102_SRF_0.22-3', 'TARA_065_SRF_0.1-0.22', 'TARA_009_SRF_0.22-1.6', 'TARA_141_SRF_0.22-3', 'TARA_045_SRF_0.22-1.6', 'TARA_042_SRF_lt-0.22', 'TARA_152_SRF_0.22-3']

def get_result(dataset='dog', method='Maxbin2', binning_mode = 'single_sample'):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample
    """
    if dataset == 'dog':
        sample_list = dog_list
    if dataset == 'human':
        sample_list = human_list
    if dataset == 'tara':
        sample_list = tara_list
    result = {}
    if method == 'VAMB' and binning_mode == 'multi_sample':
        result = {'high quality': [], 'medium quality': [], 'low quality': []}
        binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/multi_sample/VAMB_multi.csv'.format(dataset),index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                binning_result['Contamination'].astype(float) < float(contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        result['high quality'].extend(high_quality.index.tolist())

        medium_quality = binning_result[(binning_result['Completeness'].astype(float) <= float(90)) & (binning_result['Completeness'].astype(float) > float(50)) &(binning_result['Contamination'].astype(float) < float(contamination * 100))]
        result['medium quality'].extend(medium_quality.index.tolist())

        low_quality = [temp for temp in binning_result.index.tolist() if temp not in result['high quality'] and temp not in result['medium quality']]
        result['low quality'].extend(low_quality)

        return result

    else:
        for sample in sample_list:
            result[sample] = {'high quality':[], 'medium quality':[], 'low quality':[]}
            binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/{1}/{2}/{3}/result.csv'.format(dataset,binning_mode, sample, method),
                                         index_col=0)
            high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                                                      binning_result['pass.GUNC'] == True)]
            result[sample]['high quality'].extend(high_quality.index.tolist())

            medium_quality = binning_result[(binning_result['Completeness'].astype(float) <= float(90)) & (
                        binning_result['Completeness'].astype(float) > float(50)) & (
                                                        binning_result['Contamination'].astype(float) < float(
                                                    contamination * 100))]
            result[sample]['medium quality'].extend(medium_quality.index.tolist())

            low_quality = [temp for temp in binning_result.index.tolist() if
                           temp not in result[sample]['high quality'] and temp not in result[sample]['medium quality']]
            result[sample]['low quality'].extend(low_quality)
        return result

def get_num_high_quality(dataset = 'dog', method = 'Maxbin2', binning_mode = 'single_sample'):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample
    """
    result = get_result(dataset, method, binning_mode)
    if method == 'VAMB' and binning_mode == 'multi_sample':
        num_hq = len(result['high quality'])
        return num_hq
    else:
        num_hq = 0
        for sample in result:
            num_hq += len(result[sample]['high quality'])
        return num_hq

def plot_high_quality_comparison():
    num_dog_maxbin2_single = get_num_high_quality()
    num_dog_vamb_single = get_num_high_quality(method='VAMB')
    num_dog_metabat2_single = get_num_high_quality(method='Metabat2')
    num_dog_s3n2bin_single = get_num_high_quality(method='S3N2Bin')

    num_human_maxbin2_single = get_num_high_quality(dataset='human')
    num_human_vamb_single = get_num_high_quality(dataset='human', method='VAMB')
    num_human_metabat2_single = get_num_high_quality(dataset='human', method='Metabat2')
    num_human_s3n2bin_single = get_num_high_quality(dataset='human', method='S3N2Bin')

    num_tara_maxbin2_single = get_num_high_quality(dataset='tara')
    num_tara_vamb_single = get_num_high_quality(dataset='tara', method='VAMB')
    num_tara_metabat2_single = get_num_high_quality(dataset='tara', method='Metabat2')
    num_tara_s3n2bin_single = get_num_high_quality(dataset='tara', method='S3N2Bin')

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_s3n2bin_single]]),columns = ['Maxbin2','VAMB','Metabat2','$S^3N^2Bin$'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize=(3,4))
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.legend(loc = 'upper left',fontsize = 8)
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_s3n2bin_single]]),columns = ['Maxbin2','VAMB','Metabat2','$S^3N^2Bin$'], index=['Human gut'])
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False)
    ax.set_yticks(ticks=[0,300,600,900,1200,1500])
    ax.set_yticklabels(labels=[0,300,600,900,1200,1500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_s3n2bin_single]]),columns = ['Maxbin2','VAMB','Metabat2','$S^3N^2Bin$'], index=['Tara'])
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False)
    ax.set_yticks(ticks=[0,100,200,300,400])
    ax.set_yticklabels(labels=[0,100,200,300,400],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Tara'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    plt.show()


    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample')
    num_dog_s3n2bin_multi = get_num_high_quality(method='S3N2Bin',binning_mode='multi_sample')

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample')
    num_human_s3n2bin_multi = get_num_high_quality(dataset='human', method='S3N2Bin',binning_mode='multi_sample')

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample')
    num_tara_s3n2bin_multi = get_num_high_quality(dataset='tara', method='S3N2Bin', binning_mode='multi_sample')

    subset = pd.DataFrame(np.array([[num_dog_vamb_mulit,num_dog_s3n2bin_multi]]),columns = ['VAMB','$S^3N^2Bin$'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize = (2,4))
    ax.set_yticks(ticks=[0,800,1600,2400,3200])
    ax.set_yticklabels(labels=[0,800,1600,2400,3200],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.legend(loc='lower left', fontsize=6)
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_vamb_multi, num_human_s3n2bin_multi]]), columns=['VAMB', '$S^3N^2Bin$'],
                          index=['Human gut'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False)
    ax.set_yticks(ticks=[0, 300, 600, 900, 1200, 1500])
    ax.set_yticklabels(labels=[0, 300, 600, 900, 1200, 1500], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15, color='black', rotation=360)
    ax.set_ylabel('High quality genomes', fontsize=15, color='black')
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_vamb_multi,num_tara_s3n2bin_multi]]), columns=['VAMB', '$S^3N^2Bin$'],
                          index=['Tara'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False)
    ax.set_yticks(ticks=[0, 120, 240, 360, 480])
    ax.set_yticklabels(labels=[0, 120, 240, 360, 480], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Tara'], fontsize=15, color='black', rotation=360)
    ax.set_ylabel('High quality genomes', fontsize=15, color='black')
    plt.show()

def get_taxi_list(bac_path,  arr_path = None):
    bac = pd.read_csv(bac_path,
                      '\t').values

    species_list = []
    genus_list = []
    family_list = []

    for temp in bac:
        split = temp[1].split(';')

        species = split[-1]
        genus = split[-2]
        family = split[-3]

        if species != 's__':
            species_list.append(species)
        if genus != 'g__':
            genus_list.append(genus)
        if family != 'f__':
            family_list.append(family)

    if arr_path is not None:
        arr = pd.read_csv(arr_path,
                          '\t').values
        for temp in arr:
            split = temp[1].split(';')

            species = split[-1]
            genus = split[-2]
            family = split[-3]

            if species != 's__':
                species_list.append(species)
            if genus != 'g__':
                genus_list.append(genus)
            if family != 'f__':
                family_list.append(family)
    return set(family_list), set(genus_list), set(species_list)

def plot_multi_venn_comparison():
    base_path = 'Results/Real/gtdbtk_annotations/'
    ###Dog
    dog_VAMB_family, dog_VAMB_genus, dog_VAMB_species = get_taxi_list(base_path + 'multi_sample/dog/VAMB/gtdbtk.bac120.summary.tsv')
    dog_S3N2Bin_family, dog_S3N2Bin_genus, dog_S3N2Bin_species = get_taxi_list(base_path + 'multi_sample/dog/S3N2Bin/gtdbtk.bac120.summary.tsv')
    dog_S3N2Bin_family_single, dog_S3N2Bin_genus_single, dog_S3N2Bin_species_single = get_taxi_list(base_path + 'single_sample/dog/S3N2Bin/gtdbtk.bac120.summary.tsv')

    ### human
    human_VAMB_family, human_VAMB_genus, human_VAMB_species = get_taxi_list(base_path + 'multi_sample/human/VAMB/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/VAMB/gtdbtk.ar122.summary.tsv')
    human_S3N2Bin_family, human_S3N2Bin_genus, human_S3N2Bin_species = get_taxi_list(base_path + 'multi_sample/human/S3N2Bin/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/S3N2Bin/gtdbtk.ar122.summary.tsv')
    human_S3N2Bin_family_single, human_S3N2Bin_genus_single, human_S3N2Bin_species_single = get_taxi_list(base_path + 'single_sample/human/S3N2Bin/gtdbtk.bac120.summary.tsv',base_path + 'single_sample/human/S3N2Bin/gtdbtk.ar122.summary.tsv')

    ### tara
    tara_VAMB_family, tara_VAMB_genus, tara_VAMB_species = get_taxi_list(base_path + 'multi_sample/tara/VAMB/gtdbtk.bac120.summary.tsv')
    tara_S3N2Bin_family, tara_S3N2Bin_genus, tara_S3N2Bin_species = get_taxi_list(base_path + 'multi_sample/tara/S3N2Bin/gtdbtk.bac120.summary.tsv')
    tara_S3N2Bin_family_single, tara_S3N2Bin_genus_single, tara_S3N2Bin_species_single = get_taxi_list(base_path + 'single_sample/tara/S3N2Bin/gtdbtk.bac120.summary.tsv')

    out = venn3_unweighted([set(dog_VAMB_family), set(dog_S3N2Bin_family), set(dog_S3N2Bin_family_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Dog gut)", fontsize=20, alpha=1.0, color='black')
    #plt.savefig('multi_sample_dog_family.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_genus), set(dog_S3N2Bin_genus), set(dog_S3N2Bin_genus_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Dog gut)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_dog_genus.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_species), set(dog_S3N2Bin_species), set(dog_S3N2Bin_species_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Dog gut)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_dog_species.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_family), set(human_S3N2Bin_family), set(human_S3N2Bin_family_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Human gut)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_human_family.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_genus), set(human_S3N2Bin_genus), set(human_S3N2Bin_genus_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Human gut)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_human_genus.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_species), set(human_S3N2Bin_species), set(human_S3N2Bin_species_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Human gut)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_human_species.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_family), set(tara_S3N2Bin_family), set(tara_S3N2Bin_family_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Tara)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_tara_family.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_genus), set(tara_S3N2Bin_genus), set(tara_S3N2Bin_genus_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Tara)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_tara_genus.jpg')
    #plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_species), set(tara_S3N2Bin_species), set(tara_S3N2Bin_species_single)],
                           set_labels=('VAMB', 'S3N2Bin_multi', 'S3N2Bin_single'), normalize_to=5.0)
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Tara)", fontsize=20, alpha=1.0, color='black')

    #plt.savefig('multi_sample_tara_species.jpg')
    #plt.close()
    plt.show()

def plot_per_sample_comparison(dataset = 'dog'):
    """
    dataset: dog, human, tara
    """
    S3N2Bin_result = get_result(dataset, method='S3N2Bin')
    Metabat2_result =  get_result(dataset, method='Metabat2')

    result = {}

    for sample in S3N2Bin_result:
        if (len(S3N2Bin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality'])) not in result:
            result[(len(S3N2Bin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] = 1
        else:
            result[(len(S3N2Bin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] += 1



    if dataset != 'tara':
        data = np.zeros(shape=(len(result), 3))

        for i, temp in enumerate(result):
            data[i][0] = temp[0]
            data[i][1] = temp[1]
            data[i][2] = result[temp]

        data = pd.DataFrame(data).astype(int)

        data.columns = ['x', 'y', 'num']

        plt.scatter(data.x, data.y,
                    c=data.num, s=(data.num ** 2) * 60)

        plt.colorbar(shrink=0.5)
        plt.xlabel("$S^3N^2$Bin")
        plt.ylabel("Metabat2")
        if dataset == 'dog':
            plt.plot([0, 25], [0, 25], c='black')
            plt.title('Dog gut'.format(dataset), fontsize=15)
        else:
            plt.plot([0, 40], [0, 40], c='black')
            plt.title('Human gut'.format(dataset), fontsize=15)
        plt.show()
        #plt.savefig('/home1/pansj/plot/dog_compare.jpg', dpi=300)

    else:
        data = []
        for i, temp in enumerate(result):
            if temp[0] == 0 and temp[1] == 0:
                continue
            data.append([temp[0], temp[1], result[temp]])
        data = np.array(data)
        data = pd.DataFrame(data).astype(int)

        data.columns = ['x', 'y', 'num']

        plt.scatter(data.x, data.y,
                    c=data.num, s=(data.num) * 60)

        plt.colorbar()
        plt.xlabel("$S^3N^2$Bin")
        plt.ylabel("Metabat2")
        plt.plot([0, 25], [0, 25], c='black')
        plt.title('Tara', fontsize=15)
        plt.show()
        #plt.savefig('/home1/pansj/plot/tara_compare.jpg', dpi=300)

def get_overlap(dataset = 'dog'):
    if dataset == 'dog':
        sample_list = dog_list
    if dataset == 'human':
        sample_list = human_list
    if dataset == 'tara':
        sample_list = tara_list

    S3N2Bin_hq_list = []
    S3N2Bin_mq_list = []
    S3N2Bin_others_list = []

    Metabat2_hq_list = []
    Metabat2_mq_list = []
    Metabat2_others_list = []

    S3N2Bin_bin_dict = {}
    Metabat2_bin_dict = {}

    for sample in sample_list:
        S3N2Bin_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/S3N2Bin/result.csv'.format(dataset,sample),index_col=0)
        hq = S3N2Bin_binning_result[(S3N2Bin_binning_result['Completeness'].astype(float) > float(90)) & (
                S3N2Bin_binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                           S3N2Bin_binning_result['pass.GUNC'] == True)]
        all_values = S3N2Bin_binning_result.values
        all_bin = S3N2Bin_binning_result.index.tolist()
        for bin_index, checkm_index in zip(all_bin, all_values):
            S3N2Bin_bin_dict[sample + '_' + bin_index] = (checkm_index[10], checkm_index[11])

        mq = S3N2Bin_binning_result[(S3N2Bin_binning_result['Completeness'].astype(float) > float(50)) & (S3N2Bin_binning_result['Completeness'].astype(float) <= float(90)) & (
                    S3N2Bin_binning_result['Contamination'].astype(float) < float(contamination * 100))]

        S3N2Bin_all_bin = S3N2Bin_binning_result.index.tolist()
        S3N2Bin_hq = hq.index.tolist()
        S3N2Bin_mq = mq.index.tolist()
        S3N2Bin_lq = [temp for temp in S3N2Bin_all_bin if temp not in S3N2Bin_hq and temp not in S3N2Bin_mq]

        assert len(S3N2Bin_hq) + len(S3N2Bin_mq) + len(S3N2Bin_lq) == len(S3N2Bin_all_bin)

        Metabat2_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/Metabat2/result.csv'.format(dataset,sample),index_col=0)
        hq = Metabat2_binning_result[(Metabat2_binning_result['Completeness'].astype(float) > float(90)) & (
                Metabat2_binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                                            Metabat2_binning_result['pass.GUNC'] == True)]
        all_values = Metabat2_binning_result.values
        all_bin = Metabat2_binning_result.index.tolist()
        for bin_index, checkm_index in zip(all_bin, all_values):
            Metabat2_bin_dict[sample + '_' + bin_index] = (checkm_index[10], checkm_index[11])

        mq = Metabat2_binning_result[(Metabat2_binning_result['Completeness'].astype(float) > float(50)) & (
                    Metabat2_binning_result['Completeness'].astype(float) <= float(90)) & (
                                            Metabat2_binning_result['Contamination'].astype(float) < float(
                                        contamination * 100))]
        Metabat2_all_bin = Metabat2_binning_result.index.tolist()
        Metabat2_hq = hq.index.tolist()
        Metabat2_mq = mq.index.tolist()
        Metabat2_lq = [temp for temp in Metabat2_all_bin if
                           temp not in Metabat2_hq and temp not in Metabat2_mq]
        assert len(Metabat2_hq) + len(Metabat2_mq) + len(Metabat2_lq) == len(Metabat2_all_bin)

        data = pd.read_csv('Results/Real/Mash_dist/{0}/{1}/Mash_S3N2Bin_Metabat2.txt'.format(dataset, sample),header=None,sep='\t')
        data.columns = ['genome_id_1', 'genome_id_2', 'dis', 'p_value', 'Matching-hashes']
        data = data[data['dis'] <= 0.01].values
        for temp in data:
            S3N2Bin_bin = temp[0].split('/')[-1][:-3]
            Metabat2_bin = temp[1].split('/')[-1][:-3]
            if S3N2Bin_bin in S3N2Bin_hq and Metabat2_bin in Metabat2_hq:
                S3N2Bin_hq_list.append(sample + '_' + S3N2Bin_bin)

            if S3N2Bin_bin in S3N2Bin_hq and Metabat2_bin in Metabat2_mq:
                S3N2Bin_mq_list.append(sample + '_' + S3N2Bin_bin)

            if S3N2Bin_bin in S3N2Bin_hq and Metabat2_bin in Metabat2_lq:
                S3N2Bin_others_list.append(sample + '_' + S3N2Bin_bin)

            if S3N2Bin_bin in S3N2Bin_hq and Metabat2_bin in Metabat2_hq:
                Metabat2_hq_list.append(sample + '_' + Metabat2_bin)

            if S3N2Bin_bin in S3N2Bin_mq and Metabat2_bin in Metabat2_hq:
                Metabat2_mq_list.append(sample + '_' + Metabat2_bin)

            if S3N2Bin_bin in S3N2Bin_lq and Metabat2_bin in Metabat2_hq:
                Metabat2_others_list.append(sample + '_' + Metabat2_bin)
    return S3N2Bin_hq_list, S3N2Bin_mq_list, S3N2Bin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list, S3N2Bin_bin_dict, Metabat2_bin_dict



def plot_overlap_F1(dataset = 'dog'):
    S3N2Bin_hq_list, S3N2Bin_mq_list, S3N2Bin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list,S3N2Bin_bin_dict, Metabat2_bin_dict = get_overlap(dataset)

    S3N2Bin_recall = []
    S3N2Bin_precision = []
    S3N2Bin_F1 = []

    Metabat2_recall = []
    Metabat2_precision = []
    Metabat2_F1 = []
    data = []

    for bin in S3N2Bin_hq_list:
        recall = float(S3N2Bin_bin_dict[bin][0]) / 100
        precision = 1 - float(S3N2Bin_bin_dict[bin][1]) / 100
        F1 = 2 * recall * precision / (recall + precision)

        S3N2Bin_recall.append(recall)
        S3N2Bin_precision.append(precision)
        S3N2Bin_F1.append(F1)
        data.append(['Recall', recall, '$S^3N^2Bin$'])
        data.append(['Precision', precision, '$S^3N^2Bin$'])
        data.append(['F1-score', F1, '$S^3N^2Bin$'])

    for bin in Metabat2_hq_list:
        recall = float(Metabat2_bin_dict[bin][0]) / 100
        precision = 1 - float(Metabat2_bin_dict[bin][1]) / 100
        F1 = 2 * recall * precision / (recall + precision)

        Metabat2_recall.append(recall)
        Metabat2_precision.append(precision)
        Metabat2_F1.append(F1)

        data.append(['Recall', recall, 'Metabat2'])
        data.append(['Precision', precision, 'Metabat2'])
        data.append(['F1-score', F1, 'Metabat2'])

    data = pd.DataFrame(np.array(data), columns=['metrics', 'value', 'Method'])
    data[['value']] = data[['value']].astype(float)
    ax = sns.boxplot(x="metrics", y="value", hue="Method",
                     data=data, showfliers=False)

    ax.set_yticks(ticks=[0.9, 0.92, 0.94, 0.96, 0.98, 1.0])
    ax.set_yticklabels(labels=[0.9, 0.92, 0.94, 0.96, 0.98, 1.0], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Recall', 'Precision', 'F1-score'],
                       minor=False, fontsize=15, color='black')
    if dataset == 'dog':
        ax.set_title('{}'.format('Dog gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'human':
        ax.set_title('{}'.format('Human gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'tara':
        ax.set_title('{}'.format('Tara'), fontsize=15, alpha=1.0, color='black')
    ax.set_ylabel('Value', fontsize=15, color='black')
    ax.set_xlabel('Metrics', fontsize=15, color='black')
    plt.show()
    #plt.savefig('/home1/pansj/plot/dog_F1_distribution.jpg', dpi=300, bbox_inches='tight')

def get_taxa_list(bac_path, arr_path = None):
    bac = pd.read_csv(bac_path,
        '\t').values
    if arr_path is not None:
        arr = pd.read_csv(arr_path,
                          '\t').values
        bac = np.concatenate((bac,arr),axis=0)
    species_list = []
    genus_list = []
    family_list = []
    order_list = []
    Class_list = []
    phylum_list = []
    domain_list = []

    for temp in bac:
        split = temp[1].split(';')

        species = split[-1]
        genus = split[-2]
        family = split[-3]
        order = split[-4]
        Class = split[-5]
        phylum = split[-6]
        domain = split[-7]

        if species != 's__':
            species_list.append(species)
        if genus != 'g__':
            genus_list.append(genus)
        if family != 'f__':
            family_list.append(family)
        if order != 'o__':
            order_list.append(order)
        if Class != 'c__':
            Class_list.append(Class)
        if phylum != 'p__':
            phylum_list.append(phylum)
        if domain != 'd__':
            domain_list.append(domain)

    return set(domain_list),set(phylum_list),set(Class_list),set(order_list), set(family_list), set(genus_list), set(species_list)

def plot_all_taxi_overlap(dataset = 'dog'):
    if dataset != 'human':
        S3N2Bin_domain_list, S3N2Bin_phylum_list, S3N2Bin_Class_list, S3N2Bin_order_list, S3N2Bin_family_list, S3N2Bin_genus_list, S3N2Bin_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.bac120.summary.tsv'.format(dataset))

        Metabat2_domain_list, Metabat2_phylum_list, Metabat2_Class_list, Metabat2_order_list, Metabat2_family_list, Metabat2_genus_list, Metabat2_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
    else:
        S3N2Bin_domain_list, S3N2Bin_phylum_list, S3N2Bin_Class_list, S3N2Bin_order_list, S3N2Bin_family_list, S3N2Bin_genus_list, S3N2Bin_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.ar122.summary.tsv'.format(dataset))

        Metabat2_domain_list, Metabat2_phylum_list, Metabat2_Class_list, Metabat2_order_list, Metabat2_family_list, Metabat2_genus_list, Metabat2_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

    # Both; S3N2Bin_distinct; Metabat2_distinct
    domain = []
    domain.append(len(list(S3N2Bin_domain_list.intersection(Metabat2_domain_list))))
    domain.append(len(list(S3N2Bin_domain_list.difference(Metabat2_domain_list))))
    domain.append(len(list(Metabat2_domain_list.difference(S3N2Bin_domain_list))))

    phylum = []
    phylum.append(len(list(S3N2Bin_phylum_list.intersection(Metabat2_phylum_list))))
    phylum.append(len(list(S3N2Bin_phylum_list.difference(Metabat2_phylum_list))))
    phylum.append(len(list(Metabat2_phylum_list.difference(S3N2Bin_phylum_list))))

    Class = []
    Class.append(len(list(S3N2Bin_Class_list.intersection(Metabat2_Class_list))))
    Class.append(len(list(S3N2Bin_Class_list.difference(Metabat2_Class_list))))
    Class.append(len(list(Metabat2_Class_list.difference(S3N2Bin_Class_list))))

    order = []
    order.append(len(list(S3N2Bin_order_list.intersection(Metabat2_order_list))))
    order.append(len(list(S3N2Bin_order_list.difference(Metabat2_order_list))))
    order.append(len(list(Metabat2_order_list.difference(S3N2Bin_order_list))))

    family = []
    family.append(len(list(S3N2Bin_family_list.intersection(Metabat2_family_list))))
    family.append(len(list(S3N2Bin_family_list.difference(Metabat2_family_list))))
    family.append(len(list(Metabat2_family_list.difference(S3N2Bin_family_list))))

    genus = []
    genus.append(len(list(S3N2Bin_genus_list.intersection(Metabat2_genus_list))))
    genus.append(len(list(S3N2Bin_genus_list.difference(Metabat2_genus_list))))
    genus.append(len(list(Metabat2_genus_list.difference(S3N2Bin_genus_list))))

    species = []
    species.append(len(list(S3N2Bin_species_list.intersection(Metabat2_species_list))))
    species.append(len(list(S3N2Bin_species_list.difference(Metabat2_species_list))))
    species.append(len(list(Metabat2_species_list.difference(S3N2Bin_species_list))))

    subset = np.zeros((7, 3))
    for temp in range(3):
        subset[0][temp] = species[temp]
        subset[1][temp] = genus[temp]
        subset[2][temp] = family[temp]
        subset[3][temp] = order[temp]
        subset[4][temp] = Class[temp]
        subset[5][temp] = phylum[temp]
        subset[6][temp] = domain[temp]

    subset = pd.DataFrame(subset, index=['Species', 'Genus', 'Family', 'Order', 'Class', 'phylum', 'domain'],
                          columns=['Both', 'S3N2Bin_only', 'Metabat2_only'])
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False)

    ax.legend(['Both', '$S^3N^2$Bin_only', 'Metabat2_only'],
              loc='upper right', fontsize=10)
    ax.set_yticks(ticks=[0, 20, 40, 60, 80, 100])
    ax.set_yticklabels(labels=[0, 20, 40, 60, 80, 100], fontsize=12, color='black')

    ax.set_xticklabels(labels=['Species', 'Genus', 'Family', 'Order', 'Class', 'phylum', 'domain'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_ylabel('Taxas', fontsize=15, color='black')
    if dataset == 'dog':
        ax.set_title('{}'.format('Dog gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'human':
        ax.set_title('{}'.format('Human gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'tara':
        ax.set_title('{}'.format('Tara'), fontsize=15, alpha=1.0, color='black')
    plt.show()
    #plt.savefig('/home1/pansj/plot/dog_annotation.jpg', dpi=300, bbox_inches='tight')

def get_known_unknown(bac_path,  arr_path = None):
    bac = pd.read_csv(bac_path,
        '\t').values

    known_list = []
    unknown_list = []
    for temp in bac:
        split = temp[1].split(';')
        species = split[-1]

        if species != 's__':
            known_list.append(species)
        else:
            unknown_list.append(species)

    if arr_path is not None:
        arr = pd.read_csv(arr_path,
                          '\t').values
        for temp in arr:
            split = temp[1].split(';')
            species = split[-1]

            if species != 's__':
                known_list.append(species)
            else:
                unknown_list.append(species)

    return known_list, unknown_list

def plot_comparison_known_unknown(dataset = 'dog', y_label = None):
    if dataset != 'human':
        S3N2Bin_known, S3N2Bin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.bac120.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
    else:
        S3N2Bin_known, S3N2Bin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/S3N2Bin/gtdbtk.ar122.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

    subset = np.array([[len(Metabat2_known), len(Metabat2_unknown)], [len(S3N2Bin_known), len(S3N2Bin_unknown)]])
    print(subset)
    subset = pd.DataFrame(subset, index=['Metabat2', 'S3N2Bin'],
                          columns=['Known', 'Unknown'])
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False, figsize=(3, 4))

    ax.legend(['Known', 'Unknown'],
              loc='upper left', fontsize=10)
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')

    ax.set_xticklabels(labels=['Metabat2', '$S^3N^2$Bin'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_ylabel('Bins', fontsize=15, color='black')

    if dataset == 'dog':
        ax.set_title('{}'.format('Dog gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'human':
        ax.set_title('{}'.format('Human gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'tara':
        ax.set_title('{}'.format('Tara'), fontsize=15, alpha=1.0, color='black')
    plt.show()
    #plt.savefig('/home1/pansj/plot/human_unknown.jpg', dpi=300, bbox_inches='tight')

def plot_overlap_comparison(dataset = 'dog'):
    S3N2Bin_hq_list, S3N2Bin_mq_list, S3N2Bin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list,_ , _ = get_overlap(dataset)
    S3N2Bin_high_quality = get_num_high_quality(dataset, 'S3N2Bin')
    Metabat2_high_quality = get_num_high_quality(dataset, 'Metabat2')

    ### High quaity in both; Medium quality in other; Low quality in other; Miss in other
    S3N2Bin_result = []
    Metabat2_result = []
    S3N2Bin_result.append(len(S3N2Bin_hq_list))
    S3N2Bin_result.append(len(S3N2Bin_mq_list))
    S3N2Bin_result.append(len(S3N2Bin_others_list))
    S3N2Bin_result.append(S3N2Bin_high_quality - len(S3N2Bin_hq_list) - len(S3N2Bin_mq_list) - len(S3N2Bin_others_list))

    Metabat2_result.append(len(Metabat2_hq_list))
    Metabat2_result.append(len(Metabat2_mq_list))
    Metabat2_result.append(len(Metabat2_others_list))
    Metabat2_result.append(Metabat2_high_quality - len(Metabat2_hq_list) - len(Metabat2_mq_list) - len(Metabat2_others_list))

    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    x = np.char.array(['NC in both({})'.format(S3N2Bin_result[0]),'MQ in other({})'.format(S3N2Bin_result[1]),'LQ in other({})'.format(S3N2Bin_result[2]),'Miss in other(544)'.format(S3N2Bin_result[3])])
    y = np.array([S3N2Bin_result[0],S3N2Bin_result[1],S3N2Bin_result[2],S3N2Bin_result[3]])
    patches,texts = plt.pie(y, startangle=90,colors=colors,textprops={'fontsize': 12,'color':'black'})
    porcent = 100. * y / y.sum()
    labels = ['{0} - {1:1.2f} %'.format(i, j) for i, j in zip(x, porcent)]
    centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.axis('equal')
    plt.title('{}'.format('$S^3N^2$Bin'), fontsize=15, alpha=1.0,color = 'black')

    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.1, 1.),
               fontsize=9)
    plt.show()
    # plt.savefig('/home1/pansj/plot/dog_overlap_s3n2bin.jpg',dpi=300,bbox_inches = 'tight')
    # plt.close()

    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    x = np.char.array(['NC in both({})'.format(Metabat2_result[0]),'MQ in other({})'.format(Metabat2_result[1]),'LQ in other({})'.format(Metabat2_result[2]),'Miss in other({})'.format(Metabat2_result[3])])
    y = np.array([Metabat2_result[0],Metabat2_result[1],Metabat2_result[2],Metabat2_result[3]])
    patches,texts = plt.pie(y, startangle=90,colors=colors,textprops={'fontsize': 12,'color':'black'})
    porcent = 100. * y / y.sum()
    labels = ['{0} - {1:1.2f} %'.format(i, j) for i, j in zip(x, porcent)]
    centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.axis('equal')
    plt.title('{}'.format('Metabat2'), fontsize=15, alpha=1.0,color = 'black')

    plt.legend(patches, labels, loc='center left', bbox_to_anchor=(-0.1, 1.),
               fontsize=9)
    plt.show()
    # plt.savefig('/home1/pansj/plot/dog_overlap_metabat2.jpg',dpi=300,bbox_inches = 'tight')
    # plt.close()

def plot_sample_comparison_F1(sample):
    S3N2Bin_hq, S3N2Bin_mq, S3N2Bin_others, Metabat2_hq, Metabat2_mq, Metabat2_others,S3N2Bin_bin_dict , Metabat2_bin_dict = get_overlap('human')

    S3N2Bin_result = get_result('human','S3N2Bin')[sample]['high quality']
    Metabat2_result = get_result('human','Metabat2')[sample]['high quality']
    sample_list = [sample]
    for sample in sample_list:
        data = pd.read_csv(
            'Results/Real/FastANI/human/{}/fastANI_compare_whole'.format(sample),
            sep='\t', header=None)

        data.columns = ['genome_id_1', 'genome_id_2', 'ANI', 'A', 'B']
        data = data[data['ANI'] > 95]
        genome_group = data.groupby('genome_id_1')

        Metabat2_F1_list = []
        S3N2Bin_F1_list = []

        for genome, group in genome_group:
            group = group.values
            S3N2Bin_bin = group[0][0].split('/')[-1][:-3]
            if S3N2Bin_bin not in S3N2Bin_result:
                continue

            if len(group) != 0:
                recall = float(S3N2Bin_bin_dict[sample + '_' + S3N2Bin_bin][0]) / 100
                precision = 1 - float(S3N2Bin_bin_dict[sample + '_' + S3N2Bin_bin][1]) / 100
                S3N2Bin_F1 = 2 * recall * precision / (recall + precision)
                for temp in group:
                    Metabat2_bin = temp[1].split('/')[-1][:-3]
                    #
                    # if Metabat2_bin in Metabat2_result:
                    #     continue
                    recall = float(Metabat2_bin_dict[sample + '_' + Metabat2_bin][0]) / 100
                    precision = 1 - float(Metabat2_bin_dict[sample + '_' + Metabat2_bin][1]) / 100
                    if precision < 0:
                        precision = 0
                    Metabat2_F1 = 2 * recall * precision / (recall + precision)

                    Metabat2_F1_list.append(Metabat2_F1)
                    S3N2Bin_F1_list.append(S3N2Bin_F1)

        zipped = zip(Metabat2_F1_list, S3N2Bin_F1_list)
        sort_zipped = sorted(zipped, key=lambda x: (x[0]))
        result = zip(*sort_zipped)
        Metabat2_F1_list, S3N2Bin_F1_list = [list(x) for x in result]
        x = list((range(1, len(S3N2Bin_F1_list) + 1)))

        plt.plot(x, Metabat2_F1_list, label='Metabat2')
        plt.plot(x, S3N2Bin_F1_list, label='$S^3N^2$Bin')

        plt.xticks(color='black')
        plt.yticks(color='black')
        plt.legend()
        plt.title('Human gut({})'.format(sample), fontsize=15)

        plt.ylabel('F1-score', fontsize=15, color='black')
        plt.show()
        #plt.savefig('/home1/pansj/plot/human_gut_CCMD56948710ST.jpg', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    ### bar plot high quality genomes comparison
    plot_high_quality_comparison()

    ### venn plot multi annotation comparison
    plot_multi_venn_comparison()

    ### per sample comparison
    plot_per_sample_comparison()
    plot_per_sample_comparison('human')
    plot_per_sample_comparison('tara')

    ### recall, precision, F1-score box plot
    plot_overlap_F1()
    plot_overlap_F1('human')
    plot_overlap_F1('tara')

    ### bar plot the overlap of annotation in all taxi
    plot_all_taxi_overlap()
    plot_all_taxi_overlap(dataset='human')
    plot_all_taxi_overlap(dataset='tara')

    ## comparison of known and unknown species
    plot_comparison_known_unknown(y_label=[0,500,1000,1500,2000,2500])
    plot_comparison_known_unknown(dataset='human', y_label=[0,300,600,900,1200,1500])
    plot_comparison_known_unknown(dataset='tara', y_label=[0,100,200,300,400])

    ### compare the overlap of the high quality bins of S3N2Bin and Metabat2
    plot_overlap_comparison()
    plot_overlap_comparison(dataset='human')
    plot_overlap_comparison(dataset='tara')

    ### plot per sample comparasion of the improvements from Metabat2 to S3N2Bin(high quality) with fastANI > 95%
    plot_sample_comparison_F1('CCMD56948710ST')




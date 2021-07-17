"""
This script is used to reproduct the plot of the real datasets
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3_unweighted
from scipy.stats import wilcoxon

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

contamination = 0.05
dog_list = ['SAMN06172456', 'SAMN06172425', 'SAMN06172487', 'SAMN06172450', 'SAMN06172459', 'SAMN06172479', 'SAMN06172435', 'SAMN06172414', 'SAMN06172409', 'SAMEA103957796', 'SAMN06172442', 'SAMN06172500', 'SAMN06172437', 'SAMN06172413', 'SAMN06172514', 'SAMN06172403', 'SAMN06172471', 'SAMN06172490', 'SAMN06172448', 'SAMN06172504', 'SAMN06172457', 'SAMN06172441', 'SAMN06172422', 'SAMN06172408', 'SAMN06172429', 'SAMN06172420', 'SAMN06172503', 'SAMN06172410', 'SAMN06172458', 'SAMN06172493', 'SAMEA103957794', 'SAMN06172402', 'SAMN06172515', 'SAMN06172462', 'SAMN06172421', 'SAMN06172411', 'SAMN06172511', 'SAMN06172516', 'SAMN06172465', 'SAMN06172419', 'SAMN06172517', 'SAMN06172510', 'SAMN06172418', 'SAMN06172424', 'SAMN06172427', 'SAMN06172453', 'SAMN06172491', 'SAMN06172496', 'SAMN06172513', 'SAMN06172461', 'SAMN06172449', 'SAMN06172426', 'SAMN06172452', 'SAMN06172522', 'SAMN06172400', 'SAMN06172405', 'SAMN06172521', 'SAMN06172407', 'SAMN06172455', 'SAMN06172446', 'SAMN06172467', 'SAMN06172499', 'SAMN06172474', 'SAMN06172412', 'SAMN06172468', 'SAMN06172478', 'SAMN06172423', 'SAMN06172447', 'SAMN06172415', 'SAMN06172523', 'SAMN06172417', 'SAMN06172497', 'SAMN06172498', 'SAMN06172489', 'SAMN06172436', 'SAMN06172432', 'SAMN06172406', 'SAMN06172488', 'SAMN06172502', 'SAMN06172401', 'SAMN06172434', 'SAMN06172416', 'SAMN06172445', 'SAMN06172431', 'SAMN06172438', 'SAMN06172473', 'SAMN06172486', 'SAMN06172472', 'SAMN06172428', 'SAMEA103957793', 'SAMEA103957795', 'SAMN06172443', 'SAMN06172475', 'SAMN06172520', 'SAMN06172495', 'SAMN06172440', 'SAMN06172430', 'SAMN06172481', 'SAMN06172524', 'SAMN06172519', 'SAMN06172454', 'SAMN06172404', 'SAMN06172460', 'SAMN06172433', 'SAMN06172469', 'SAMN06172451', 'SAMN06172476', 'SAMN06172492', 'SAMN06172484', 'SAMN06172509', 'SAMN06172506', 'SAMN06172518', 'SAMN06172477', 'SAMN06172470', 'SAMN06172482', 'SAMN06172512', 'SAMN06172494', 'SAMN06172485', 'SAMN06172508', 'SAMN06172466', 'SAMN06172507', 'SAMN06172444', 'SAMN06172505', 'SAMN06172464', 'SAMN06172439', 'SAMN06172501', 'SAMN06172483', 'SAMN06172463', 'SAMN06172480']

human_list = ['CCMD41521570ST', 'CCMD75147712ST', 'CCMD18579000ST', 'CCMD53508245ST', 'CCMD19168690ST', 'CCMD52117727ST', 'CCMD42956136ST', 'CCMD79349503ST', 'CCMD89306485ST', 'CCMD76409700ST', 'CCMD31134579ST', 'CCMD71242853ST', 'CCMD89107682ST', 'CCMD76222476ST', 'CCMD10032470ST', 'CCMD17410933ST', 'CCMD38158721ST', 'CCMD35081859ST', 'CCMD54057834ST', 'CCMD28738636ST', 'CCMD98702133ST', 'CCMD30626189ST', 'CCMD32965613ST', 'CCMD53522274ST', 'CCMD37575804ST', 'CCMD68973846ST', 'CCMD25475945ST', 'CCMD65406197ST', 'CCMD21703880ST', 'CCMD50300306ST', 'CCMD51228890ST', 'CCMD59540613ST', 'CCMD49942357ST', 'CCMD95431029ST', 'CCMD41202658ST', 'CCMD15562448ST', 'CCMD21593359ST', 'CCMD92404903ST', 'CCMD50538120ST', 'CCMD49461418ST', 'CCMD72690923ST', 'CCMD85481373ST', 'CCMD39882286ST', 'CCMD18829815ST', 'CCMD51154251ST', 'CCMD85661207ST', 'CCMD71915439ST', 'CCMD39157124ST', 'CCMD22852639ST', 'CCMD35801800ST', 'CCMD27463710ST', 'CCMD59583015ST', 'CCMD89967135ST', 'CCMD52145360ST', 'CCMD95676152ST', 'CCMD45004878ST', 'CCMD67373733ST', 'CCMD99929634ST', 'CCMD89643949ST', 'CCMD26625622ST', 'CCMD23541216ST', 'CCMD31009081ST', 'CCMD99440714ST', 'CCMD66848156ST', 'CCMD65222621ST', 'CCMD98531134ST', 'CCMD45812507ST', 'CCMD46727384ST', 'CCMD73128545ST', 'CCMD30627121ST', 'CCMD50529145ST', 'CCMD98198513ST', 'CCMD93755960ST', 'CCMD35633353ST', 'CCMD56948710ST', 'CCMD27867141ST', 'CCMD32288175ST', 'CCMD29706695ST', 'CCMD72666896ST', 'CCMD10191450ST', 'CCMD49025643ST', 'CCMD74592084ST']

tara_list = ['TARA_041_SRF_0.1-0.22', 'TARA_038_SRF_0.22-1.6', 'TARA_076_SRF_0.22-3', 'TARA_023_SRF_0.22-1.6', 'TARA_042_SRF_0.22-1.6', 'TARA_124_SRF_0.22-3', 'TARA_124_SRF_0.22-0.45', 'TARA_066_SRF_lt-0.22', 'TARA_057_SRF_0.22-3', 'TARA_124_SRF_0.45-0.8', 'TARA_004_SRF_0.22-1.6', 'TARA_018_SRF_0.22-1.6', 'TARA_070_SRF_0.22-0.45', 'TARA_034_SRF_lt-0.22', 'TARA_064_SRF_0.22-3', 'TARA_125_SRF_0.22-0.45', 'TARA_111_SRF_0.22-3', 'TARA_122_SRF_0.22-0.45', 'TARA_145_SRF_0.22-3', 'TARA_099_SRF_0.22-3', 'TARA_038_SRF_lt-0.22', 'TARA_082_SRF_0.22-3', 'TARA_041_SRF_lt-0.22', 'TARA_146_SRF_0.22-3', 'TARA_151_SRF_0.22-3', 'TARA_123_SRF_0.22-3', 'TARA_110_SRF_0.22-3', 'TARA_150_SRF_0.22-3', 'TARA_072_SRF_lt-0.22', 'TARA_085_SRF_0.22-3', 'TARA_098_SRF_0.22-3', 'TARA_078_SRF_0.22-3', 'TARA_149_SRF_0.22-3', 'TARA_094_SRF_0.22-3', 'TARA_068_SRF_0.22-3', 'TARA_148_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_018_SRF_lt-0.22', 'TARA_138_SRF_0.22-3', 'TARA_093_SRF_0.22-3', 'TARA_041_SRF_0.22-1.6', 'TARA_122_SRF_0.22-3', 'TARA_078_SRF_0.45-0.8', 'TARA_070_SRF_lt-0.22', 'TARA_065_SRF_lt-0.22', 'TARA_122_SRF_0.1-0.22', 'TARA_036_SRF_0.22-1.6', 'TARA_031_SRF_0.22-1.6', 'TARA_142_SRF_0.22-3', 'TARA_124_SRF_0.1-0.22', 'TARA_036_SRF_lt-0.22', 'TARA_065_SRF_0.22-3', 'TARA_067_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_109_SRF_0.22-3', 'TARA_068_SRF_lt-0.22', 'TARA_109_SRF_lt-0.22', 'TARA_064_SRF_lt-0.22', 'TARA_048_SRF_0.22-1.6', 'TARA_034_SRF_0.22-1.6', 'TARA_070_SRF_0.45-0.8', 'TARA_025_SRF_lt-0.22', 'TARA_133_SRF_0.22-3', 'TARA_096_SRF_0.22-3', 'TARA_038_SRF_0.1-0.22', 'TARA_007_SRF_0.22-1.6', 'TARA_048_SRF_0.1-0.22', 'TARA_140_SRF_0.22-3', 'TARA_034_SRF_0.1-0.22', 'TARA_067_SRF_0.22-3', 'TARA_125_SRF_0.45-0.8', 'TARA_030_SRF_0.22-1.6', 'TARA_031_SRF_lt-0.22', 'TARA_032_SRF_0.22-1.6', 'TARA_070_SRF_0.22-3', 'TARA_132_SRF_0.22-3', 'TARA_076_SRF_lt-0.22', 'TARA_125_SRF_0.1-0.22', 'TARA_123_SRF_0.45-0.8', 'TARA_078_SRF_lt-0.22', 'TARA_068_SRF_0.45-0.8', 'TARA_068_SRF_0.22-0.45', 'TARA_067_SRF_0.22-0.45', 'TARA_100_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_137_SRF_0.22-3', 'TARA_076_SRF_0.22-0.45', 'TARA_125_SRF_0.22-3', 'TARA_078_SRF_0.22-0.45', 'TARA_076_SRF_0.45-0.8', 'TARA_084_SRF_0.22-3', 'TARA_032_SRF_lt-0.22', 'TARA_025_SRF_0.22-1.6', 'TARA_062_SRF_0.22-3', 'TARA_066_SRF_0.22-3', 'TARA_036_SRF_0.1-0.22', 'TARA_056_SRF_0.22-3', 'TARA_072_SRF_0.22-3', 'TARA_128_SRF_0.22-3', 'TARA_052_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_123_SRF_0.22-0.45', 'TARA_102_SRF_0.22-3', 'TARA_065_SRF_0.1-0.22', 'TARA_009_SRF_0.22-1.6', 'TARA_141_SRF_0.22-3', 'TARA_045_SRF_0.22-1.6', 'TARA_042_SRF_lt-0.22', 'TARA_152_SRF_0.22-3']

PRJNA504891_list = ['SRR8784379', 'SRR8180449', 'SRR8784372', 'SRR8784373', 'SRR8784385', 'SRR8784383', 'SRR8784375', 'SRR8784384', 'SRR8784360', 'SRR8784395','SRR8784376', 'SRR8784393', 'SRR8784363', 'SRR8784378', 'SRR8784361', 'SRR8784357', 'SRR8784394', 'SRR8784370', 'SRR8784397', 'SRR8784374', 'SRR8784364', 'SRR8784390', 'SRR8784381','SRR8784391', 'SRR8784358', 'SRR8784380', 'SRR8784387', 'SRR8784386', 'SRR8784365', 'SRR8784359', 'SRR8784382', 'SRR8784396','SRR8784377', 'SRR8784362', 'SRR8784389', 'SRR8180448','SRR8180450', 'SRR8784366','SRR8784371', 'SRR8784353', 'SRR8784354', 'SRR8180446', 'SRR8784388', 'SRR8784367', 'SRR8784392', 'SRR8784369', 'SRR8180447', 'SRR8784356','SRR8784355', 'SRR8784368']

PRJNA290729_list = ['SAMN03922475', 'SAMN03922449', 'SAMN03922521', 'SAMN03922488', 'SAMN03922526', 'SAMN03922468', 'SAMN03922500', 'SAMN03922512', 'SAMN03922494', 'SAMN03922450', 'SAMN03922479', 'SAMN03922484', 'SAMN03922492', 'SAMN03922470', 'SAMN03922480', 'SAMN03922505', 'SAMN03922516', 'SAMN03922527', 'SAMN03922513', 'SAMN03922472', 'SAMN03922504', 'SAMN03922523', 'SAMN03922528', 'SAMN03922510', 'SAMN03922497', 'SAMN03922518', 'SAMN03922465', 'SAMN03922517', 'SAMN03922522', 'SAMN03922458', 'SAMN03922511', 'SAMN03922531', 'SAMN03922462', 'SAMN03922457', 'SAMN03922507', 'SAMN03922509', 'SAMN03922474', 'SAMN03922529', 'SAMN03922456', 'SAMN03922506', 'SAMN03922464', 'SAMN03922453', 'SAMN03922463', 'SAMN03922466', 'SAMN03922539', 'SAMN03922503', 'SAMN03922477', 'SAMN03922495', 'SAMN03922451', 'SAMN03922538', 'SAMN03922461', 'SAMN03922532', 'SAMN03922476', 'SAMN03922469', 'SAMN03922540', 'SAMN03922533', 'SAMN03922530', 'SAMN03922536', 'SAMN03922519', 'SAMN03922471', 'SAMN03922489', 'SAMN03922524', 'SAMN03922496', 'SAMN03922467', 'SAMN03922520', 'SAMN03922483', 'SAMN03922452', 'SAMN03922508', 'SAMN03922486', 'SAMN03922473', 'SAMN03922515', 'SAMN03922455', 'SAMN03922534', 'SAMN03922490', 'SAMN03922498', 'SAMN03922525', 'SAMN03922487', 'SAMN03922482', 'SAMN03922501', 'SAMN03922537', 'SAMN03922535', 'SAMN03922485', 'SAMN03922459', 'SAMN03922499', 'SAMN03922514', 'SAMN03922454', 'SAMN03922493', 'SAMN03922478', 'SAMN03922491', 'SAMN03922481', 'SAMN03922502', 'SAMN03922460']

def get_result(dataset='dog', method='Maxbin2', binning_mode = 'single_sample', checkm_only = False):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample

    checkm_only: if just using checkm or using checkm and GUNC
    """
    if dataset == 'dog':
        sample_list = dog_list
    elif dataset == 'human':
        sample_list = human_list
    elif dataset == 'tara':
        sample_list = tara_list
    elif dataset == 'PRJNA290729':
        sample_list = PRJNA290729_list
    elif dataset == 'PRJNA504891':
        sample_list = PRJNA504891_list
    else:
        raise KeyError(f"Unknown dataset {dataset}")

    result = {}
    if method == 'VAMB' and binning_mode == 'multi_sample':
        result = {'high quality': []}
        binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/multi_sample/VAMB_multi.csv'.format(dataset),index_col=0)
        if not checkm_only:
            high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        else:
            high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(contamination * 100))]
        high_quality = high_quality.index.tolist()
        result['high quality'].extend(high_quality)


        return result
    else:
        for sample in sample_list:
            result[sample] = {'high quality':[]}
            binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/{1}/{2}/{3}/result.csv'.format(dataset,binning_mode, sample, method), index_col=0)
            if not checkm_only:
                high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                        binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                                                          binning_result['pass.GUNC'] == True)]
            else:
                high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                        binning_result['Contamination'].astype(float) < float(contamination * 100))]
            high_quality = high_quality.index.tolist()
            result[sample]['high quality'].extend(high_quality)

        return result



def get_num_high_quality(dataset = 'dog', method = 'Maxbin2', binning_mode = 'single_sample', checkm_only = False):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample
    """
    result = get_result(dataset, method, binning_mode, checkm_only)
    if method == 'VAMB' and binning_mode == 'multi_sample':
        num_hq = len(result['high quality'])
        return num_hq
    else:
        num_hq = 0
        for sample in result:
            num_hq += len(result[sample]['high quality'])
        return num_hq

def tranfer_multi():
    x1 = [1,2,3,4,5,6]
    y_human_multi = [144,168.4,178.2,188.6,190.2,191.2]
    y_human_err_max = [9, 9.6, 3.8, 6.4, 6.8, 6.8]
    y_human_err_min = [13, 4.4, 3.2, 8.6, 4.2, 6.2]
    y_human_semibin = [165] * 6
    y_human_metabat2 = [133] * 6

    x2 = [1,2,3,4,5,6]
    y_dog_multi = [165.8, 182.2,190.6,199,202.4,204.6]
    y_dog_err_max = [7.2, 8.8, 7.4, 7, 2.6, 5.4]
    y_dog_err_min = [7.8, 6.2, 7.6, 10, 1.4, 7.6]
    y_dog_semibin = [194] * 6
    y_dog_metabat2 = [120] * 6

    x3 = [1,2,3,4,5,6]
    y_tara_multi = [21.6,23.4,24.6,23.8,25.2,25.6]
    y_tara_err_max = [1.4, 0.6, 1.4, 1.2, 1.8, 2.4]
    y_tara_err_min = [1.6, 0.4, 1.6, 0.8, 1.2, 1.6]
    y_tara_semibin = [26] * 6
    y_tara_metabat2 = [20] * 6

    fig, ax = plt.subplots(figsize = (2,4))
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    line_width = 1
    plt.errorbar(x1,y_human_multi,yerr=[y_human_err_min, y_human_err_max],label='SemiBin(pre-train)',color='#034e7b',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_human_semibin, label='SemiBin',color='#0570b0',linewidth=line_width)
    plt.plot(x1, y_human_metabat2, label='Metabat2',color='#74a9cf',linewidth=line_width)

    plt.errorbar(x1,y_dog_multi,yerr=[y_dog_err_min, y_dog_err_max],label='SemiBin(pre-train)',color='#005a32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_dog_semibin, label='SemiBin',color='#238443',linewidth=line_width)
    plt.plot(x1, y_dog_metabat2, label='Metabat2',color='#78c679',linewidth=line_width)

    plt.errorbar(x1,y_tara_multi,yerr=[y_tara_err_min, y_tara_err_max],label='SemiBin(pre-train)',color='#8c2d04',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_tara_semibin, label='SemiBin',color='#cc4c02',linewidth=line_width)
    plt.plot(x1, y_tara_metabat2, label='Metabat2',color='#fe9929',linewidth=line_width)

    ax.set_xticks([1,2,3,4,5,6])
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_xlabel('Number of samples', fontsize=15, color='black')
    ax.set_xticklabels([1,3,5,10,15,20],color = 'black')
    plt.yticks([0,25,50,75,100,125,150,175,200], color='black')

    plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize=8)
    plt.savefig('transfer_multiple.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_high_quality_comparison():
    num_dog_maxbin2_single = get_num_high_quality()
    num_dog_vamb_single = get_num_high_quality(method='VAMB')
    num_dog_metabat2_single = get_num_high_quality(method='Metabat2')
    num_dog_semibin_single = get_num_high_quality(method='S3N2Bin')
    num_dog_semibin_pretrain_single = get_num_high_quality(method='SemiBin_pretrain')

    num_human_maxbin2_single = get_num_high_quality(dataset='human')
    num_human_vamb_single = get_num_high_quality(dataset='human', method='VAMB')
    num_human_metabat2_single = get_num_high_quality(dataset='human', method='Metabat2')
    num_human_semibin_single = get_num_high_quality(dataset='human', method='S3N2Bin')
    num_human_semibin_pretrain_single = get_num_high_quality(dataset='human', method='SemiBin_pretrain')

    num_tara_maxbin2_single = get_num_high_quality(dataset='tara')
    num_tara_vamb_single = get_num_high_quality(dataset='tara', method='VAMB')
    num_tara_metabat2_single = get_num_high_quality(dataset='tara', method='Metabat2')
    num_tara_semibin_single = get_num_high_quality(dataset='tara', method='S3N2Bin')
    num_tara_semibin_pretrain_single = get_num_high_quality(dataset='tara', method='SemiBin_pretrain')

    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample')
    num_dog_semibin_multi = get_num_high_quality(method='S3N2Bin',binning_mode='multi_sample')

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample')
    num_human_semibin_multi = get_num_high_quality(dataset='human', method='S3N2Bin',binning_mode='multi_sample')

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample')
    num_tara_semibin_multi = get_num_high_quality(dataset='tara', method='S3N2Bin', binning_mode='multi_sample')

    # print(num_dog_maxbin2_single, num_dog_vamb_single, num_dog_metabat2_single, num_dog_semibin_single,
    #       num_dog_semibin_pretrain_single)
    # print(num_human_maxbin2_single, num_human_vamb_single, num_human_metabat2_single, num_human_semibin_single,
    #       num_human_semibin_pretrain_single)
    # print(num_tara_maxbin2_single, num_tara_vamb_single, num_tara_metabat2_single, num_tara_semibin_single,
    #       num_tara_semibin_pretrain_single)
    #
    # print(num_dog_vamb_mulit,num_dog_semibin_multi)
    # print(num_human_vamb_multi, num_human_semibin_multi)
    # print(num_tara_vamb_multi,num_tara_semibin_multi)

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single, num_dog_semibin_pretrain_single,num_dog_vamb_mulit,num_dog_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000,3500])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500,3000,3500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_dog_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)
    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single,num_human_vamb_multi,num_human_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Human gut'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,300,600,900,1200,1500])
    ax.set_yticklabels(labels=[0,300,600,900,1200,1500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality bins', fontsize=15,color = 'black')
    plt.savefig('Real_human_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)
    subset = pd.DataFrame(np.array([[num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single,num_tara_semibin_pretrain_single, num_tara_vamb_multi,num_tara_semibin_multi]]),columns = ['Maxbin2','VAMB(single)','Metabat2','SemiBin(single)','SemiBin(pre-train)','VAMB(multi)', 'SemiBin(multi)'], index=['Tara'])
    ax = subset.plot(kind='bar',figsize=(2,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0','#33a02c','#1f78b4'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15,color = 'black',rotation = 360)
    ax.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize = 8)
    plt.savefig('Real_tara_hq.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()
    print(subset)

def plot_checkm_high_quality_comparison():
    num_dog_maxbin2_single = get_num_high_quality(checkm_only=True)
    num_dog_vamb_single = get_num_high_quality(method='VAMB',checkm_only=True)
    num_dog_metabat2_single = get_num_high_quality(method='Metabat2',checkm_only=True)
    num_dog_semibin_single = get_num_high_quality(method='S3N2Bin',checkm_only=True)
    num_dog_semibin_pretrain_single = get_num_high_quality(method='SemiBin_pretrain',checkm_only=True)

    num_human_maxbin2_single = get_num_high_quality(dataset='human',checkm_only=True)
    num_human_vamb_single = get_num_high_quality(dataset='human', method='VAMB',checkm_only=True)
    num_human_metabat2_single = get_num_high_quality(dataset='human', method='Metabat2',checkm_only=True)
    num_human_semibin_single = get_num_high_quality(dataset='human', method='S3N2Bin',checkm_only=True)
    num_human_semibin_pretrain_single = get_num_high_quality(dataset='human', method='SemiBin_pretrain',checkm_only=True)

    num_tara_maxbin2_single = get_num_high_quality(dataset='tara',checkm_only=True)
    num_tara_vamb_single = get_num_high_quality(dataset='tara', method='VAMB',checkm_only=True)
    num_tara_metabat2_single = get_num_high_quality(dataset='tara', method='Metabat2',checkm_only=True)
    num_tara_semibin_single = get_num_high_quality(dataset='tara', method='S3N2Bin',checkm_only=True)
    num_tara_semibin_pretrain_single = get_num_high_quality(dataset='tara', method='SemiBin_pretrain',checkm_only=True)


    print(num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single, num_dog_semibin_pretrain_single)
    print(num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single)
    print(num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single,num_tara_semibin_pretrain_single)

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single,num_dog_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0'])
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500,3000],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_dog_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Human gut'])
    ax = subset.plot(kind='bar',figsize=(3,4), color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0'])
    ax.set_yticks(ticks=[0,400,800,1200,1600])
    ax.set_yticklabels(labels=[0,400,800,1200,1600],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality bins', fontsize=15,color = 'black')
    ax.legend(loc = 'upper left',fontsize = 8)
    plt.savefig('Real_human_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single, num_tara_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Tara'])
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3','#1d91c0'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Tara'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_tara_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_dog_semibin_multi = get_num_high_quality(method='S3N2Bin',binning_mode='multi_sample',checkm=True)

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_human_semibin_multi = get_num_high_quality(dataset='human', method='S3N2Bin',binning_mode='multi_sample',checkm_only=True)

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample',checkm_only=True)
    num_tara_semibin_multi = get_num_high_quality(dataset='tara', method='S3N2Bin', binning_mode='multi_sample',checkm_only=True)

    print(num_dog_vamb_mulit,num_dog_semibin_multi)
    print(num_human_vamb_multi, num_human_semibin_multi)
    print(num_tara_vamb_multi,num_tara_semibin_multi)

    subset = pd.DataFrame(np.array([[num_dog_vamb_mulit,num_dog_semibin_multi]]),columns = ['VAMB','SemiBin'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize = (2,4), color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3'])
    ax.set_yticks(ticks=[0,800,1600,2400,3200,4000])
    ax.set_yticklabels(labels=[0,800,1600,2400,3200,4000],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.legend(loc='lower left', fontsize=6)
    plt.savefig('Real_dog_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_vamb_multi, num_human_semibin_multi]]), columns=['VAMB', 'SemiBin'],
                          index=['Human gut'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3'])
    ax.set_yticks(ticks=[0, 400, 800, 1200, 1600, 2000])
    ax.set_yticklabels(labels=[0, 400, 800, 1200, 1600, 2000], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15, color='black', rotation=360)
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    plt.savefig('Real_human_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_vamb_multi,num_tara_semibin_multi]]), columns=['VAMB', 'SemiBin'],
                          index=['Tara'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3'])
    ax.set_yticks(ticks=[0, 150, 300, 450, 600, 750])
    ax.set_yticklabels(labels=[0, 150, 300, 450, 600, 750], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15, color='black', rotation=360)
    plt.savefig('Real_tara_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
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
    dog_SemiBin_family, dog_SemiBin_genus, dog_SemiBin_species = get_taxi_list(base_path + 'multi_sample/dog/S3N2Bin/gtdbtk.bac120.summary.tsv')
    dog_SemiBin_family_single, dog_SemiBin_genus_single, dog_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/dog/S3N2Bin/gtdbtk.bac120.summary.tsv')

    ### human
    human_VAMB_family, human_VAMB_genus, human_VAMB_species = get_taxi_list(base_path + 'multi_sample/human/VAMB/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/VAMB/gtdbtk.ar122.summary.tsv')
    human_SemiBin_family, human_SemiBin_genus, human_SemiBin_species = get_taxi_list(base_path + 'multi_sample/human/S3N2Bin/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/S3N2Bin/gtdbtk.ar122.summary.tsv')
    human_SemiBin_family_single, human_SemiBin_genus_single, human_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/human/S3N2Bin/gtdbtk.bac120.summary.tsv',base_path + 'single_sample/human/S3N2Bin/gtdbtk.ar122.summary.tsv')

    ### tara
    tara_VAMB_family, tara_VAMB_genus, tara_VAMB_species = get_taxi_list(base_path + 'multi_sample/tara/VAMB/gtdbtk.bac120.summary.tsv')
    tara_SemiBin_family, tara_SemiBin_genus, tara_SemiBin_species = get_taxi_list(base_path + 'multi_sample/tara/S3N2Bin/gtdbtk.bac120.summary.tsv')
    tara_SemiBin_family_single, tara_SemiBin_genus_single, tara_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/tara/S3N2Bin/gtdbtk.bac120.summary.tsv')

    out = venn3_unweighted([set(dog_VAMB_family), set(dog_SemiBin_family), set(dog_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Dog gut)", fontsize=20, alpha=1.0, color='black')
    plt.savefig('multi_sample_dog_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_genus), set(dog_SemiBin_genus), set(dog_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_species), set(dog_SemiBin_species), set(dog_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_family), set(human_SemiBin_family), set(human_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_genus), set(human_SemiBin_genus), set(human_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_species), set(human_SemiBin_species), set(human_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_family), set(tara_SemiBin_family), set(tara_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_genus), set(tara_SemiBin_genus), set(tara_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_species), set(tara_SemiBin_species), set(tara_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#1f78b4','#33a02c','#ff7f00'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

def plot_per_sample_comparison(dataset = 'dog', output = None):
    """
    dataset: dog, human, tara
    """
    SemiBin_result = get_result(dataset, method='SemiBin_pretrain')
    Metabat2_result = get_result(dataset, method='Metabat2')

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormap = ['#084594','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7']
    newcmp = LinearSegmentedColormap.from_list('cmps', colormap)
    result = {}

    num_SemiBin = []
    num_Metabat2 = []

    for sample in SemiBin_result:
        num_SemiBin.append(len(SemiBin_result[sample]['high quality']))
        num_Metabat2.append(len(Metabat2_result[sample]['high quality']))
    print(wilcoxon(num_Metabat2, num_SemiBin))
    for sample in SemiBin_result:
        if (len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality'])) not in result:
            result[(len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] = 1
        else:
            result[(len(SemiBin_result[sample]['high quality']), len(Metabat2_result[sample]['high quality']))] += 1

    if dataset != 'tara':
        data = np.zeros(shape=(len(result), 3))

        for i, temp in enumerate(result):
            data[i][0] = temp[0]
            data[i][1] = temp[1]
            data[i][2] = result[temp]

        data = pd.DataFrame(data).astype(int)

        data.columns = ['x', 'y', 'num']

        plt.scatter(data.x, data.y,
                    c=data.num, s=(data.num ** 2) * 60, cmap=newcmp)

        plt.colorbar(shrink=0.5)
        plt.xlabel("SemiBin")
        plt.ylabel("Metabat2")
        if dataset == 'dog':
            plt.plot([0, 35], [0, 35], c='black')
            plt.title('Dog gut'.format(dataset), fontsize=15)
        else:
            plt.plot([0, 40], [0, 40], c='black')
            plt.title('Human gut'.format(dataset), fontsize=15)
        plt.show()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()


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
                    c=data.num, s=(data.num) * 60, cmap=newcmp)

        plt.colorbar()
        plt.xlabel("SemiBin")
        plt.ylabel("Metabat2")
        plt.plot([0, 25], [0, 25], c='black')
        plt.title('Ocean', fontsize=15)
        plt.show()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()

def get_overlap(dataset = 'dog'):
    if dataset == 'dog':
        sample_list = dog_list
    if dataset == 'human':
        sample_list = human_list
    if dataset == 'tara':
        sample_list = tara_list

    SemiBin_hq_list = []
    SemiBin_mq_list = []
    SemiBin_others_list = []

    Metabat2_hq_list = []
    Metabat2_mq_list = []
    Metabat2_others_list = []

    SemiBin_bin_dict = {}
    Metabat2_bin_dict = {}

    for sample in sample_list:
        SemiBin_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/SemiBin_pretrain/result.csv'.format(dataset,sample),index_col=0)

        all_values = SemiBin_binning_result.values
        all_bin = SemiBin_binning_result.index.tolist()
        for bin_index, checkm_index in zip(all_bin, all_values):
            SemiBin_bin_dict[sample + '_' + bin_index] = (checkm_index[10], checkm_index[11])
        SemiBin_all_bin = SemiBin_binning_result.index.tolist()
        hq = SemiBin_binning_result[(SemiBin_binning_result['Completeness'].astype(float) > float(90)) & (SemiBin_binning_result['Contamination'].astype(float) < float(contamination * 100)) & (SemiBin_binning_result['pass.GUNC'] == True)]
        SemiBin_hq = hq.index.tolist()

        mq = SemiBin_binning_result[(SemiBin_binning_result['Completeness'].astype(float) >= float(50)) & ( SemiBin_binning_result['Contamination'].astype(float) < float(0.1 * 100))]
        SemiBin_mq = [temp for temp in mq.index.tolist() if temp not in SemiBin_hq]

        SemiBin_others = [temp for temp in SemiBin_all_bin if temp not in SemiBin_hq and temp not in SemiBin_mq]

        assert len(SemiBin_hq) + len(SemiBin_mq) + len(SemiBin_others) == len(SemiBin_all_bin)

        Metabat2_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/Metabat2/result.csv'.format(dataset,sample),index_col=0)

        all_values = Metabat2_binning_result.values
        all_bin = Metabat2_binning_result.index.tolist()
        Metabat2_all_bin = Metabat2_binning_result.index.tolist()

        for bin_index, checkm_index in zip(all_bin, all_values):
            Metabat2_bin_dict[sample + '_' + bin_index] = (checkm_index[10], checkm_index[11])

        hq = Metabat2_binning_result[(Metabat2_binning_result['Completeness'].astype(float) > float(90)) & (
                Metabat2_binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                                            Metabat2_binning_result['pass.GUNC'] == True)]
        Metabat2_hq = hq.index.tolist()

        mq = Metabat2_binning_result[(Metabat2_binning_result['Completeness'].astype(float) >= float(50)) & (
                                            Metabat2_binning_result['Contamination'].astype(float) < float(0.1 * 100))]
        Metabat2_mq = [temp for temp in mq.index.tolist() if temp not in Metabat2_hq]

        Metabat2_others = [temp for temp in Metabat2_all_bin if
                       temp not in Metabat2_hq and temp not in Metabat2_mq]


        assert len(Metabat2_hq) + len(Metabat2_mq) + len(Metabat2_others) == len(Metabat2_all_bin)

        data = pd.read_csv('Results/Real/Mash_dist/{0}/{1}/Mash_SemiBin_pretrain_Metabat2.txt'.format(dataset, sample),header=None,sep='\t')
        data.columns = ['genome_id_1', 'genome_id_2', 'dis', 'p_value', 'Matching-hashes']
        data = data[data['dis'] <= 0.01].values
        for temp in data:
            SemiBin_bin = temp[0].split('/')[-1][:-3]
            Metabat2_bin = temp[1].split('/')[-1][:-3]
            if SemiBin_bin in SemiBin_hq and Metabat2_bin in Metabat2_hq:
                SemiBin_hq_list.append(sample + '_' + SemiBin_bin)

            if SemiBin_bin in SemiBin_hq and Metabat2_bin in Metabat2_mq:
                SemiBin_mq_list.append(sample + '_' + SemiBin_bin)

            if SemiBin_bin in SemiBin_hq and Metabat2_bin in Metabat2_others:
                SemiBin_others_list.append(sample + '_' + SemiBin_bin)

            if SemiBin_bin in SemiBin_hq and Metabat2_bin in Metabat2_hq:
                Metabat2_hq_list.append(sample + '_' + Metabat2_bin)

            if SemiBin_bin in SemiBin_mq and Metabat2_bin in Metabat2_hq:
                Metabat2_mq_list.append(sample + '_' + Metabat2_bin)

            if SemiBin_bin in SemiBin_others and Metabat2_bin in Metabat2_hq:
                Metabat2_others_list.append(sample + '_' + Metabat2_bin)
    return SemiBin_hq_list, SemiBin_mq_list, SemiBin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list, SemiBin_bin_dict, Metabat2_bin_dict

def plot_sankey_overlap(dataset = 'dog', output = None):
    SemiBin_hq_list, SemiBin_mq_list, SemiBin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list,_ , _ = get_overlap(dataset)

    SemiBin_high_quality = get_num_high_quality(dataset, 'SemiBin_pretrain')
    Metabat2_high_quality = get_num_high_quality(dataset, 'Metabat2')
    SemiBin_other = SemiBin_high_quality - len(SemiBin_hq_list) - len(SemiBin_mq_list) - len(SemiBin_others_list)
    Metabat2_other = Metabat2_high_quality - len(Metabat2_hq_list) - len(Metabat2_mq_list) - len(Metabat2_others_list)

    import plotly.graph_objects as go

    unique_list = ["HQ", "MQ", "LQ", "Miss", "HQ", "MQ", "LQ", "Miss"]
    sources = [0, 0, 0, 0, 4, 4, 4, 4]
    targets = [4, 5, 6, 7, 0, 1, 2, 3]
    values = [len(Metabat2_hq_list), len(Metabat2_mq_list), len(Metabat2_others_list), Metabat2_other,len(SemiBin_hq_list), len(SemiBin_mq_list), len(SemiBin_others_list), SemiBin_other, ]
    print(len(Metabat2_hq_list), len(Metabat2_mq_list), len(Metabat2_others_list), Metabat2_other,len(SemiBin_hq_list), len(SemiBin_mq_list), len(SemiBin_others_list), SemiBin_other,)

    print(Metabat2_high_quality - len(Metabat2_hq_list), (Metabat2_high_quality - len(Metabat2_hq_list)) / Metabat2_high_quality)

    print(SemiBin_high_quality - len(SemiBin_hq_list), (SemiBin_high_quality - len(SemiBin_hq_list)) / SemiBin_high_quality)

    layout = go.Layout(autosize=True, margin={'l': 0, 'r': 0, 't': 0, 'b': 0})
    # plotly setup
    if dataset != 'dog':
        x_ = [0.1, 0.1, 0.1, 0.1, 0.4, 0.4, 0.4, 0.4]
        y_ = [0.15, 0.2, 0.3, 0.4, 0.1, 0.3, 0.4, 0.48]
    else:
        x_ = [0.1, 0.1, 0.1, 0.1, 0.4, 0.4, 0.4, 0.4]
        y_ = [0.15, 0.2, 0.3, 0.4, 0.1, 0.3, 0.4, 0.52]
    fig = go.Figure(layout=layout, data=[go.Sankey(
        arrangement='snap',
        orientation="h",
        valueformat=".0f",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0),
            label=unique_list,
            color=['#fe9929', '#fec44f', '#fee391', '#fff7bc', '#fe9929', '#fec44f', '#fee391', '#fff7bc'],
            # x=nodified[0],
            # y=nodified[1]
            x = x_,
            y = y_
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=['#fe9929', '#fec44f', '#fee391', '#fff7bc','#fe9929', '#fec44f', '#fee391', '#fff7bc']
        ))])


    fig.write_image(output)

def plot_overlap_F1(dataset = 'dog'):
    SemiBin_hq_list, SemiBin_mq_list, SemiBin_others_list, Metabat2_hq_list, Metabat2_mq_list, Metabat2_others_list,SemiBin_bin_dict, Metabat2_bin_dict = get_overlap(dataset)

    SemiBin_recall = []
    SemiBin_precision = []
    SemiBin_F1 = []

    Metabat2_recall = []
    Metabat2_precision = []
    Metabat2_F1 = []
    data = []
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

    for bin in SemiBin_hq_list:
        recall = float(SemiBin_bin_dict[bin][0]) / 100
        precision = 1 - float(SemiBin_bin_dict[bin][1]) / 100
        F1 = 2 * recall * precision / (recall + precision)

        SemiBin_recall.append(recall)
        SemiBin_precision.append(precision)
        SemiBin_F1.append(F1)
        data.append(['Recall', recall, 'SemiBin'])
        data.append(['Precision', precision, 'SemiBin'])
        data.append(['F1-score', F1, 'SemiBin'])

    print('Recall:', wilcoxon(Metabat2_recall,SemiBin_recall))
    print('Precision:', wilcoxon(Metabat2_precision,SemiBin_precision))
    print('F1:', wilcoxon(Metabat2_F1,SemiBin_F1))

    data = pd.DataFrame(np.array(data), columns=['metrics', 'value', 'Method'])
    data[['value']] = data[['value']].astype(float)

    ax = sns.boxplot(x="metrics", y="value", hue="Method",
                     data=data, fliersize=0, palette= ['#ec7014','#1b9e77'],linewidth=0.1)
    # sns.stripplot(data=data, x="metrics",dodge=True, y="value", hue="Method", size=2, palette= ['#a6cee3','#1f78b4'],)



    ax.set_yticks(ticks=[0.9, 0.92, 0.94, 0.96, 0.98, 1.0])
    ax.set_yticklabels(labels=[0.9, 0.92, 0.94, 0.96, 0.98, 1.0], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Recall', 'Precision', 'F1-score'],
                       minor=False, fontsize=15, color='black')
    if dataset == 'dog':
        ax.set_title('{}'.format('Dog gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'human':
        ax.set_title('{}'.format('Human gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'tara':
        ax.set_title('{}'.format('Ocean'), fontsize=15, alpha=1.0, color='black')
    # ax.set_ylabel('Value', fontsize=15, color='black')
    # ax.set_xlabel('Metrics', fontsize=15, color='black')
    plt.title('HQ in both')
    plt.show()
    plt.savefig('{}_F1_distribution.pdf'.format(dataset), dpi=300, bbox_inches='tight')
    plt.close()

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

def plot_all_taxi_overlap(dataset = 'dog', output = None,y_label = None):
    if dataset != 'human':
        SemiBin_domain_list, SemiBin_phylum_list, SemiBin_Class_list, SemiBin_order_list, SemiBin_family_list, SemiBin_genus_list, SemiBin_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset))

        Metabat2_domain_list, Metabat2_phylum_list, Metabat2_Class_list, Metabat2_order_list, Metabat2_family_list, Metabat2_genus_list, Metabat2_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
    else:
        SemiBin_domain_list, SemiBin_phylum_list, SemiBin_Class_list, SemiBin_order_list, SemiBin_family_list, SemiBin_genus_list, SemiBin_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.ar122.summary.tsv'.format(dataset))

        Metabat2_domain_list, Metabat2_phylum_list, Metabat2_Class_list, Metabat2_order_list, Metabat2_family_list, Metabat2_genus_list, Metabat2_species_list = get_taxa_list('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

    # Both; S3N2Bin_distinct; Metabat2_distinct
    domain = []
    domain.append(len(list(SemiBin_domain_list.intersection(Metabat2_domain_list))))
    domain.append(len(list(SemiBin_domain_list.difference(Metabat2_domain_list))))
    domain.append(len(list(Metabat2_domain_list.difference(SemiBin_domain_list))))

    phylum = []
    phylum.append(len(list(SemiBin_phylum_list.intersection(Metabat2_phylum_list))))
    phylum.append(len(list(SemiBin_phylum_list.difference(Metabat2_phylum_list))))
    phylum.append(len(list(Metabat2_phylum_list.difference(SemiBin_phylum_list))))

    Class = []
    Class.append(len(list(SemiBin_Class_list.intersection(Metabat2_Class_list))))
    Class.append(len(list(SemiBin_Class_list.difference(Metabat2_Class_list))))
    Class.append(len(list(Metabat2_Class_list.difference(SemiBin_Class_list))))

    order = []
    order.append(len(list(SemiBin_order_list.intersection(Metabat2_order_list))))
    order.append(len(list(SemiBin_order_list.difference(Metabat2_order_list))))
    order.append(len(list(Metabat2_order_list.difference(SemiBin_order_list))))

    family = []
    family.append(len(list(SemiBin_family_list.intersection(Metabat2_family_list))))
    family.append(len(list(SemiBin_family_list.difference(Metabat2_family_list))))
    family.append(len(list(Metabat2_family_list.difference(SemiBin_family_list))))

    genus = []
    genus.append(len(list(SemiBin_genus_list.intersection(Metabat2_genus_list))))
    genus.append(len(list(SemiBin_genus_list.difference(Metabat2_genus_list))))
    genus.append(len(list(Metabat2_genus_list.difference(SemiBin_genus_list))))

    species = []
    species.append(len(list(SemiBin_species_list.intersection(Metabat2_species_list))))
    species.append(len(list(SemiBin_species_list.difference(Metabat2_species_list))))
    species.append(len(list(Metabat2_species_list.difference(SemiBin_species_list))))

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
                          columns=['Both', 'SemiBin(pre-train) only', 'Metabat2 only'])
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False, color = ['#a6cee3','#fdbf6f', '#b2df8a', ])

    ax.legend(['Both', 'SemiBin(pre-train) only', 'Metabat2 only'],
              loc='upper right', fontsize=10)
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')

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
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


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

def plot_comparison_known_unknown(dataset = 'dog', y_label = None,output = None):
    if dataset != 'human':
        SemiBin_known, SemiBin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
    else:
        SemiBin_known, SemiBin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.ar122.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

    subset = np.array([[len(Metabat2_known), len(Metabat2_unknown)], [len(SemiBin_known), len(SemiBin_unknown)]])
    print(subset)
    subset = pd.DataFrame(subset, index=['Metabat2', 'SemiBin(pre-train)'],
                          columns=['Known', 'Unknown'])
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False, figsize=(3, 4), color = ['#a6cee3','#fdbf6f', '#b2df8a', ])

    ax.legend(['Known', 'Unknown'],
              loc='upper left', fontsize=10)
    ax.set_yticks(ticks=y_label)
    ax.set_yticklabels(labels=y_label, fontsize=12, color='black')

    ax.set_xticklabels(labels=['Metabat2', 'SemiBin'], rotation=50,
                       minor=False, fontsize=15, color='black')
    ax.set_ylabel('Number of species', fontsize=15, color='black')

    if dataset == 'dog':
        ax.set_title('{}'.format('Dog gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'human':
        ax.set_title('{}'.format('Human gut'), fontsize=15, alpha=1.0, color='black')
    if dataset == 'tara':
        ax.set_title('{}'.format('Tara'), fontsize=15, alpha=1.0, color='black')
    plt.show()
    plt.show()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    #plt.savefig('/home1/pansj/plot/human_unknown.jpg', dpi=300, bbox_inches='tight')

def plot_transfer():
    df = pd.DataFrame([[131, 170, 58, 144, 165, 159, 131, 122, 145, 99, 103, 124],
                       [120, 194, 72, 134, 150, 166, 175, 166, 191, 127, 122, 156],
                       [35, 58, 26, 38, 41, 43, 37, 41, 49, 49, 56, 52]], index=['Human gut', 'Dog gut', 'Tara'],
                      columns=['Metabat2', 'SemiBin', 'No_Semi', 'Human_low', 'Human_medium', 'Human_high',
                               'Dog_low', 'Dog_medium', 'Dog_high', 'Tara_low', 'Tara_medium', 'Tara_high', ])

    df = df[['Tara_high', 'Tara_medium', 'Tara_low', 'Dog_high', 'Dog_medium', 'Dog_low', 'Human_high', 'Human_medium',
             'Human_low', 'No_Semi', 'SemiBin', 'Metabat2']].T
    print(df)
    fig, ax = plt.subplots()
    N = df.index.size

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormap = ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7']
    colormap.reverse()
    newcmp1 = LinearSegmentedColormap.from_list('cmps', colormap)

    colormap = ['#005a32', '#238b45', '#41ab5d', '#74c476', '#a1d99b', '#c7e9c0', '#e5f5e0', ]
    colormap.reverse()
    newcmp2 = LinearSegmentedColormap.from_list('cmps', colormap)

    colormap = ['#8c2d04', '#d94801', '#f16913', '#fd8d3c', '#fdae6b', '#fdd0a2', '#fee6ce', ]
    colormap.reverse()
    newcmp3 = LinearSegmentedColormap.from_list('cmps', colormap)

    # first heatmap
    im1 = ax.imshow(np.vstack([df['Human gut'].astype(float)]).T, aspect='auto', extent=[-0.5, 0.5, -0.5, N - 0.5],
                    origin='lower', cmap=newcmp1)
    # second heatmap
    im2 = ax.imshow(np.vstack([df['Dog gut'].astype(float)]).T, aspect='auto', extent=[0.5, 1.5, -0.5, N - 0.5],
                    origin='lower', cmap=newcmp2)

    im3 = ax.imshow(np.vstack([df['Tara'].astype(float)]).T, aspect='auto', extent=[1.5, 2.5, -0.5, N - 0.5],
                    origin='lower', cmap=newcmp3)

    cbar1 = fig.colorbar(im1, ax=ax, label='Human gut')
    cbar2 = fig.colorbar(im2, ax=ax, label='Dog gut')
    cbar3 = fig.colorbar(im3, ax=ax, label='Tara')

    ax.set_xlim(-0.5, 2.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Human gut', 'Dog gut', 'Tara'])
    ax.set_yticks(range(N))
    ax.set_yticklabels(df.index)
    # ax.set_ylabel('Countries')
    fig.tight_layout()
    plt.show()
    plt.savefig('transfer_results.pdf', dpi=300, bbox_inches='tight')

def plot_extra_bar():
    result_Metabat2_non_western = get_result('PRJNA504891', 'Metabat2')
    result_SemiBin_non_western = get_result('PRJNA504891', 'SemiBin')
    result_SemiBin_pretrain_same_non_western = get_result('PRJNA504891', 'SemiBin_pretrain_same')
    result_SemiBin_pretrain_out_non_western = get_result('PRJNA504891', 'SemiBin_pretrain_out')

    result_Metabat2_western = get_result('PRJNA290729', 'Metabat2')
    result_SemiBin_western = get_result('PRJNA290729', 'SemiBin')
    result_SemiBin_pretrain_same_western = get_result('PRJNA290729', 'SemiBin_pretrain_same')
    result_SemiBin_pretrain_out_western = get_result('PRJNA290729', 'SemiBin_pretrain_out')

    hq_Metabat2_non_western = []
    hq_SemiBin_non_western = []
    hq_SemiBin_pretrain_same_non_western = []
    hq_SemiBin_pretrain_out_non_western = []

    hq_Metabat2_western = []
    hq_SemiBin_western  = []
    hq_SemiBin_pretrain_same_western = []
    hq_SemiBin_pretrain_out_western = []

    for sample in PRJNA504891_list:
        hq_Metabat2_non_western.append(len(result_Metabat2_non_western[sample]['high quality']))
        hq_SemiBin_non_western.append(len(result_SemiBin_non_western[sample]['high quality']))
        hq_SemiBin_pretrain_same_non_western.append(len(result_SemiBin_pretrain_same_non_western[sample]['high quality']))
        hq_SemiBin_pretrain_out_non_western.append(len(result_SemiBin_pretrain_out_non_western[sample]['high quality']))

    for sample in PRJNA290729_list:
        hq_Metabat2_western.append(len(result_Metabat2_western[sample]['high quality']))
        hq_SemiBin_western.append(len(result_SemiBin_western[sample]['high quality']))
        hq_SemiBin_pretrain_same_western.append(len(result_SemiBin_pretrain_same_western[sample]['high quality']))
        hq_SemiBin_pretrain_out_western.append(len(result_SemiBin_pretrain_out_western[sample]['high quality']))
    print('non-western:', wilcoxon(hq_Metabat2_non_western, hq_SemiBin_pretrain_out_non_western))
    print(len(hq_Metabat2_non_western))
    print('western:', wilcoxon(hq_Metabat2_western, hq_SemiBin_pretrain_out_western))
    print(len(hq_Metabat2_western))
    print((np.sum(hq_SemiBin_pretrain_out_western) - np.sum(hq_Metabat2_western)) / np.sum(hq_Metabat2_western))
    print((np.sum(hq_SemiBin_pretrain_out_non_western) - np.sum(hq_Metabat2_non_western)) / np.sum(hq_Metabat2_non_western))


    subset = pd.DataFrame(np.array([[np.sum(hq_Metabat2_non_western),np.sum(hq_SemiBin_non_western), np.sum(hq_SemiBin_pretrain_same_non_western), np.sum(hq_SemiBin_pretrain_out_non_western)]]),columns = ['Metabat2','SemiBin','SemiBin(pre-train; internal)','SemiBin(pre-train; external)'], index=['Non-western human gut'])
    ax = subset.plot(kind='bar',width = 0.3,color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Non-western human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality bins', fontsize=15,color = 'black')
    #ax.legend(loc='lower left', fontsize=6)
    plt.show()
    plt.savefig('holdout_Non_western_bar.pdf', dpi=300, bbox_inches='tight')
    print(subset)
    subset = pd.DataFrame(np.array([[np.sum(hq_Metabat2_western),np.sum(hq_SemiBin_western), np.sum(hq_SemiBin_pretrain_same_western), np.sum(hq_SemiBin_pretrain_out_western)]]),columns = ['Metabat2','SemiBin','SemiBin(pre-train; internal)','SemiBin(pre-train; external)'], index=['Western human gut'])
    ax = subset.plot(kind='bar',width = 0.3,color=['#fb9a99','#b2df8a','#fdbf6f','#a6cee3'])
    ax.set_yticks(ticks=[0,100,200,300,400,500,600])
    ax.set_yticklabels(labels=[0,100,200,300,400,500,600],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Western human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality bins', fontsize=15,color = 'black')
    #ax.legend(loc='lower left', fontsize=6)
    plt.show()
    plt.savefig('holdout_western_bar.pdf', dpi=300, bbox_inches='tight')
    print(subset)


def plot_extra_per_sample():
    result_Metabat2_non_western = get_result('PRJNA504891', 'Metabat2')
    result_SemiBin_pretrain_out_non_western = get_result('PRJNA504891', 'SemiBin_pretrain_out')

    result_Metabat2_western = get_result('PRJNA290729', 'Metabat2')
    result_SemiBin_pretrain_out_western = get_result('PRJNA290729', 'SemiBin_pretrain_out')

    hq_Metabat2_non_western = []
    hq_SemiBin_pretrain_out_non_western = []

    hq_Metabat2_western = []
    hq_SemiBin_pretrain_out_western = []

    for sample in PRJNA504891_list:
        hq_Metabat2_non_western.append(len(result_Metabat2_non_western[sample]['high quality']))
        hq_SemiBin_pretrain_out_non_western.append(len(result_SemiBin_pretrain_out_non_western[sample]['high quality']))

    for sample in PRJNA290729_list:
        hq_Metabat2_western.append(len(result_Metabat2_western[sample]['high quality']))
        hq_SemiBin_pretrain_out_western.append(len(result_SemiBin_pretrain_out_western[sample]['high quality']))

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormap = ['#084594','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7']
    newcmp = LinearSegmentedColormap.from_list('cmps', colormap)
    data_non_western = np.zeros(shape=(len(hq_Metabat2_non_western), 3))
    data_western = np.zeros(shape=(len(hq_Metabat2_western), 3))

    per_sample_result_non_Western = {}
    per_sample_result_Western = {}

    for i in range(len(hq_Metabat2_non_western)):
        if (hq_SemiBin_pretrain_out_non_western[i], hq_Metabat2_non_western[i]) not in per_sample_result_non_Western:
            per_sample_result_non_Western[(hq_SemiBin_pretrain_out_non_western[i], hq_Metabat2_non_western[i])] = 1
        else:
            per_sample_result_non_Western[(hq_SemiBin_pretrain_out_non_western[i], hq_Metabat2_non_western[i])] += 1

    for i, temp in enumerate(per_sample_result_non_Western):
        data_non_western[i][0] = temp[0]
        data_non_western[i][1] = temp[1]
        data_non_western[i][2] = per_sample_result_non_Western[temp]
    data_non_western = pd.DataFrame(data_non_western).astype(int)
    data_non_western.columns = ['x', 'y', 'num']


    for i in range(len(hq_Metabat2_western)):
        if (hq_SemiBin_pretrain_out_western[i], hq_Metabat2_western[i]) not in per_sample_result_Western:
            per_sample_result_Western[(hq_SemiBin_pretrain_out_western[i], hq_Metabat2_western[i])] = 1
        else:
            per_sample_result_Western[(hq_SemiBin_pretrain_out_western[i], hq_Metabat2_western[i])] += 1

    for i, temp in enumerate(per_sample_result_Western):
        data_western[i][0] = temp[0]
        data_western[i][1] = temp[1]
        data_western[i][2] = per_sample_result_Western[temp]
    data_western = pd.DataFrame(data_western).astype(int)
    data_western.columns = ['x', 'y', 'num']

    plt.scatter(data_non_western.x, data_non_western.y,
                c=data_non_western.num, s=(data_non_western.num ** 2) * 60, cmap=newcmp)

    plt.colorbar(shrink=0.5)
    plt.xlabel("SemiBin(pretrain; external)")
    plt.ylabel("Metabat2")

    plt.plot([0, 35], [0, 35], c='black')
    plt.title('Non-western human gut', fontsize=15)
    plt.show()
    plt.savefig('holdout_non_western_persample.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    plt.scatter(data_western.x, data_western.y,
                c=data_western.num, s=(data_western.num ** 2) * 60, cmap=newcmp)

    plt.colorbar(shrink=0.5)
    plt.xlabel("SemiBin(pretrain; external)")
    plt.ylabel("Metabat2")

    plt.plot([0, 15], [0, 15], c='black')
    plt.title('Western human gut', fontsize=15)
    plt.show()
    plt.savefig('holdout_western_persample.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def CAT_mmseqs():
    subset = pd.DataFrame(np.array([[1320,1497],[1834,2415]]),columns = ['CAT','mmseqs'], index=['Human gut','Dog gut'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#b2df8a', '#a6cee3'],figsize=(3,4))
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut','Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Single-sample binning', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_CAT_mmseqs_single.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[1652,1968],[2934,3448]]),columns = ['CAT','mmseqs'], index=['Human gut','Dog gut'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#b2df8a', '#a6cee3'],figsize=(3,4))
    ax.set_yticks(ticks=[0,800,1600,2400,3200])
    ax.set_yticklabels(labels=[0,800,1600,2400,3200],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut','Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High quality genomes', fontsize=15,color = 'black')
    ax.set_title('Multi-sample binning', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_CAT_mmseqs_multi.pdf', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':


    #tranfer_multi()
    #
    # ## bar plot high quality genomes comparison
    #plot_high_quality_comparison()
    #plot_checkm_high_quality_comparison()
    # ### venn plot multi annotation comparison
    # plot_multi_venn_comparison()
    #
    # ### per sample comparison
    # plot_per_sample_comparison(output='dog_compare_persample.pdf')
    # plot_per_sample_comparison('human',output='human_compare_persample.pdf')
    # plot_per_sample_comparison('tara',output='tara_compare_persample.pdf')
    #
    #
    plot_sankey_overlap(dataset='human',output='human_sankey.pdf')
    plot_sankey_overlap(output='dog_sankey.pdf')
    plot_sankey_overlap(dataset='tara',output='tara_sankey.pdf')
    #
    #
    ### recall, precision, F1-score box plot
    # plot_overlap_F1()
    # plot_overlap_F1('human')
    # plot_overlap_F1('tara')
    #
    # ### bar plot the overlap of annotation in all taxi
    # plot_all_taxi_overlap(output='dog_taxi_overlap.pdf',y_label=[0,20,40,60,80,100])
    # plot_all_taxi_overlap(dataset='human',output='human_taxi_overlap.pdf',y_label=[0,100,200,300,400])
    # plot_all_taxi_overlap(dataset='tara',output='tara_taxi_overlap.pdf',y_label=[0,50,100,150,200,250])
    #
    # ## comparison of known and unknown species
    # plot_comparison_known_unknown(y_label=[0,500,1000,1500,2000,2500], output='dog_taxi_known_unknown.pdf')
    # plot_comparison_known_unknown(dataset='human', y_label=[0,300,600,900,1200,1500], output='human_taxi_known_unknown.pdf')
    # plot_comparison_known_unknown(dataset='tara', y_label=[0,100,200,300,400], output='tara_taxi_known_unknown.pdf')
    #
    #plot_transfer()
    #
    # plot_extra_bar()
    #
    # plot_extra_per_sample()
    #
    #
    # CAT_mmseqs()

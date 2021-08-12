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

transfer_single_human = ['CCMD41521570ST', 'CCMD21593359ST', 'CCMD89107682ST', 'CCMD18579000ST', 'CCMD76222476ST', 'CCMD22852639ST', 'CCMD99440714ST', 'CCMD38158721ST', 'CCMD50300306ST', 'CCMD98198513ST']

transfer_single_dog = ['SAMN06172505', 'SAMN06172402', 'SAMN06172415', 'SAMN06172428', 'SAMN06172503', 'SAMN06172429', 'SAMN06172448', 'SAMN06172456', 'SAMN06172449', 'SAMN06172418']

transfer_single_tara = ['TARA_125_SRF_0.22-3', 'TARA_102_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_122_SRF_0.1-0.22', 'TARA_066_SRF_lt-0.22', 'TARA_067_SRF_0.22-3', 'TARA_123_SRF_0.45-0.8', 'TARA_093_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_038_SRF_0.22-1.6']

transfer_multi_human = ['CCMD54057834ST', 'CCMD99440714ST', 'CCMD10191450ST', 'CCMD92404903ST', 'CCMD53522274ST', 'CCMD50300306ST', 'CCMD89107682ST', 'CCMD95431029ST', 'CCMD27463710ST', 'CCMD49025643ST']

transfer_multi_dog = ['SAMN06172506', 'SAMN06172436', 'SAMN06172403', 'SAMN06172520', 'SAMN06172512', 'SAMN06172433', 'SAMN06172493', 'SAMN06172422', 'SAMN06172523', 'SAMN06172407']

transfer_multi_tara = ['TARA_045_SRF_0.22-1.6', 'TARA_041_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_140_SRF_0.22-3', 'TARA_138_SRF_0.22-3', 'TARA_034_SRF_lt-0.22', 'TARA_066_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_056_SRF_0.22-3', 'TARA_034_SRF_0.22-1.6']

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
    elif dataset == 'dog_transfer_single':
        sample_list = transfer_single_dog
        dataset = 'dog'
    elif dataset == 'human_transfer_single':
        sample_list = transfer_single_human
        dataset = 'human'
    elif dataset == 'tara_transfer_single':
        sample_list = transfer_single_tara
        dataset = 'tara'
    elif dataset == 'dog_transfer_multi':
        sample_list = transfer_multi_dog
        dataset = 'dog'
    elif dataset == 'human_transfer_multi':
        sample_list = transfer_multi_human
        dataset = 'human'
    elif dataset == 'tara_transfer_multi':
        sample_list = transfer_multi_tara
        dataset = 'tara'
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


def get_results_table(dataset, method, binning_mode='single_sample', checkm_only=False):
    r = get_result(dataset,
            method,
            binning_mode,
            checkm_only)
    if method == 'VAMB' and binning_mode == 'multi_sample':
        from collections import Counter
        hq = r['high quality']
        r = pd.Series(Counter(b.rsplit('C', 1)[0] for b in hq)).reset_index().rename(columns={'index':'sample', 0: 'nr_hq'})
    else:
        r = pd.DataFrame([(k,len(v['high quality'])) for k,v in r.items()], columns=['sample', 'nr_hq'])
    if binning_mode == 'multi_sample':
        method += '_multi'
    r['dataset'] = dataset
    r['method'] = method
    r['binning_mode'] = binning_mode
    r['checkm_only'] = checkm_only
    return r


def get_num_high_quality(dataset = 'dog', method = 'Maxbin2', binning_mode = 'single_sample', checkm_only = False):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample
    """
    return get_results_table(dataset, method, binning_mode, checkm_only)['nr_hq'].sum()


def get_num_high_quality_pretrain(dataset, num_run, num_sample):
    if dataset == 'human':
        sample_list = transfer_multi_human
    if dataset == 'dog':
        sample_list = transfer_multi_dog
    if dataset == 'tara':
        sample_list = transfer_multi_tara
    num_high_quality = 0
    for sample in sample_list:
        binning_result = pd.read_csv('Results/Real/CheckM_GUNC/Pretrain/{3}/{0}/{1}/{2}/result.csv'.format(num_run, num_sample, sample,dataset),index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                binning_result['Contamination'].astype(float) < float(contamination * 100)) & (
                                                  binning_result['pass.GUNC'] == True)]
        num_high_quality += high_quality.shape[0]
    return num_high_quality


def tranfer_multi():
    #human
    human_result = {1:[],3:[],5:[],10:[],15:[],20:[]}
    dog_result = {1:[],3:[],5:[],10:[],15:[],20:[]}
    tara_result = {1:[],3:[],5:[],10:[],15:[],20:[]}


    for j in [1,3,5,10,15,20]:
        for i in [1, 2, 3, 4, 5]:
            human_result[j].append(get_num_high_quality_pretrain('human', i, j))
            dog_result[j].append(get_num_high_quality_pretrain('dog', i, j))
            tara_result[j].append(get_num_high_quality_pretrain('tara', i, j))

    human_SemiBin = get_num_high_quality('human_transfer_multi', 'SemiBin')
    human_Metabat2 = get_num_high_quality('human_transfer_multi', 'Metabat2')

    dog_SemiBin = get_num_high_quality('dog_transfer_multi', 'SemiBin')
    dog_Metabat2 = get_num_high_quality('dog_transfer_multi', 'Metabat2')

    tara_SemiBin = get_num_high_quality('tara_transfer_multi', 'SemiBin')
    tara_Metabat2 = get_num_high_quality('tara_transfer_multi', 'Metabat2')

    x1 = [1, 2, 3, 4, 5, 6]
    y_human_multi = []
    y_human_err_max  = []
    y_human_err_min = []
    y_human_semibin = [human_SemiBin] * 6
    y_human_metabat2 = [human_Metabat2] * 6

    y_dog_multi = []
    y_dog_err_max = []
    y_dog_err_min = []
    y_dog_semibin = [dog_SemiBin] * 6
    y_dog_metabat2 = [dog_Metabat2] * 6

    y_tara_multi = []
    y_tara_err_max = []
    y_tara_err_min = []
    y_tara_semibin = [tara_SemiBin] * 6
    y_tara_metabat2 = [tara_Metabat2] * 6

    for i in [1,3,5,10,15,20]:
        y_human_multi.append(np.mean(human_result[i]))
        y_human_err_max.append(np.max(human_result[i]) - np.mean(human_result[i]))
        y_human_err_min.append(np.mean(human_result[i]) - np.min(human_result[i]))

        y_dog_multi.append(np.mean(dog_result[i]))
        y_dog_err_max.append(np.max(dog_result[i]) - np.mean(dog_result[i]))
        y_dog_err_min.append(np.mean(dog_result[i]) - np.min(dog_result[i]))

        y_tara_multi.append(np.mean(tara_result[i]))
        y_tara_err_max.append(np.max(tara_result[i]) - np.mean(tara_result[i]))
        y_tara_err_min.append(np.mean(tara_result[i]) - np.min(tara_result[i]))

    fig, ax = plt.subplots(figsize = (2,4))
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    line_width = 1
    plt.errorbar(x1,y_human_multi,yerr=[y_human_err_min, y_human_err_max],label='SemiBin(pretrain)',color='#034e7b',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_human_semibin, label='SemiBin',color='#0570b0',linewidth=line_width)
    plt.plot(x1, y_human_metabat2, label='Metabat2',color='#74a9cf',linewidth=line_width)

    plt.errorbar(x1,y_dog_multi,yerr=[y_dog_err_min, y_dog_err_max],label='SemiBin(pretrain)',color='#005a32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_dog_semibin, label='SemiBin',color='#238443',linewidth=line_width)
    plt.plot(x1, y_dog_metabat2, label='Metabat2',color='#78c679',linewidth=line_width)

    plt.errorbar(x1,y_tara_multi,yerr=[y_tara_err_min, y_tara_err_max],label='SemiBin(pretrain)',color='#8c2d04',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
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


def plot_bar_per_sample_com(dataset, diff_label = None, num_label = None):
    num_single = pd.concat([
            get_results_table(dataset=dataset, method='Maxbin2'),
            get_results_table(dataset=dataset, method='VAMB'),
            get_results_table(dataset=dataset, method='Metabat2'),
            get_results_table(dataset=dataset, method='SemiBin'),
            get_results_table(dataset=dataset, method='SemiBin_pretrain'),
        ])

    counts_single = pd.pivot(num_single[['sample', 'nr_hq', 'method']], values=['nr_hq'], index='sample', columns='method')
    counts_single  = counts_single.T.xs('nr_hq').T
    num_multi = pd.concat([
            get_results_table(dataset=dataset, method='VAMB', binning_mode='multi_sample'),
            get_results_table(dataset=dataset, method='SemiBin', binning_mode='multi_sample'),
        ])
    counts_multi = pd.pivot(num_multi[['sample', 'nr_hq', 'method']], values=['nr_hq'], index='sample', columns='method')
    counts_multi  = counts_multi.T.xs('nr_hq').T
    counts_multi = counts_multi.fillna(0)

    counts = pd.concat((counts_single,counts_multi),axis=1)

    Maxbin2 = counts['Maxbin2'].values
    VAMB = counts['VAMB'].values
    Metabat2 = counts['Metabat2'].values
    SemiBin = counts['SemiBin'].values
    SemiBin_pretrain = counts['SemiBin_pretrain'].values
    print('Compared to Maxbin2 wilcoxon:', wilcoxon(Maxbin2, SemiBin_pretrain))
    print('Compared to VAMB wilcoxon:', wilcoxon(VAMB, SemiBin_pretrain))
    print('Compared to Metabat2 wilcoxon:', wilcoxon(Metabat2, SemiBin_pretrain))

    print('SemiBin compared to Metabat2: {0}({1}) improvement'.format(np.sum(SemiBin) - np.sum(Metabat2), (np.sum(SemiBin) - np.sum(Metabat2))/ np.sum(Metabat2)))
    print('SemiBin_pretrain compared to Metabat2: {0}({1}) improvement'.format(np.sum(SemiBin_pretrain) - np.sum(Metabat2), (np.sum(SemiBin_pretrain) - np.sum(Metabat2))/ np.sum(Metabat2)))
    print('SemiBin_pretrain compared to SemiBin: {0}({1}) improvement'.format(np.sum(SemiBin_pretrain) - np.sum(SemiBin), (np.sum(SemiBin_pretrain) - np.sum(SemiBin))/ np.sum(SemiBin)))

    VAMB_multi = counts['VAMB_multi'].values
    SemiBin_multi = counts['SemiBin_multi'].values
    print('SemiBin_multi compared to VAMB_multi: {0}({1}) improvement'.format(np.sum(SemiBin_multi) - np.sum(VAMB_multi), (np.sum(SemiBin_multi) - np.sum(VAMB_multi))/ np.sum(VAMB_multi)))
    print('Compared to VAMB_multi wilcoxon:', wilcoxon(VAMB_multi, SemiBin_multi))

    diff_single = counts_single.T - counts_single.SemiBin_pretrain

    diff_multi = counts_multi.T - counts_multi.SemiBin_multi

    diff = pd.concat((diff_single,diff_multi))
    fig,axes = plt.subplots(1, 2, sharex='col',figsize = (9,3))

    val = diff.T
    val = val[['Maxbin2', 'VAMB', 'Metabat2', 'SemiBin', 'SemiBin_pretrain','VAMB_multi','SemiBin_multi',]]
    ax = sns.swarmplot(y='method', x='value', data=pd.melt(val), ax=axes[1], size=2, palette=['#e7298a','#6a51a3','#ec7014','#41AB5D','#005A32','#9e9ac8','#a1d99b'])
    if diff_label is not None:
        ax.set_xticks(ticks=diff_label)
        ax.set_xticklabels(labels=diff_label)

    v = pd.DataFrame({'total': counts.sum()})
    v = v.reindex(index = ['Maxbin2',  'VAMB','Metabat2','SemiBin', 'SemiBin_pretrain','VAMB_multi','SemiBin_multi',])
    ax = sns.barplot(data=v.reset_index(), x='total', y='method', ax=axes[0], palette=['#e7298a','#6a51a3','#ec7014','#41AB5D','#005A32','#9e9ac8','#a1d99b'])

    if num_label is not None:
        ax.set_xticks(ticks=num_label)
        ax.set_xticklabels(labels=num_label)

    ax.set_yticklabels(labels=[])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.despine()
    fig.tight_layout()
    fig.savefig('bar_per_sample_comparsion_plot_{}.pdf'.format(dataset), dpi=300)


def plot_checkm_high_quality_comparison():
    num_dog_maxbin2_single = get_num_high_quality(checkm_only=True)
    num_dog_vamb_single = get_num_high_quality(method='VAMB',checkm_only=True)
    num_dog_metabat2_single = get_num_high_quality(method='Metabat2',checkm_only=True)
    num_dog_semibin_single = get_num_high_quality(method='SemiBin',checkm_only=True)
    num_dog_semibin_pretrain_single = get_num_high_quality(method='SemiBin_pretrain',checkm_only=True)

    num_human_maxbin2_single = get_num_high_quality(dataset='human',checkm_only=True)
    num_human_vamb_single = get_num_high_quality(dataset='human', method='VAMB',checkm_only=True)
    num_human_metabat2_single = get_num_high_quality(dataset='human', method='Metabat2',checkm_only=True)
    num_human_semibin_single = get_num_high_quality(dataset='human', method='SemiBin',checkm_only=True)
    num_human_semibin_pretrain_single = get_num_high_quality(dataset='human', method='SemiBin_pretrain',checkm_only=True)

    num_tara_maxbin2_single = get_num_high_quality(dataset='tara',checkm_only=True)
    num_tara_vamb_single = get_num_high_quality(dataset='tara', method='VAMB',checkm_only=True)
    num_tara_metabat2_single = get_num_high_quality(dataset='tara', method='Metabat2',checkm_only=True)
    num_tara_semibin_single = get_num_high_quality(dataset='tara', method='SemiBin',checkm_only=True)
    num_tara_semibin_pretrain_single = get_num_high_quality(dataset='tara', method='SemiBin_pretrain',checkm_only=True)


    print(num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single, num_dog_semibin_pretrain_single)
    print(num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single)
    print(num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single,num_tara_semibin_pretrain_single)

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single,num_dog_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Dog gut'])
    print(subset)
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500,3000],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_dog_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Human gut'])
    print(subset)
    ax = subset.plot(kind='bar',figsize=(3,4), color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,400,800,1200,1600])
    ax.set_yticklabels(labels=[0,400,800,1200,1600],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.legend(loc = 'upper left',fontsize = 8)
    plt.savefig('Real_human_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single, num_tara_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Tara'])
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_tara_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_dog_semibin_multi = get_num_high_quality(method='SemiBin',binning_mode='multi_sample',checkm_only=True)

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_human_semibin_multi = get_num_high_quality(dataset='human', method='SemiBin',binning_mode='multi_sample',checkm_only=True)

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample',checkm_only=True)
    num_tara_semibin_multi = get_num_high_quality(dataset='tara', method='SemiBin', binning_mode='multi_sample',checkm_only=True)

    print(num_dog_vamb_mulit,num_dog_semibin_multi)
    print(num_human_vamb_multi, num_human_semibin_multi)
    print(num_tara_vamb_multi,num_tara_semibin_multi)

    subset = pd.DataFrame(np.array([[num_dog_vamb_mulit,num_dog_semibin_multi]]),columns = ['VAMB','SemiBin'], index=['Dog gut'])
    ax = subset.plot(kind='bar',figsize = (2,4), color=['#41AB5D','#7570b3'])
    ax.set_yticks(ticks=[0,800,1600,2400,3200,4000])
    ax.set_yticklabels(labels=[0,800,1600,2400,3200,4000],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.legend(loc='lower left', fontsize=6)
    plt.savefig('Real_dog_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_vamb_multi, num_human_semibin_multi]]), columns=['VAMB', 'SemiBin'],
                          index=['Human gut'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#41AB5D','#7570b3'])
    ax.set_yticks(ticks=[0, 400, 800, 1200, 1600, 2000])
    ax.set_yticklabels(labels=[0, 400, 800, 1200, 1600, 2000], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Human gut'], fontsize=15, color='black', rotation=360)
    ax.set_ylabel('High-quality bins', fontsize=15, color='black')
    plt.savefig('Real_human_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_tara_vamb_multi,num_tara_semibin_multi]]), columns=['VAMB', 'SemiBin'],
                          index=['Tara'])
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#41AB5D','#7570b3'])
    ax.set_yticks(ticks=[0, 150, 300, 450, 600, 750])
    ax.set_yticklabels(labels=[0, 150, 300, 450, 600, 750], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15, color='black', rotation=360)
    plt.savefig('Real_tara_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    print('Dog single ', (num_dog_semibin_pretrain_single - num_dog_metabat2_single)/num_dog_metabat2_single)
    print('human single ', (num_human_semibin_pretrain_single - num_human_metabat2_single) / num_human_metabat2_single)
    print('ocean single ', (num_tara_semibin_pretrain_single - num_tara_metabat2_single) / num_tara_metabat2_single)

    print('Dog multi', (num_dog_semibin_multi - num_dog_vamb_mulit)/num_dog_vamb_mulit)
    print('human multi', (num_human_semibin_multi - num_human_vamb_multi) / num_human_vamb_multi)
    print('tara multi', (num_tara_semibin_multi - num_tara_vamb_multi) / num_tara_vamb_multi)


def get_taxi_list(bac_path,  arr_path = None):
    bac = pd.read_csv(bac_path, '\t').values
    print(bac.shape)
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
        print(arr.shape)
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
    print('dog')
    dog_VAMB_family, dog_VAMB_genus, dog_VAMB_species = get_taxi_list(base_path + 'multi_sample/dog/VAMB/gtdbtk.bac120.summary.tsv')
    dog_SemiBin_family, dog_SemiBin_genus, dog_SemiBin_species = get_taxi_list(base_path + 'multi_sample/dog/SemiBin/gtdbtk.bac120.summary.tsv')
    dog_SemiBin_family_single, dog_SemiBin_genus_single, dog_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/dog/SemiBin/gtdbtk.bac120.summary.tsv')

    ### human
    print('human')
    human_VAMB_family, human_VAMB_genus, human_VAMB_species = get_taxi_list(base_path + 'multi_sample/human/VAMB/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/VAMB/gtdbtk.ar122.summary.tsv')
    human_SemiBin_family, human_SemiBin_genus, human_SemiBin_species = get_taxi_list(base_path + 'multi_sample/human/SemiBin/gtdbtk.bac120.summary.tsv',base_path + 'multi_sample/human/SemiBin/gtdbtk.ar122.summary.tsv')
    human_SemiBin_family_single, human_SemiBin_genus_single, human_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/human/SemiBin/gtdbtk.bac120.summary.tsv',base_path + 'single_sample/human/SemiBin/gtdbtk.ar122.summary.tsv')

    ### tara
    print('tara')
    tara_VAMB_family, tara_VAMB_genus, tara_VAMB_species = get_taxi_list(base_path + 'multi_sample/tara/VAMB/gtdbtk.bac120.summary.tsv')
    tara_SemiBin_family, tara_SemiBin_genus, tara_SemiBin_species = get_taxi_list(base_path + 'multi_sample/tara/SemiBin/gtdbtk.bac120.summary.tsv')
    tara_SemiBin_family_single, tara_SemiBin_genus_single, tara_SemiBin_species_single = get_taxi_list(base_path + 'single_sample/tara/SemiBin/gtdbtk.bac120.summary.tsv')

    out = venn3_unweighted([set(dog_VAMB_family), set(dog_SemiBin_family), set(dog_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Dog gut)", fontsize=20, alpha=1.0, color='black')
    plt.savefig('multi_sample_dog_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_genus), set(dog_SemiBin_genus), set(dog_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(dog_VAMB_species), set(dog_SemiBin_species), set(dog_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_family), set(human_SemiBin_family), set(human_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_genus), set(human_SemiBin_genus), set(human_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(human_VAMB_species), set(human_SemiBin_species), set(human_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_family), set(tara_SemiBin_family), set(tara_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_genus), set(tara_SemiBin_genus), set(tara_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    out = venn3_unweighted([set(tara_VAMB_species), set(tara_SemiBin_species), set(tara_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()


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
    # print(len(SemiBin_hq_list), len(set(SemiBin_hq_list)))
    # print(len(SemiBin_mq_list), len(set(SemiBin_mq_list)))
    # print(len(SemiBin_others_list), len(set(SemiBin_others_list)))
    # print(len(Metabat2_hq_list), len(set(Metabat2_hq_list)))
    # print(len(Metabat2_mq_list), len(set(Metabat2_mq_list)))
    # print(len(Metabat2_others_list), len(set(Metabat2_others_list)))

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
        data.append(['Recall', recall, 'SemiBin(pretrain)'])
        data.append(['Precision', precision, 'SemiBin(pretrain)'])
        data.append(['F1-score', F1, 'SemiBin(pretrain)'])

    print('Recall:', wilcoxon(Metabat2_recall,SemiBin_recall))
    print('Precision:', wilcoxon(Metabat2_precision,SemiBin_precision))
    print('F1:', wilcoxon(Metabat2_F1,SemiBin_F1))

    data = pd.DataFrame(np.array(data), columns=['metrics', 'value', 'Method'])
    data[['value']] = data[['value']].astype(float)

    plt.figure(figsize=(4, 4))
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
    print(bac.shape)
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
    # print(len(SemiBin_domain_list), len(SemiBin_phylum_list), len(SemiBin_Class_list), len(SemiBin_order_list), len(SemiBin_family_list), len(SemiBin_genus_list), len(SemiBin_species_list))
    # print(len(Metabat2_domain_list), len(Metabat2_phylum_list), len(Metabat2_Class_list), len(Metabat2_order_list), len(Metabat2_family_list), len(Metabat2_genus_list), len(Metabat2_species_list))

    # Both; S3N2Bin_distinct; Metabat2_distinct
    domain = []
    domain.append(len(list(SemiBin_domain_list.intersection(Metabat2_domain_list))))
    domain.append(len(list(SemiBin_domain_list.difference(Metabat2_domain_list))))
    domain.append(len(list(Metabat2_domain_list.difference(SemiBin_domain_list))))
    # print(len(SemiBin_domain_list))
    # print(len(Metabat2_domain_list))
    phylum = []
    phylum.append(len(list(SemiBin_phylum_list.intersection(Metabat2_phylum_list))))
    phylum.append(len(list(SemiBin_phylum_list.difference(Metabat2_phylum_list))))
    phylum.append(len(list(Metabat2_phylum_list.difference(SemiBin_phylum_list))))
    domain.append(len(list(Metabat2_domain_list.difference(SemiBin_domain_list))))
    # print(len(SemiBin_phylum_list))
    # print(len(Metabat2_phylum_list))
    Class = []
    Class.append(len(list(SemiBin_Class_list.intersection(Metabat2_Class_list))))
    Class.append(len(list(SemiBin_Class_list.difference(Metabat2_Class_list))))
    Class.append(len(list(Metabat2_Class_list.difference(SemiBin_Class_list))))
    # print(len(SemiBin_Class_list))
    # print(len(Metabat2_Class_list))
    order = []
    order.append(len(list(SemiBin_order_list.intersection(Metabat2_order_list))))
    order.append(len(list(SemiBin_order_list.difference(Metabat2_order_list))))
    order.append(len(list(Metabat2_order_list.difference(SemiBin_order_list))))
    # print(len(SemiBin_order_list))
    # print(len(Metabat2_order_list))
    family = []
    family.append(len(list(SemiBin_family_list.intersection(Metabat2_family_list))))
    family.append(len(list(SemiBin_family_list.difference(Metabat2_family_list))))
    family.append(len(list(Metabat2_family_list.difference(SemiBin_family_list))))
    # print(len(SemiBin_family_list))
    # print(len(Metabat2_family_list))
    genus = []
    genus.append(len(list(SemiBin_genus_list.intersection(Metabat2_genus_list))))
    genus.append(len(list(SemiBin_genus_list.difference(Metabat2_genus_list))))
    genus.append(len(list(Metabat2_genus_list.difference(SemiBin_genus_list))))
    # print(len(SemiBin_genus_list))
    # print(len(Metabat2_genus_list))
    species = []
    species.append(len(list(SemiBin_species_list.intersection(Metabat2_species_list))))
    species.append(len(list(SemiBin_species_list.difference(Metabat2_species_list))))
    species.append(len(list(Metabat2_species_list.difference(SemiBin_species_list))))
    # print(len(SemiBin_species_list))
    # print(len(Metabat2_species_list))
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
                          columns=['Both', 'SemiBin(pretrain) only', 'Metabat2 only'])
    print(subset)
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False, color = ['#7570b3','#1b9e77', '#ec7014', ])

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
    print(bac.shape)
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
        print(arr.shape)
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
    subset = pd.DataFrame(subset, index=['Metabat2', 'SemiBin(pretrain)'],
                          columns=['Known', 'Unknown'])
    print(subset)
    ax = subset.plot(kind="bar", stacked=True,
                     legend=False, figsize=(3, 4), color = ['#1b9e77', '#ec7014',])

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
    human_result = []
    dog_result = []
    tara_result = []

    ### human
    human_result.append(get_num_high_quality('human_transfer_single', 'Metabat2'))
    human_result.append(get_num_high_quality('human_transfer_single', 'SemiBin'))

    for method in ['no_semi', 'human_low', 'human_medium', 'human_high', 'dog_low', 'dog_medium', 'dog_high', 'tara_low',
                   'tara_medium', 'tara_high']:
        human_result.append(get_num_high_quality('human_transfer_single', method))

    dog_result.append(get_num_high_quality('dog_transfer_single', 'Metabat2'))
    dog_result.append(get_num_high_quality('dog_transfer_single', 'SemiBin'))

    for method in ['no_semi', 'human_low', 'human_medium', 'human_high', 'dog_low', 'dog_medium', 'dog_high', 'tara_low',
                   'tara_medium', 'tara_high']:
        dog_result.append(get_num_high_quality('dog_transfer_single', method))

    tara_result.append(get_num_high_quality('tara_transfer_single', 'Metabat2'))
    tara_result.append(get_num_high_quality('tara_transfer_single', 'SemiBin'))

    for method in ['no_semi', 'human_low', 'human_medium', 'human_high', 'dog_low', 'dog_medium', 'dog_high', 'tara_low',
                   'tara_medium', 'tara_high']:
        tara_result.append(get_num_high_quality('tara_transfer_single', method))

    df = pd.DataFrame([human_result,
                       dog_result,
                       tara_result], index=['Human gut', 'Dog gut', 'Tara'],
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
    result_Metabat2_african = get_result('PRJNA504891', 'Metabat2')
    result_SemiBin_african = get_result('PRJNA504891', 'SemiBin')
    result_SemiBin_pretrain_same_african = get_result('PRJNA504891', 'SemiBin_pretrain_same')
    result_SemiBin_pretrain_out_african = get_result('PRJNA504891', 'SemiBin_pretrain_out')

    result_Metabat2_german = get_result('PRJNA290729', 'Metabat2')
    result_SemiBin_german = get_result('PRJNA290729', 'SemiBin')
    result_SemiBin_pretrain_same_german = get_result('PRJNA290729', 'SemiBin_pretrain_same')
    result_SemiBin_pretrain_out_german = get_result('PRJNA290729', 'SemiBin_pretrain_out')

    hq_Metabat2_african = []
    hq_SemiBin_african = []
    hq_SemiBin_pretrain_same_african = []
    hq_SemiBin_pretrain_out_african = []

    hq_Metabat2_german = []
    hq_SemiBin_german  = []
    hq_SemiBin_pretrain_same_german = []
    hq_SemiBin_pretrain_out_german = []

    for sample in PRJNA504891_list:
        hq_Metabat2_african.append(len(result_Metabat2_african[sample]['high quality']))
        hq_SemiBin_african.append(len(result_SemiBin_african[sample]['high quality']))
        hq_SemiBin_pretrain_same_african.append(len(result_SemiBin_pretrain_same_african[sample]['high quality']))
        hq_SemiBin_pretrain_out_african.append(len(result_SemiBin_pretrain_out_african[sample]['high quality']))

    for sample in PRJNA290729_list:
        hq_Metabat2_german.append(len(result_Metabat2_german[sample]['high quality']))
        hq_SemiBin_german.append(len(result_SemiBin_german[sample]['high quality']))
        hq_SemiBin_pretrain_same_german.append(len(result_SemiBin_pretrain_same_german[sample]['high quality']))
        hq_SemiBin_pretrain_out_german.append(len(result_SemiBin_pretrain_out_german[sample]['high quality']))

    print('African:', wilcoxon(hq_Metabat2_african, hq_SemiBin_pretrain_out_african))
    print(len(hq_Metabat2_african))
    print((np.sum(hq_SemiBin_pretrain_out_african) - np.sum(hq_Metabat2_african)) / np.sum(hq_Metabat2_african))
    print(np.sum(hq_Metabat2_african))
    print(np.sum(hq_SemiBin_african))
    print(np.sum(hq_SemiBin_pretrain_same_african))
    print(np.sum(hq_SemiBin_pretrain_out_african))

    print('German:', wilcoxon(hq_Metabat2_german, hq_SemiBin_pretrain_out_german))
    print(len(hq_Metabat2_german))
    print((np.sum(hq_SemiBin_pretrain_out_german) - np.sum(hq_Metabat2_german)) / np.sum(hq_Metabat2_german))
    print(np.sum(hq_Metabat2_german))
    print(np.sum(hq_SemiBin_german))
    print(np.sum(hq_SemiBin_pretrain_same_german))
    print(np.sum(hq_SemiBin_pretrain_out_german))


    subset = pd.DataFrame(np.array([[np.sum(hq_Metabat2_african),np.sum(hq_SemiBin_african), np.sum(hq_SemiBin_pretrain_same_african), np.sum(hq_SemiBin_pretrain_out_african)]]),columns = ['Metabat2','SemiBin','SemiBin(pre-train; internal)','SemiBin(pre-train; external)'], index=['African human gut'])
    ax = subset.plot(kind='bar',width = 0.3,color=['#e7298a', '#7570b3', '#d95f02', '#1b9e77'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['African human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    #ax.legend(loc='lower left', fontsize=6)
    plt.show()
    plt.savefig('holdout_African_bar.pdf', dpi=300, bbox_inches='tight')
    print(subset)
    subset = pd.DataFrame(np.array([[np.sum(hq_Metabat2_german),np.sum(hq_SemiBin_german), np.sum(hq_SemiBin_pretrain_same_german), np.sum(hq_SemiBin_pretrain_out_german)]]),columns = ['Metabat2','SemiBin','SemiBin(pre-train; internal)','SemiBin(pre-train; external)'], index=['German human gut'])
    ax = subset.plot(kind='bar',width = 0.3,color=['#e7298a', '#7570b3', '#d95f02', '#1b9e77'])
    ax.set_yticks(ticks=[0,100,200,300,400,500,600])
    ax.set_yticklabels(labels=[0,100,200,300,400,500,600],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['German human gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    #ax.legend(loc='lower left', fontsize=6)
    plt.show()
    plt.savefig('holdout_German_bar.pdf', dpi=300, bbox_inches='tight')
    print(subset)


def plot_extra_per_sample():
    result_Metabat2_african = get_result('PRJNA504891', 'Metabat2')
    result_SemiBin_pretrain_out_african = get_result('PRJNA504891', 'SemiBin_pretrain_out')

    result_Metabat2_german = get_result('PRJNA290729', 'Metabat2')
    result_SemiBin_pretrain_out_german = get_result('PRJNA290729', 'SemiBin_pretrain_out')

    hq_Metabat2_african = []
    hq_SemiBin_pretrain_out_african = []

    hq_Metabat2_german = []
    hq_SemiBin_pretrain_out_german = []

    for sample in PRJNA504891_list:
        hq_Metabat2_african.append(len(result_Metabat2_african[sample]['high quality']))
        hq_SemiBin_pretrain_out_african.append(len(result_SemiBin_pretrain_out_african[sample]['high quality']))

    for sample in PRJNA290729_list:
        hq_Metabat2_german.append(len(result_Metabat2_german[sample]['high quality']))
        hq_SemiBin_pretrain_out_german.append(len(result_SemiBin_pretrain_out_german[sample]['high quality']))

    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormap = ['#084594','#2171b5','#4292c6','#6baed6','#9ecae1','#c6dbef','#deebf7']
    newcmp = LinearSegmentedColormap.from_list('cmps', colormap)
    data_african = np.zeros(shape=(len(hq_Metabat2_african), 3))
    data_german = np.zeros(shape=(len(hq_Metabat2_german), 3))

    per_sample_result_african = {}
    per_sample_result_german = {}

    for i in range(len(hq_Metabat2_african)):
        if (hq_SemiBin_pretrain_out_african[i], hq_Metabat2_african[i]) not in per_sample_result_african:
            per_sample_result_african[(hq_SemiBin_pretrain_out_african[i], hq_Metabat2_african[i])] = 1
        else:
            per_sample_result_african[(hq_SemiBin_pretrain_out_african[i], hq_Metabat2_african[i])] += 1

    for i, temp in enumerate(per_sample_result_african):
        data_african[i][0] = temp[0]
        data_african[i][1] = temp[1]
        data_african[i][2] = per_sample_result_african[temp]
    data_african = pd.DataFrame(data_african).astype(int)
    data_african.columns = ['x', 'y', 'num']
    print(data_african)

    for i in range(len(hq_Metabat2_german)):
        if (hq_SemiBin_pretrain_out_german[i], hq_Metabat2_german[i]) not in per_sample_result_german:
            per_sample_result_german[(hq_SemiBin_pretrain_out_german[i], hq_Metabat2_german[i])] = 1
        else:
            per_sample_result_german[(hq_SemiBin_pretrain_out_german[i], hq_Metabat2_german[i])] += 1

    for i, temp in enumerate(per_sample_result_german):
        data_german[i][0] = temp[0]
        data_german[i][1] = temp[1]
        data_german[i][2] = per_sample_result_german[temp]
    data_german = pd.DataFrame(data_german).astype(int)
    data_german.columns = ['x', 'y', 'num']
    print(data_german)

    plt.scatter(data_african.x, data_african.y,
                c=data_african.num, s=(data_african.num ** 2) * 60, cmap=newcmp)

    plt.colorbar(shrink=0.5)
    plt.xlabel("SemiBin(pretrain; external)")
    plt.ylabel("Metabat2")

    plt.plot([0, 35], [0, 35], c='black')
    plt.title('African human gut', fontsize=15)
    plt.show()
    plt.savefig('holdout_african_persample.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    plt.scatter(data_german.x, data_german.y,
                c=data_german.num, s=(data_german.num ** 2) * 60, cmap=newcmp)

    plt.colorbar(shrink=0.5)
    plt.xlabel("SemiBin(pretrain; external)")
    plt.ylabel("Metabat2")

    plt.plot([0, 15], [0, 15], c='black')
    plt.title('German human gut', fontsize=15)
    plt.show()
    plt.savefig('holdout_german_persample.pdf', dpi=300, bbox_inches='tight')
    plt.close()


def CAT_mmseqs():
    result_human_single_mmseqs2 = get_num_high_quality('human', method='SemiBin',)
    result_human_single_CAT = get_num_high_quality('human', method='SemiBin_CAT')

    result_dog_single_mmseqs2 = get_num_high_quality('dog', method='SemiBin')
    result_dog_single_CAT = get_num_high_quality('dog', method='SemiBin_CAT')

    result_human_multi_mmseqs2 = get_num_high_quality('human', method='SemiBin', binning_mode='multi_sample')
    result_human_multi_CAT = get_num_high_quality('human', method='SemiBin_CAT', binning_mode='multi_sample')

    result_dog_multi_mmseqs2 = get_num_high_quality('dog', method='SemiBin', binning_mode='multi_sample')
    result_dog_multi_CAT = get_num_high_quality('dog', method='SemiBin_CAT', binning_mode='multi_sample')

    subset = pd.DataFrame(np.array([[result_human_single_CAT,result_human_single_mmseqs2],[result_dog_single_CAT,result_dog_single_mmseqs2]]),columns = ['CAT','mmseqs'], index=['Human gut','Dog gut'])
    print(subset)

    ax = subset.plot(kind='bar',width = 0.6,color = ['#7570b3', '#1b9e77'],figsize=(3,4))
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut','Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.set_title('Single-sample binning', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_CAT_mmseqs_single.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    subset = pd.DataFrame(np.array([[result_human_multi_CAT,result_human_multi_mmseqs2],[result_dog_multi_CAT,result_dog_multi_mmseqs2]]),columns = ['CAT','mmseqs'], index=['Human gut','Dog gut'])
    print(subset)

    ax = subset.plot(kind='bar',width = 0.6,color = ['#7570b3', '#1b9e77'],figsize=(3,4))
    ax.set_yticks(ticks=[0,800,1600,2400,3200])
    ax.set_yticklabels(labels=[0,800,1600,2400,3200],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut','Dog gut'], fontsize=15,color = 'black',rotation = 360)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.set_title('Multi-sample binning', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_CAT_mmseqs_multi.pdf', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # print('human')
    plot_bar_per_sample_com('human',[-20,-10,0,10],[0,300,600,900,1200,1500])
    # print('dog')
    plot_bar_per_sample_com('dog',[-30,-20,-10,0,10],[0,500,1000,1500,2000,2500,3000,3500])
    # print('tara')
    plot_bar_per_sample_com('tara',[-20,-15,-10,-5,0,5], [0,100,200,300,400,500])

    tranfer_multi()

    ### bar plot high quality genomes comparison
    plot_checkm_high_quality_comparison()

    ### venn plot multi annotation comparison
    plot_multi_venn_comparison()

    #alternative_real_compare_test()
    plot_sankey_overlap(dataset='human',output='human_sankey.pdf')
    plot_sankey_overlap(output='dog_sankey.pdf')
    plot_sankey_overlap(dataset='tara',output='tara_sankey.pdf')

    ### recall, precision, F1-score box plot

    plot_overlap_F1('human')
    plot_overlap_F1()
    plot_overlap_F1('tara')

    # ### bar plot the overlap of annotation in all taxi
    plot_all_taxi_overlap(output='dog_taxi_overlap.pdf',y_label=[0,20,40,60,80,100])
    plot_all_taxi_overlap(dataset='human',output='human_taxi_overlap.pdf',y_label=[0,100,200,300,400])
    plot_all_taxi_overlap(dataset='tara',output='tara_taxi_overlap.pdf',y_label=[0,50,100,150,200,250])

    # ## comparison of known and unknown species
    plot_comparison_known_unknown(y_label=[0,500,1000,1500,2000,2500], output='dog_taxi_known_unknown.pdf')
    plot_comparison_known_unknown(dataset='human', y_label=[0,300,600,900,1200,1500], output='human_taxi_known_unknown.pdf')
    plot_comparison_known_unknown(dataset='tara', y_label=[0,100,200,300,400], output='tara_taxi_known_unknown.pdf')

    plot_transfer()

    plot_extra_bar()

    plot_extra_per_sample()

    CAT_mmseqs()

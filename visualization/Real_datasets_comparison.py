"""
This script is used to reproduct the plot of the real datasets
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib_venn import venn3_unweighted
from scipy.stats import wilcoxon

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

contamination = 0.05
dog_list = ['SAMN06172456', 'SAMN06172425', 'SAMN06172487', 'SAMN06172450', 'SAMN06172459', 'SAMN06172479', 'SAMN06172435', 'SAMN06172414', 'SAMN06172409', 'SAMEA103957796', 'SAMN06172442', 'SAMN06172500', 'SAMN06172437', 'SAMN06172413', 'SAMN06172514', 'SAMN06172403', 'SAMN06172471', 'SAMN06172490', 'SAMN06172448', 'SAMN06172504', 'SAMN06172457', 'SAMN06172441', 'SAMN06172422', 'SAMN06172408', 'SAMN06172429', 'SAMN06172420', 'SAMN06172503', 'SAMN06172410', 'SAMN06172458', 'SAMN06172493', 'SAMEA103957794', 'SAMN06172402', 'SAMN06172515', 'SAMN06172462', 'SAMN06172421', 'SAMN06172411', 'SAMN06172511', 'SAMN06172516', 'SAMN06172465', 'SAMN06172419', 'SAMN06172517', 'SAMN06172510', 'SAMN06172418', 'SAMN06172424', 'SAMN06172427', 'SAMN06172453', 'SAMN06172491', 'SAMN06172496', 'SAMN06172513', 'SAMN06172461', 'SAMN06172449', 'SAMN06172426', 'SAMN06172452', 'SAMN06172522', 'SAMN06172400', 'SAMN06172405', 'SAMN06172521', 'SAMN06172407', 'SAMN06172455', 'SAMN06172446', 'SAMN06172467', 'SAMN06172499', 'SAMN06172474', 'SAMN06172412', 'SAMN06172468', 'SAMN06172478', 'SAMN06172423', 'SAMN06172447', 'SAMN06172415', 'SAMN06172523', 'SAMN06172417', 'SAMN06172497', 'SAMN06172498', 'SAMN06172489', 'SAMN06172436', 'SAMN06172432', 'SAMN06172406', 'SAMN06172488', 'SAMN06172502', 'SAMN06172401', 'SAMN06172434', 'SAMN06172416', 'SAMN06172445', 'SAMN06172431', 'SAMN06172438', 'SAMN06172473', 'SAMN06172486', 'SAMN06172472', 'SAMN06172428', 'SAMEA103957793', 'SAMEA103957795', 'SAMN06172443', 'SAMN06172475', 'SAMN06172520', 'SAMN06172495', 'SAMN06172440', 'SAMN06172430', 'SAMN06172481', 'SAMN06172524', 'SAMN06172519', 'SAMN06172454', 'SAMN06172404', 'SAMN06172460', 'SAMN06172433', 'SAMN06172469', 'SAMN06172451', 'SAMN06172476', 'SAMN06172492', 'SAMN06172484', 'SAMN06172509', 'SAMN06172506', 'SAMN06172518', 'SAMN06172477', 'SAMN06172470', 'SAMN06172482', 'SAMN06172512', 'SAMN06172494', 'SAMN06172485', 'SAMN06172508', 'SAMN06172466', 'SAMN06172507', 'SAMN06172444', 'SAMN06172505', 'SAMN06172464', 'SAMN06172439', 'SAMN06172501', 'SAMN06172483', 'SAMN06172463', 'SAMN06172480']

human_list = ['CCMD41521570ST', 'CCMD75147712ST', 'CCMD18579000ST', 'CCMD53508245ST', 'CCMD19168690ST', 'CCMD52117727ST', 'CCMD42956136ST', 'CCMD79349503ST', 'CCMD89306485ST', 'CCMD76409700ST', 'CCMD31134579ST', 'CCMD71242853ST', 'CCMD89107682ST', 'CCMD76222476ST', 'CCMD10032470ST', 'CCMD17410933ST', 'CCMD38158721ST', 'CCMD35081859ST', 'CCMD54057834ST', 'CCMD28738636ST', 'CCMD98702133ST', 'CCMD30626189ST', 'CCMD32965613ST', 'CCMD53522274ST', 'CCMD37575804ST', 'CCMD68973846ST', 'CCMD25475945ST', 'CCMD65406197ST', 'CCMD21703880ST', 'CCMD50300306ST', 'CCMD51228890ST', 'CCMD59540613ST', 'CCMD49942357ST', 'CCMD95431029ST', 'CCMD41202658ST', 'CCMD15562448ST', 'CCMD21593359ST', 'CCMD92404903ST', 'CCMD50538120ST', 'CCMD49461418ST', 'CCMD72690923ST', 'CCMD85481373ST', 'CCMD39882286ST', 'CCMD18829815ST', 'CCMD51154251ST', 'CCMD85661207ST', 'CCMD71915439ST', 'CCMD39157124ST', 'CCMD22852639ST', 'CCMD35801800ST', 'CCMD27463710ST', 'CCMD59583015ST', 'CCMD89967135ST', 'CCMD52145360ST', 'CCMD95676152ST', 'CCMD45004878ST', 'CCMD67373733ST', 'CCMD99929634ST', 'CCMD89643949ST', 'CCMD26625622ST', 'CCMD23541216ST', 'CCMD31009081ST', 'CCMD99440714ST', 'CCMD66848156ST', 'CCMD65222621ST', 'CCMD98531134ST', 'CCMD45812507ST', 'CCMD46727384ST', 'CCMD73128545ST', 'CCMD30627121ST', 'CCMD50529145ST', 'CCMD98198513ST', 'CCMD93755960ST', 'CCMD35633353ST', 'CCMD56948710ST', 'CCMD27867141ST', 'CCMD32288175ST', 'CCMD29706695ST', 'CCMD72666896ST', 'CCMD10191450ST', 'CCMD49025643ST', 'CCMD74592084ST']

tara_list = ['TARA_041_SRF_0.1-0.22', 'TARA_038_SRF_0.22-1.6', 'TARA_076_SRF_0.22-3', 'TARA_023_SRF_0.22-1.6', 'TARA_042_SRF_0.22-1.6', 'TARA_124_SRF_0.22-3', 'TARA_124_SRF_0.22-0.45', 'TARA_066_SRF_lt-0.22', 'TARA_057_SRF_0.22-3', 'TARA_124_SRF_0.45-0.8', 'TARA_004_SRF_0.22-1.6', 'TARA_018_SRF_0.22-1.6', 'TARA_070_SRF_0.22-0.45', 'TARA_034_SRF_lt-0.22', 'TARA_064_SRF_0.22-3', 'TARA_125_SRF_0.22-0.45', 'TARA_111_SRF_0.22-3', 'TARA_122_SRF_0.22-0.45', 'TARA_145_SRF_0.22-3', 'TARA_099_SRF_0.22-3', 'TARA_038_SRF_lt-0.22', 'TARA_082_SRF_0.22-3', 'TARA_041_SRF_lt-0.22', 'TARA_146_SRF_0.22-3', 'TARA_151_SRF_0.22-3', 'TARA_123_SRF_0.22-3', 'TARA_110_SRF_0.22-3', 'TARA_150_SRF_0.22-3', 'TARA_072_SRF_lt-0.22', 'TARA_085_SRF_0.22-3', 'TARA_098_SRF_0.22-3', 'TARA_078_SRF_0.22-3', 'TARA_149_SRF_0.22-3', 'TARA_094_SRF_0.22-3', 'TARA_068_SRF_0.22-3', 'TARA_148_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_018_SRF_lt-0.22', 'TARA_138_SRF_0.22-3', 'TARA_093_SRF_0.22-3', 'TARA_041_SRF_0.22-1.6', 'TARA_122_SRF_0.22-3', 'TARA_078_SRF_0.45-0.8', 'TARA_070_SRF_lt-0.22', 'TARA_065_SRF_lt-0.22', 'TARA_122_SRF_0.1-0.22', 'TARA_036_SRF_0.22-1.6', 'TARA_031_SRF_0.22-1.6', 'TARA_142_SRF_0.22-3', 'TARA_124_SRF_0.1-0.22', 'TARA_036_SRF_lt-0.22', 'TARA_065_SRF_0.22-3', 'TARA_067_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_109_SRF_0.22-3', 'TARA_068_SRF_lt-0.22', 'TARA_109_SRF_lt-0.22', 'TARA_064_SRF_lt-0.22', 'TARA_048_SRF_0.22-1.6', 'TARA_034_SRF_0.22-1.6', 'TARA_070_SRF_0.45-0.8', 'TARA_025_SRF_lt-0.22', 'TARA_133_SRF_0.22-3', 'TARA_096_SRF_0.22-3', 'TARA_038_SRF_0.1-0.22', 'TARA_007_SRF_0.22-1.6', 'TARA_048_SRF_0.1-0.22', 'TARA_140_SRF_0.22-3', 'TARA_034_SRF_0.1-0.22', 'TARA_067_SRF_0.22-3', 'TARA_125_SRF_0.45-0.8', 'TARA_030_SRF_0.22-1.6', 'TARA_031_SRF_lt-0.22', 'TARA_032_SRF_0.22-1.6', 'TARA_070_SRF_0.22-3', 'TARA_132_SRF_0.22-3', 'TARA_076_SRF_lt-0.22', 'TARA_125_SRF_0.1-0.22', 'TARA_123_SRF_0.45-0.8', 'TARA_078_SRF_lt-0.22', 'TARA_068_SRF_0.45-0.8', 'TARA_068_SRF_0.22-0.45', 'TARA_067_SRF_0.22-0.45', 'TARA_100_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_137_SRF_0.22-3', 'TARA_076_SRF_0.22-0.45', 'TARA_125_SRF_0.22-3', 'TARA_078_SRF_0.22-0.45', 'TARA_076_SRF_0.45-0.8', 'TARA_084_SRF_0.22-3', 'TARA_032_SRF_lt-0.22', 'TARA_025_SRF_0.22-1.6', 'TARA_062_SRF_0.22-3', 'TARA_066_SRF_0.22-3', 'TARA_036_SRF_0.1-0.22', 'TARA_056_SRF_0.22-3', 'TARA_072_SRF_0.22-3', 'TARA_128_SRF_0.22-3', 'TARA_052_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_123_SRF_0.22-0.45', 'TARA_102_SRF_0.22-3', 'TARA_065_SRF_0.1-0.22', 'TARA_009_SRF_0.22-1.6', 'TARA_141_SRF_0.22-3', 'TARA_045_SRF_0.22-1.6', 'TARA_042_SRF_lt-0.22', 'TARA_152_SRF_0.22-3']

PRJNA504891_list = ['SRR8784379', 'SRR8180449', 'SRR8784372', 'SRR8784373', 'SRR8784385', 'SRR8784383', 'SRR8784375', 'SRR8784384', 'SRR8784360', 'SRR8784395','SRR8784376', 'SRR8784393', 'SRR8784363', 'SRR8784378', 'SRR8784361', 'SRR8784357', 'SRR8784394', 'SRR8784370', 'SRR8784397', 'SRR8784374', 'SRR8784364', 'SRR8784390', 'SRR8784381','SRR8784391', 'SRR8784358', 'SRR8784380', 'SRR8784387', 'SRR8784386', 'SRR8784365', 'SRR8784359', 'SRR8784382', 'SRR8784396','SRR8784377', 'SRR8784362', 'SRR8784389', 'SRR8180448','SRR8180450', 'SRR8784366','SRR8784371', 'SRR8784353', 'SRR8784354', 'SRR8180446', 'SRR8784388', 'SRR8784367', 'SRR8784392', 'SRR8784369', 'SRR8180447', 'SRR8784356','SRR8784355', 'SRR8784368']

PRJNA290729_list = ['SAMN03922475', 'SAMN03922449', 'SAMN03922521', 'SAMN03922488', 'SAMN03922526', 'SAMN03922468', 'SAMN03922500', 'SAMN03922512', 'SAMN03922494', 'SAMN03922450', 'SAMN03922479', 'SAMN03922484', 'SAMN03922492', 'SAMN03922470', 'SAMN03922480', 'SAMN03922505', 'SAMN03922516', 'SAMN03922527', 'SAMN03922513', 'SAMN03922472', 'SAMN03922504', 'SAMN03922523', 'SAMN03922528', 'SAMN03922510', 'SAMN03922497', 'SAMN03922518', 'SAMN03922465', 'SAMN03922517', 'SAMN03922522', 'SAMN03922458', 'SAMN03922511', 'SAMN03922531', 'SAMN03922462', 'SAMN03922457', 'SAMN03922507', 'SAMN03922509', 'SAMN03922474', 'SAMN03922529', 'SAMN03922456', 'SAMN03922506', 'SAMN03922464', 'SAMN03922453', 'SAMN03922463', 'SAMN03922466', 'SAMN03922539', 'SAMN03922503', 'SAMN03922477', 'SAMN03922495', 'SAMN03922451', 'SAMN03922538', 'SAMN03922461', 'SAMN03922532', 'SAMN03922476', 'SAMN03922469', 'SAMN03922540', 'SAMN03922533', 'SAMN03922530', 'SAMN03922536', 'SAMN03922519', 'SAMN03922471', 'SAMN03922489', 'SAMN03922524', 'SAMN03922496', 'SAMN03922467', 'SAMN03922520', 'SAMN03922483', 'SAMN03922452', 'SAMN03922508', 'SAMN03922486', 'SAMN03922473', 'SAMN03922515', 'SAMN03922455', 'SAMN03922534', 'SAMN03922490', 'SAMN03922498', 'SAMN03922525', 'SAMN03922487', 'SAMN03922482', 'SAMN03922501', 'SAMN03922537', 'SAMN03922535', 'SAMN03922485', 'SAMN03922459', 'SAMN03922499', 'SAMN03922514', 'SAMN03922454', 'SAMN03922493', 'SAMN03922478', 'SAMN03922491', 'SAMN03922481', 'SAMN03922502', 'SAMN03922460']

soil_list = ['SAMN06268061', 'SAMN06268063', 'SAMN06264885', 'SAMN06267090',
                 'SAMN06267080', 'SAMN06266457', 'SAMN05421921', 'SAMN06264649',
                 'SAMN06264650', 'SAMN07631258', 'SAMN06264630', 'SAMN06267102',
                 'SAMN06267104', 'SAMN06264384', 'SAMN06266487', 'SAMN06266460',
                 'SAMN06266447', 'SAMN06264634', 'SAMN06266424', 'SAMN06268167',
                 'SAMN06266446', 'SAMN06267099', 'SAMN06266484', 'SAMN06266459',
                 'SAMN06267079', 'SAMN06267094', 'SAMN06264884', 'SAMN06266490',
                 'SAMN06266453', 'SAMN06264385', 'SAMN06264631', 'SAMN06266388',
                 'SAMN06266336', 'SAMN06266479', 'SAMN06266485', 'SAMN06267092',
                 'SAMN07631257', 'SAMN06266454', 'SAMN06266483', 'SAMN06268059',
                 'SAMN06267098', 'SAMN06268058', 'SAMN06268170', 'SAMN06266423',
                 'SAMN06264948', 'SAMN06267083', 'SAMN06264648', 'SAMN07631255',
                 'SAMN06266478', 'SAMN06267085', 'SAMN06266448', 'SAMN06267101',
                 'SAMN06268168', 'SAMN06267097', 'SAMN06266475', 'SAMN06266450',
                 'SAMN06264881', 'SAMN06264635', 'SAMN06267088', 'SAMN06266458',
                 'SAMN06267095', 'SAMN06264383', 'SAMN06266461', 'SAMN06266449',
                 'SAMN06267096', 'SAMN06267087', 'SAMN06267103', 'SAMN06266486',
                 'SAMN06267084', 'SAMN06264882', 'SAMN06266387', 'SAMN06266473',
                 'SAMN05421920', 'SAMN06266491', 'SAMN05421922', 'SAMN06267086',
                 'SAMN06264947', 'SAMN06268062', 'SAMN07631256', 'SAMN06267100',
                 'SAMN06267091', 'SAMN06268166', 'SAMN06264632', 'SAMN06264883',
                 'SAMN06266481', 'SAMN06266482', 'SAMN06267089', 'SAMN06267093',
                 'SAMN06264633', 'SAMN06266456', 'SAMN06267081', 'SAMN06266474',
                 'SAMN06267105', 'SAMN05421524', 'SAMN05421649', 'SAMN06268169',
                 'SAMN06267082', 'SAMN06268060', 'SAMN06266477', 'SAMN06266489',
                 'SAMN06266455']

transfer_single_human = ['CCMD41521570ST', 'CCMD21593359ST', 'CCMD89107682ST', 'CCMD18579000ST', 'CCMD76222476ST', 'CCMD22852639ST', 'CCMD99440714ST', 'CCMD38158721ST', 'CCMD50300306ST', 'CCMD98198513ST']

transfer_single_dog = ['SAMN06172505', 'SAMN06172402', 'SAMN06172415', 'SAMN06172428', 'SAMN06172503', 'SAMN06172429', 'SAMN06172448', 'SAMN06172456', 'SAMN06172449', 'SAMN06172418']

transfer_single_tara = ['TARA_125_SRF_0.22-3', 'TARA_102_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_122_SRF_0.1-0.22', 'TARA_066_SRF_lt-0.22', 'TARA_067_SRF_0.22-3', 'TARA_123_SRF_0.45-0.8', 'TARA_093_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_038_SRF_0.22-1.6']

transfer_multi_human = ['CCMD54057834ST', 'CCMD99440714ST', 'CCMD10191450ST', 'CCMD92404903ST', 'CCMD53522274ST', 'CCMD50300306ST', 'CCMD89107682ST', 'CCMD95431029ST', 'CCMD27463710ST', 'CCMD49025643ST']

transfer_multi_dog = ['SAMN06172506', 'SAMN06172436', 'SAMN06172403', 'SAMN06172520', 'SAMN06172512', 'SAMN06172433', 'SAMN06172493', 'SAMN06172422', 'SAMN06172523', 'SAMN06172407']

transfer_multi_tara = ['TARA_045_SRF_0.22-1.6', 'TARA_041_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_140_SRF_0.22-3', 'TARA_138_SRF_0.22-3', 'TARA_034_SRF_lt-0.22', 'TARA_066_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_056_SRF_0.22-3', 'TARA_034_SRF_0.22-1.6']

transfer_multi_soil = ['SAMN06267105', 'SAMN06267087', 'SAMN06267090', 'SAMN06268058', 'SAMN06266450', 'SAMN06266387', 'SAMN06267104', 'SAMN06264635', 'SAMN06264650', 'SAMN06268059']

def get_result(dataset='dog', method='Maxbin2', binning_mode = 'single_sample', checkm_only = False):
    """
    dataset: dog, human, gut
    method: Maxbin2, Metabat2, VAMB, S3N2Bin
    binning_mode: single_sample, multi_sample
    checkm_only: if just using checkm or using checkm and GUNC
    """
    if dataset == 'dog':
        sample_list = dog_list
    elif dataset == 'soil':
        sample_list = soil_list
    elif dataset == 'human':
        sample_list = human_list
    elif dataset == 'tara':
        sample_list = tara_list
    elif dataset == 'soil':
        sample_list = soil_list
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
    elif dataset == 'soil_transfer_multi':
        sample_list = transfer_multi_soil
        dataset = 'soil'
    else:
        raise KeyError(f"Unknown dataset {dataset}")

    result = {}
    if method == 'VAMB' and binning_mode == 'multi_sample':
        result = {'high quality': []}
        if dataset != 'soil':
            binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/multi_sample/VAMB_multi.csv'.format(dataset),index_col=0)
        else:
            binning_result = pd.read_csv('updated_results/Soil_benchmark/CheckM_GUNC/multi_sample/VAMB_multi.csv', index_col=0)
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
            if dataset != 'soil' and dataset != 'soil_transfer_multi':
                binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/{1}/{2}/{3}/result.csv'.format(dataset,binning_mode, sample, method), index_col=0)
            else:
                binning_result = pd.read_csv('updated_results/Soil_benchmark/CheckM_GUNC/{0}/{1}/{2}/result.csv'.format(binning_mode, sample, method), index_col=0)
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
    if dataset == 'soil':
        sample_list = transfer_multi_soil

    num_high_quality = 0
    for sample in sample_list:
        if dataset != 'soil':
            binning_result = pd.read_csv('Results/Real/CheckM_GUNC/Pretrain/{3}/{0}/{1}/{2}/result.csv'.format(num_run, num_sample, sample,dataset),index_col=0)
        else:
            binning_result = pd.read_csv('updated_results/Soil_round/{0}/{1}/{2}/result.csv'.format(num_run, num_sample, sample))
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
    soil_result = {1:[],3:[],5:[],10:[],15:[],20:[]}

    # soil_result = {1:[9, 12, 10, 9, 11],
    #                3:[13, 13, 13, 10, 11],
    #                5:[15, 13, 15, 15, 14],
    #                10:[12, 15, 15, 13, 12],
    #                15:[15, 17, 14, 15, 15],
    #                20:[13, 13, 14, 12, 14]}

    for j in [1,3,5,10,15,20]:
        for i in [1, 2, 3, 4, 5]:
            human_result[j].append(get_num_high_quality_pretrain('human', i, j))
            dog_result[j].append(get_num_high_quality_pretrain('dog', i, j))
            tara_result[j].append(get_num_high_quality_pretrain('tara', i, j))
            soil_result[j].append(get_num_high_quality_pretrain('soil', i, j))
    print(soil_result)
    human_SemiBin = get_num_high_quality('human_transfer_multi', 'SemiBin')
    human_Metabat2 = get_num_high_quality('human_transfer_multi', 'Metabat2')

    dog_SemiBin = get_num_high_quality('dog_transfer_multi', 'SemiBin')
    dog_Metabat2 = get_num_high_quality('dog_transfer_multi', 'Metabat2')

    tara_SemiBin = get_num_high_quality('tara_transfer_multi', 'SemiBin')
    tara_Metabat2 = get_num_high_quality('tara_transfer_multi', 'Metabat2')

    soil_SemiBin = get_num_high_quality('soil_transfer_multi', 'SemiBin')
    soil_Metabat2 = get_num_high_quality('soil_transfer_multi', 'Metabat2')
    print(soil_SemiBin)
    print(soil_Metabat2)

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

    y_soil_multi = []
    y_soil_err_max = []
    y_soil_err_min = []
    y_soil_semibin = [soil_SemiBin] * 6
    y_soil_metabat2 = [soil_Metabat2] * 6

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

        y_soil_multi.append(np.mean(soil_result[i]))
        y_soil_err_max.append(np.max(soil_result[i]) - np.mean(soil_result[i]))
        y_soil_err_min.append(np.mean(soil_result[i]) - np.min(soil_result[i]))
    fig, ax = plt.subplots(figsize = (4,2))
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    line_width = 1
    plt.errorbar(x1,y_human_multi,yerr=[y_human_err_min, y_human_err_max],label='SemiBin(pretrain)',color = '#005A32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_human_semibin, label='SemiBin',color='#41AB5D',linewidth=line_width)
    plt.plot(x1, y_human_metabat2, label='Metabat2',color='#ec7014',linewidth=line_width)
    ax.set_xticks([1,2,3,4,5,6])
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_xlabel('Number of samples', fontsize=15, color='black')
    ax.set_xticklabels([1,3,5,10,15,20],color = 'black')
    plt.yticks([0,25,50,75,100,125,150,175,200], color='black')
    plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize=8)
    plt.savefig('transfer_multiple_human.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize = (4,2))
    plt.errorbar(x1,y_dog_multi,yerr=[y_dog_err_min, y_dog_err_max],label='SemiBin(pretrain)',color='#005A32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_dog_semibin, label='SemiBin',color='#41AB5D',linewidth=line_width)
    plt.plot(x1, y_dog_metabat2, label='Metabat2',color='#ec7014',linewidth=line_width)
    ax.set_xticks([1,2,3,4,5,6])
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_xlabel('Number of samples', fontsize=15, color='black')
    ax.set_xticklabels([1,3,5,10,15,20],color = 'black')
    plt.yticks([0,25,50,75,100,125,150,175,200], color='black')
    plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize=8)
    plt.savefig('transfer_multiple_dog.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize = (4,2))
    plt.errorbar(x1,y_tara_multi,yerr=[y_tara_err_min, y_tara_err_max],label='SemiBin(pretrain)',color='#005A32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_tara_semibin, label='SemiBin',color='#41AB5D',linewidth=line_width)
    plt.plot(x1, y_tara_metabat2, label='Metabat2',color='#ec7014',linewidth=line_width)
    ax.set_xticks([1,2,3,4,5,6])
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_xlabel('Number of samples', fontsize=15, color='black')
    ax.set_xticklabels([1,3,5,10,15,20],color = 'black')
    plt.yticks([0,5,10,15,20,25], color='black')
    plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize=8)
    plt.savefig('transfer_multiple_ocean.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize = (4,2))
    plt.errorbar(x1,y_soil_multi,yerr=[y_soil_err_min, y_soil_err_max],label='SemiBin(pretrain)',color='#005A32',fmt='k-o',lw = line_width,elinewidth=1,ms=3,capsize=3)
    plt.plot(x1, y_soil_semibin, label='SemiBin',color='#41AB5D',linewidth=line_width)
    plt.plot(x1, y_soil_metabat2, label='Metabat2',color='#ec7014',linewidth=line_width)
    ax.set_xticks([1,2,3,4,5,6])
    ax.set_ylabel('High quality bins', fontsize=15, color='black')
    ax.set_xlabel('Number of samples', fontsize=15, color='black')
    ax.set_xticklabels([1,3,5,10,15,20],color = 'black')
    plt.yticks([0,5,10,15,20], color='black')
    plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0., fontsize=8)
    plt.savefig('transfer_multiple_soil.pdf', dpi=300, bbox_inches='tight')
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

    # diff_single = counts_single.T - counts_single.SemiBin_pretrain
    # diff_multi = counts_multi.T - counts_multi.SemiBin_multi

    diff_single = counts_single.SemiBin_pretrain - counts_single.T
    print('Others better than SemiBin:')
    print((diff_single.T['Maxbin2'] < 0).sum() + (diff_single.T['Metabat2'] < 0).sum() + (diff_single.T['VAMB'] < 0).sum())
    diff_multi = counts_multi.SemiBin_multi - counts_multi.T
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
    print(v)
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

    num_soil_maxbin2_single = get_num_high_quality(dataset='soil',checkm_only=True)
    num_soil_vamb_single = get_num_high_quality(dataset='soil', method='VAMB',checkm_only=True)
    num_soil_metabat2_single = get_num_high_quality(dataset='soil', method='Metabat2',checkm_only=True)
    num_soil_semibin_single = get_num_high_quality(dataset='soil', method='SemiBin',checkm_only=True)
    num_soil_semibin_pretrain_single = get_num_high_quality(dataset='soil', method='SemiBin_pretrain',checkm_only=True)

    print(num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single, num_dog_semibin_pretrain_single)
    print(num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single)
    print(num_tara_maxbin2_single,num_tara_vamb_single,num_tara_metabat2_single,num_tara_semibin_single,num_tara_semibin_pretrain_single)
    print(num_soil_maxbin2_single, num_soil_vamb_single, num_soil_metabat2_single, num_soil_semibin_single, num_soil_semibin_pretrain_single)

    subset = pd.DataFrame(np.array([[num_dog_maxbin2_single,num_dog_vamb_single,num_dog_metabat2_single,num_dog_semibin_single,num_dog_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pretrain)'], index=['Dog gut'])
    print(subset)
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000])
    ax.set_yticklabels(labels=[0,500,1000,1500,2000,2500,3000],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Dog gut'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_dog_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_human_maxbin2_single,num_human_vamb_single,num_human_metabat2_single,num_human_semibin_single, num_human_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pretrain)'], index=['Human gut'])
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
    print(subset)
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,100,200,300,400,500])
    ax.set_yticklabels(labels=[0,100,200,300,400,500],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_tara_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_soil_maxbin2_single,num_soil_vamb_single,num_soil_metabat2_single,num_soil_semibin_single, num_soil_semibin_pretrain_single]]),columns = ['Maxbin2','VAMB','Metabat2','SemiBin','SemiBin(pre-train)'], index=['Soil'])
    print(subset)
    ax = subset.plot(kind='bar',figsize=(3,4),legend = False, color=['#e7298a','#7570b3','#ec7014','#41AB5D','#005A32'])
    ax.set_yticks(ticks=[0,20,40,60,80,100])
    ax.set_yticklabels(labels=[0,20,40,60,80,100],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Soil'], fontsize=15,color = 'black',rotation = 360)
    plt.savefig('Real_soil_hq_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    num_dog_vamb_mulit = get_num_high_quality(method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_dog_semibin_multi = get_num_high_quality(method='SemiBin',binning_mode='multi_sample',checkm_only=True)

    num_human_vamb_multi = get_num_high_quality(dataset='human', method='VAMB',binning_mode='multi_sample',checkm_only=True)
    num_human_semibin_multi = get_num_high_quality(dataset='human', method='SemiBin',binning_mode='multi_sample',checkm_only=True)

    num_tara_vamb_multi = get_num_high_quality(dataset='tara', method='VAMB', binning_mode='multi_sample',checkm_only=True)
    num_tara_semibin_multi = get_num_high_quality(dataset='tara', method='SemiBin', binning_mode='multi_sample',checkm_only=True)

    num_soil_vamb_multi = get_num_high_quality(dataset='soil', method='VAMB', binning_mode='multi_sample',checkm_only=True)
    num_soil_semibin_multi = get_num_high_quality(dataset='soil', method='SemiBin', binning_mode='multi_sample',checkm_only=True)

    print(num_dog_vamb_mulit,num_dog_semibin_multi)
    print(num_human_vamb_multi, num_human_semibin_multi)
    print(num_tara_vamb_multi,num_tara_semibin_multi)
    print(num_soil_vamb_multi,num_soil_semibin_multi)

    subset = pd.DataFrame(np.array([[num_dog_vamb_mulit,num_dog_semibin_multi]]),columns = ['VAMB','SemiBin'], index=['Dog gut'])
    print(subset)
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
    print(subset)
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
    print(subset)
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#41AB5D','#7570b3'])
    ax.set_yticks(ticks=[0, 150, 300, 450, 600, 750])
    ax.set_yticklabels(labels=[0, 150, 300, 450, 600, 750], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Ocean'], fontsize=15, color='black', rotation=360)
    plt.savefig('Real_tara_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    subset = pd.DataFrame(np.array([[num_soil_vamb_multi,num_soil_semibin_multi]]), columns=['VAMB', 'SemiBin'],
                          index=['Soil'])
    print(subset)
    ax = subset.plot(kind='bar', figsize=(2, 4),legend = False, color=['#41AB5D','#7570b3'])
    ax.set_yticks(ticks=[0, 50, 100, 150, 200, 250])
    ax.set_yticklabels(labels=[0, 50, 100, 150, 200, 250], fontsize=12, color='black')
    ax.set_xticklabels(labels=['Soil'], fontsize=15, color='black', rotation=360)
    plt.savefig('Real_soil_hq_multi_checkm.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    plt.show()

    print('Dog single ', (num_dog_semibin_pretrain_single - num_dog_metabat2_single), (num_dog_semibin_pretrain_single - num_dog_metabat2_single)/num_dog_metabat2_single)
    print('human single ', (num_human_semibin_pretrain_single - num_human_metabat2_single), (num_human_semibin_pretrain_single - num_human_metabat2_single) / num_human_metabat2_single)
    print('ocean single ', (num_tara_semibin_pretrain_single - num_tara_metabat2_single), (num_tara_semibin_pretrain_single - num_tara_metabat2_single) / num_tara_metabat2_single)
    print('soil single ', (num_soil_semibin_pretrain_single - num_soil_metabat2_single), (num_soil_semibin_pretrain_single - num_soil_metabat2_single) / num_soil_metabat2_single)

    print('Dog multi', (num_dog_semibin_multi - num_dog_vamb_mulit), (num_dog_semibin_multi - num_dog_vamb_mulit)/num_dog_vamb_mulit)
    print('human multi', (num_human_semibin_multi - num_human_vamb_multi), (num_human_semibin_multi - num_human_vamb_multi) / num_human_vamb_multi)
    print('tara multi', (num_tara_semibin_multi - num_tara_vamb_multi), (num_tara_semibin_multi - num_tara_vamb_multi) / num_tara_vamb_multi)
    print('soil multi', (num_soil_semibin_multi - num_soil_vamb_multi), (num_soil_semibin_multi - num_soil_vamb_multi) / num_soil_vamb_multi)

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
    # for temp in family_list:
    #     if temp[0:3] != 'f__':
    #         print(temp)
    # for temp in genus_list:
    #     if temp[0:3] != 'g__':
    #         print(temp)
    # for temp in species_list:
    #     if temp[0:3] != 's__':
    #         print(temp)
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

    print('soil')
    soil_VAMB_family, soil_VAMB_genus, soil_VAMB_species = get_taxi_list( 'updated_results/Soil_benchmark/gtdbtk_annotations/multi_sample/VAMB/gtdbtk.bac120.summary.tsv','updated_results/Soil_benchmark/gtdbtk_annotations/multi_sample/VAMB/gtdbtk.ar122.summary.tsv')
    soil_SemiBin_family, soil_SemiBin_genus, soil_SemiBin_species = get_taxi_list('updated_results/Soil_benchmark/gtdbtk_annotations/multi_sample/SemiBin/gtdbtk.bac120.summary.tsv', 'updated_results/Soil_benchmark/gtdbtk_annotations/multi_sample/SemiBin/gtdbtk.ar122.summary.tsv')
    soil_SemiBin_family_single, soil_SemiBin_genus_single, soil_SemiBin_species_single = get_taxi_list('updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin/gtdbtk.bac120.summary.tsv',
'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin/gtdbtk.ar122.summary.tsv')

    out = venn3_unweighted([set(dog_VAMB_family), set(dog_SemiBin_family), set(dog_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Dog gut)", fontsize=20, alpha=1.0, color='black')
    plt.savefig('multi_sample_dog_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(dog_VAMB_genus), set(dog_SemiBin_genus), set(dog_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(dog_VAMB_species), set(dog_SemiBin_species), set(dog_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Dog gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_dog_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(human_VAMB_family), set(human_SemiBin_family), set(human_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(human_VAMB_genus), set(human_SemiBin_genus), set(human_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(human_VAMB_species), set(human_SemiBin_species), set(human_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Human gut)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_human_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(tara_VAMB_family), set(tara_SemiBin_family), set(tara_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(tara_VAMB_genus), set(tara_SemiBin_genus), set(tara_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(tara_VAMB_species), set(tara_SemiBin_species), set(tara_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Ocean)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_tara_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(soil_VAMB_family), set(soil_SemiBin_family), set(soil_SemiBin_family_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("family(Soil)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_soil_family.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(soil_VAMB_genus), set(soil_SemiBin_genus), set(soil_SemiBin_genus_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("genus(Soil)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_soil_genus.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    out = venn3_unweighted([set(soil_VAMB_species), set(soil_SemiBin_species), set(soil_SemiBin_species_single)],
                           set_labels=('VAMB', 'SemiBin_multi', 'SemiBin_single'), normalize_to=5.0,set_colors=('#ec7014', '#7570b3', '#1b9e77'))
    for text in out.subset_labels:
        text.set_fontsize(16)
    for text in out.set_labels:
        text.set_fontsize(16)
    plt.title("species(Soil)", fontsize=20, alpha=1.0, color='black')

    plt.savefig('multi_sample_soil_species.pdf', dpi=300, bbox_inches='tight')
    plt.close()




def get_overlap(dataset = 'dog'):
    if dataset == 'dog':
        sample_list = dog_list
    if dataset == 'human':
        sample_list = human_list
    if dataset == 'tara':
        sample_list = tara_list
    if dataset == 'soil':
        sample_list = soil_list

    SemiBin_hq_list = []
    SemiBin_mq_list = []
    SemiBin_others_list = []

    Metabat2_hq_list = []
    Metabat2_mq_list = []
    Metabat2_others_list = []

    SemiBin_bin_dict = {}
    Metabat2_bin_dict = {}

    for sample in sample_list:
        if dataset != 'soil':
            SemiBin_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/SemiBin_pretrain/result.csv'.format(dataset,sample),index_col=0)
        else:
            SemiBin_binning_result = pd.read_csv('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{0}/SemiBin_pretrain/result.csv'.format(sample), index_col=0)
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

        if dataset != 'soil':
            Metabat2_binning_result = pd.read_csv('Results/Real/CheckM_GUNC/{0}/single_sample/{1}/Metabat2/result.csv'.format(dataset,sample),index_col=0)
        else:
            Metabat2_binning_result = pd.read_csv('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{0}/Metabat2/result.csv'.format(sample), index_col=0)

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

        if dataset != 'soil':
            data = pd.read_csv('Results/Real/Mash_dist/{0}/{1}/Mash_SemiBin_pretrain_Metabat2.txt'.format(dataset, sample),header=None,sep='\t')
        else:
            data = pd.read_csv('updated_results/Soil_benchmark/Mash_dist/{}/Mash_SemiBin_pretrain_Metabat2.txt'.format(sample), header=None,sep='\t')
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
    SemiBin_miss = SemiBin_high_quality - len(SemiBin_hq_list) - len(SemiBin_mq_list) - len(SemiBin_others_list)
    Metabat2_miss = Metabat2_high_quality - len(Metabat2_hq_list) - len(Metabat2_mq_list) - len(Metabat2_others_list)

    import plotly.graph_objects as go

    unique_list = ["HQ-HQ", "HQ-MQ", "HQ-LQ", "HQ-Miss", "HQ-HQ", "HQ-MQ", "HQ-LQ", "HQ-Miss"]
    sources = [0, 0, 0, 0, 4, 4, 4, 4]
    targets = [4, 5, 6, 7, 0, 1, 2, 3]
    values = [len(Metabat2_hq_list), len(Metabat2_mq_list), len(Metabat2_others_list), Metabat2_miss,len(SemiBin_hq_list), len(SemiBin_mq_list), len(SemiBin_others_list), SemiBin_miss, ]

    print(len(Metabat2_hq_list), len(Metabat2_mq_list), len(Metabat2_others_list), Metabat2_miss,len(SemiBin_hq_list), len(SemiBin_mq_list), len(SemiBin_others_list), SemiBin_miss,)

    print(Metabat2_high_quality - len(Metabat2_hq_list), (Metabat2_high_quality - len(Metabat2_hq_list)) / Metabat2_high_quality)

    print(SemiBin_high_quality - len(SemiBin_hq_list), (SemiBin_high_quality - len(SemiBin_hq_list)) / SemiBin_high_quality)

    layout = go.Layout(autosize=True, margin={'l': 0, 'r': 0, 't': 0, 'b': 0})
    # plotly setup
    if dataset == 'human' or dataset == 'tara':
        x_ = [0.1, 0.1, 0.1, 0.1, 0.4, 0.4, 0.4, 0.4]
        y_ = [0.15, 0.2, 0.3, 0.4, 0.1, 0.3, 0.4, 0.48]
    elif dataset == 'dog':
        x_ = [0.1, 0.1, 0.1, 0.1, 0.4, 0.4, 0.4, 0.4]
        y_ = [0.15, 0.2, 0.3, 0.4, 0.1, 0.3, 0.4, 0.52]
    else:
        x_ = [0.1, 0.1, 0.1, 0.1, 0.4, 0.4, 0.4, 0.4]
        y_ = [0.15, 0.2, 0.3, 0.4, 0.05, 0.3, 0.4, 0.52]
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
    print(len(SemiBin_hq_list), len(set(SemiBin_hq_list)))
    print(len(SemiBin_mq_list), len(set(SemiBin_mq_list)))
    print(len(SemiBin_others_list), len(set(SemiBin_others_list)))
    print(len(Metabat2_hq_list), len(set(Metabat2_hq_list)))
    print(len(Metabat2_mq_list), len(set(Metabat2_mq_list)))
    print(len(Metabat2_others_list), len(set(Metabat2_others_list)))

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
    print(data)
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
    if dataset == 'soil':
        ax.set_title('{}'.format('Soil'), fontsize=15, alpha=1.0, color='black')
    plt.title('HQ in both')
    plt.savefig('{}_F1_distribution.pdf'.format(dataset), dpi=300, bbox_inches='tight')
    plt.close()


def get_taxa_sets(bac_path, arr_path = None):
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
        if species[0:3] != 's__' or genus[0:3] != 'g__' or family[0:3] != 'f__' or order[0:3] != 'o__' or Class[0:3] != 'c__' or phylum[0:3] != 'p__' or domain[0:3] != 'd__':
            print('error')

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


def plot_taxa_overlap():
    COLORS = {
            'both': '#7570b3',
            'semibin': '#1b9e77',
            'metabat2': '#ec7014',
            }
    fig,axes = plt.subplots(2,2, sharex=True, figsize=[6,8])
    for dataset,ax in zip(['dog', 'tara', 'human', 'soil'], axes.flat):
        if dataset == 'dog' or dataset == 'tara':
            SemiBin_taxa = get_taxa_sets('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset))

            Metabat2_taxa = get_taxa_sets('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
        elif dataset == 'human':
            SemiBin_taxa = get_taxa_sets('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.ar122.summary.tsv'.format(dataset))

            Metabat2_taxa = get_taxa_sets('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

        elif dataset == 'soil':
            SemiBin_taxa = get_taxa_sets('updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin_pretrain/gtdbtk.bac120.summary.tsv', 'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin_pretrain/gtdbtk.ar122.summary.tsv')
            Metabat2_taxa =  get_taxa_sets('updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/Metabat2/gtdbtk.bac120.summary.tsv','updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/Metabat2/gtdbtk.ar122.summary.tsv')

        assert len(SemiBin_taxa) == len(Metabat2_taxa)
        subset = pd.DataFrame([(len(sb & mb2), len(sb - mb2), len(mb2 - sb)) for sb,mb2 in zip(SemiBin_taxa, Metabat2_taxa)],
                    index=['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
                      columns=['Both', 'SemiBin(pretrain) only', 'Metabat2 only'])
        subset = subset.iloc[::-1]
        ax.bar(np.arange(len(subset))-.15, subset.Both + subset['Metabat2 only'], width=.3, color=COLORS['metabat2'], label='Metabat2 only')
        ax.bar(np.arange(len(subset))+.15, subset.Both + subset['SemiBin(pretrain) only'], width=.3, color=COLORS['semibin'], label='SemiBin(pretrain) only')
        ax.bar(np.arange(len(subset)), subset.Both, width=.6, color=COLORS['both'], label='Both')
        title = {
                'dog': 'Dog gut',
                'human': 'Human gut',
                'tara': 'Ocean',
                'soil': 'Soil',
                }[dataset]
        ax.set_title(title, fontsize=10, alpha=1.0, color='black')
        ax.set_ylabel('Number of taxa', fontsize=10, color='black')
    axes[0,0].legend(loc='upper right', fontsize=10)

    for ax in axes[1]:
        ax.set_xticks(np.arange(len(subset)))
        ax.set_xticklabels(labels=subset.index, rotation=90, minor=False, fontsize=10, color='black')
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig('overlaps.svg')


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
            if species[0:3] != 's__':
                print('error')
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
                if species[0:3] != 's__':
                    print('error')
                known_list.append(species)
            else:
                unknown_list.append(species)

    return known_list, unknown_list


def plot_comparison_known_unknown(dataset = 'dog', y_label = None,output = None):
    if dataset == 'dog' or dataset == 'tara':
        SemiBin_known, SemiBin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset))
    if dataset == 'human':
        SemiBin_known, SemiBin_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/SemiBin_pretrain/gtdbtk.ar122.summary.tsv'.format(dataset))

        Metabat2_known, Metabat2_unknown = get_known_unknown('Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.bac120.summary.tsv'.format(dataset), 'Results/Real/gtdbtk_annotations/single_sample/{0}/Metabat2/gtdbtk.ar122.summary.tsv'.format(dataset))

    if dataset == 'soil':
        SemiBin_known, SemiBin_unknown = get_known_unknown(
            'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin_pretrain/gtdbtk.bac120.summary.tsv',
            'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/SemiBin_pretrain/gtdbtk.ar122.summary.tsv')

        Metabat2_known, Metabat2_unknown = get_known_unknown(
            'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/Metabat2/gtdbtk.bac120.summary.tsv',
            'updated_results/Soil_benchmark/gtdbtk_annotations/single_sample/Metabat2/gtdbtk.ar122.summary.tsv')

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
    if dataset == 'soil':
        ax.set_title('{}'.format('Soil'), fontsize=15, alpha=1.0, color='black')
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

def seq_depth_effect():
    # human gut
    num_human_semibin_pretrain_single = get_results_table(dataset='human', method='SemiBin_pretrain', checkm_only=False)
    print(num_human_semibin_pretrain_single['nr_hq'].sum())
    sample_list = num_human_semibin_pretrain_single['sample'].values

    human_gut_meta = pd.read_csv('updated_results/sequence_depth/human_gut.txt', index_col=0, sep=',')
    sample_name = human_gut_meta['Alias'].values
    sample_base = human_gut_meta['Bases'].values
    base_dict = dict()

    for sample, num_base in zip(sample_name, sample_base):
        sample_ = sample.split('-')[0]
        if sample_ not in base_dict:
            base_dict[sample_] = num_base/1e9
        else:
            base_dict[sample_] += num_base / 1e9

    num_base_list = []

    for sample in sample_list:
        num_base_list.append(base_dict[sample])

    num_human_semibin_pretrain_single['bases'] = num_base_list
    print(num_human_semibin_pretrain_single)
    sns.lmplot(x = 'bases', y = 'nr_hq', data = num_human_semibin_pretrain_single)
    plt.title('Human gut', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('seq_effect_human.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # dog gut
    num_dog_semibin_pretrain_single = get_results_table(dataset='dog', method='SemiBin_pretrain', checkm_only=False)
    print(num_dog_semibin_pretrain_single['nr_hq'].sum())
    sample_list = num_dog_semibin_pretrain_single['sample'].values

    dog_gut_meta = pd.read_csv('updated_results/sequence_depth/dog.txt', index_col=0, sep=',')
    sample_name = dog_gut_meta['BioSample'].values
    sample_base = dog_gut_meta['Bases'].values
    base_dict = dict()

    for sample, num_base in zip(sample_name, sample_base):
        if sample not in base_dict:
            base_dict[sample] = num_base/1e9
        else:
            base_dict[sample] += num_base / 1e9

    num_base_list = []

    for sample in sample_list:
        num_base_list.append(base_dict[sample])

    num_dog_semibin_pretrain_single['bases'] = num_base_list
    print(num_dog_semibin_pretrain_single)
    sns.lmplot(x = 'bases', y = 'nr_hq', data = num_dog_semibin_pretrain_single)
    plt.title('Dog gut', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('seq_effect_dog.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # soil
    num_ocean_semibin_pretrain_single = get_results_table(dataset='tara', method='SemiBin_pretrain', checkm_only=False)
    print(num_ocean_semibin_pretrain_single['nr_hq'].sum())
    sample_list = num_ocean_semibin_pretrain_single['sample'].values
    ocean_gut_meta = pd.read_csv('updated_results/sequence_depth/ocean.txt', index_col=None, sep=',')
    ocean_gut_meta1 = pd.read_csv('updated_results/sequence_depth/ocean1.txt',index_col=None, sep=',')
    ocean_gut_meta2 = pd.read_csv('updated_results/sequence_depth/ocean2.txt',index_col=None, sep=',')
    ocean_gut_meta3 = pd.read_csv('updated_results/sequence_depth/ocean3.txt',index_col=None, sep=',')
    run_name = ocean_gut_meta['Run'].values
    run_base = ocean_gut_meta['Bases'].values

    run_name1 = ocean_gut_meta1['Run'].values
    run_base1 = ocean_gut_meta1['Bases'].values

    run_name2 = ocean_gut_meta2['Run'].values
    run_base2 = ocean_gut_meta2['Bases'].values

    run_name3 = ocean_gut_meta3['Run'].values
    run_base3 = ocean_gut_meta3['Bases'].values
    run_dict = dict()

    for run, base in zip(run_name, run_base):
        run_dict[run] = base

    for run, base in zip(run_name1, run_base1):
        run_dict[run] = base

    for run, base in zip(run_name2, run_base2):
        run_dict[run] = base

    for run, base in zip(run_name3, run_base3):
        run_dict[run] = base
    num_base_list = []
    ocean_meta = pd.read_excel('updated_results/sequence_depth/tara.xlsx', sheet_name='tableS1', engine='openpyxl', )
    for sample in sample_list:
        if 'lt' in sample:
            sample_ = sample.replace('lt', '<')
        else:
            sample_ = sample
        ocean_meta_ = ocean_meta[ocean_meta['Sample label [TARA_station#_environmental-feature_size-fraction]'] == sample_]
        run_list = ocean_meta_['INSDC run accession number(s)'].values[0]
        run_list = run_list.split('|')
        num_base = 0
        for run in run_list:
            num_base += run_dict[run]
        num_base_list.append(num_base/1e9)

    num_ocean_semibin_pretrain_single['bases'] = num_base_list
    print(num_ocean_semibin_pretrain_single)
    sns.lmplot(x = 'bases', y = 'nr_hq', data = num_ocean_semibin_pretrain_single)
    plt.title('Ocean', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('seq_effect_ocean.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    num_soil_semibin_pretrain_single = get_results_table(dataset='soil', method='SemiBin_pretrain', checkm_only=False)
    print(num_soil_semibin_pretrain_single['nr_hq'].sum())
    sample_list = num_soil_semibin_pretrain_single['sample'].values
    
    # soil
    soil_meta = pd.read_csv('updated_results/sequence_depth/soil.txt', index_col=0, sep='\t')
    sample_name = soil_meta['BioSample'].values
    sample_base = soil_meta['Bases'].values
    base_dict = dict()

    for sample, num_base in zip(sample_name, sample_base):
        if sample not in base_dict:
            base_dict[sample] = num_base/1e9
        else:
            base_dict[sample] += num_base / 1e9

    num_base_list = []

    for sample in sample_list:
        num_base_list.append(base_dict[sample])

    num_soil_semibin_pretrain_single['bases'] = num_base_list
    print(num_soil_semibin_pretrain_single)
    sns.lmplot(x = 'bases', y = 'nr_hq', data = num_soil_semibin_pretrain_single)
    plt.title('Soil', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('seq_effect_soil.pdf', dpi=300, bbox_inches='tight')
    plt.close()


def max_edge_effect():
    def get_hq_num(result):
        binning_result = pd.read_csv(result, index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(
                contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        return len(high_quality)
    human_num_200 = 0
    human_num_500 = 0
    human_num_1000 = 0
    for sample in transfer_single_human:
        human_num_200 += get_hq_num('Results/Real/CheckM_GUNC/human/single_sample/{}/SemiBin/result.csv'.format(sample))
        human_num_500 += get_hq_num('updated_results/effect_max_edge/human/500/{}/result.csv'.format(sample))
        human_num_1000 += get_hq_num('updated_results/effect_max_edge/human/1000/{}/result.csv'.format(sample))
    subset = pd.DataFrame(np.array([[human_num_200,human_num_500,human_num_1000]]),columns = ['200','500','1000'], index=['Human gut'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#1b9e77', '#7570b3','#ec7014'],figsize=(3,4))
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.set_title('Human gut', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_max_edge_human_gut.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    dog_num_200 = 0
    dog_num_500 = 0
    dog_num_1000 = 0
    for sample in transfer_single_dog:
        dog_num_200 += get_hq_num('Results/Real/CheckM_GUNC/dog/single_sample/{}/SemiBin/result.csv'.format(sample))
        dog_num_500 += get_hq_num('updated_results/effect_max_edge/dog/500/{}/result.csv'.format(sample))
        dog_num_1000 += get_hq_num('updated_results/effect_max_edge/dog/1000/{}/result.csv'.format(sample))
    subset = pd.DataFrame(np.array([[dog_num_200,dog_num_500,dog_num_1000]]),columns = ['200','500','1000'], index=['Dog gut'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#1b9e77', '#7570b3','#ec7014'],figsize=(3,4))
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.set_title('Dog gut', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_max_edge_dog_gut.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    ocean_num_200 = 0
    ocean_num_500 = 0
    ocean_num_1000 = 0
    for sample in transfer_single_tara:
        ocean_num_200 += get_hq_num('Results/Real/CheckM_GUNC/tara/single_sample/{}/SemiBin/result.csv'.format(sample))
        ocean_num_500 += get_hq_num('updated_results/effect_max_edge/ocean/500/{}/result.csv'.format(sample))
        ocean_num_1000 += get_hq_num('updated_results/effect_max_edge/ocean/1000/{}/result.csv'.format(sample))
    subset = pd.DataFrame(np.array([[ocean_num_200,ocean_num_500,ocean_num_1000]]),columns = ['200','500','1000'], index=['Ocean'])
    print(subset)
    ax = subset.plot(kind='bar',width = 0.6,color = ['#1b9e77', '#7570b3','#ec7014'],figsize=(3,4))
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    ax.set_title('Ocean', fontsize=20, alpha=1.0,color = 'black')
    plt.savefig('Real_max_edge_ocean.pdf', dpi=300, bbox_inches='tight')
    plt.close()

built_environment_test= ['SAMN03270049', 'SAMN03270027', 'SAMN03270036', 'SAMN03270052', 'SAMN03270038', 'SAMN03270056', 'SAMN03270028', 'SAMN03270031', 'SAMN03270033', 'SAMN03270043']
cat_gut_test = ['SAMEA2150826', 'SAMEA2153824', 'SAMEA2145937', 'SAMEA2144282', 'SAMEA2144745', 'SAMEA2147091', 'SAMEA2142929', 'SAMEA2155252', 'SAMEA2151095', 'SAMEA2143833']
human_oral_test = ['SAMEA2737990', 'SAMEA2737917', 'SAMEA2738035', 'SAMEA2737935', 'SAMEA2738027', 'SAMEA2738037', 'SAMEA2737954', 'SAMEA2738002', 'SAMEA2738014']
mouse_gut_test = ['SAMEA3134376', 'SAMEA3134377', 'SAMEA3134375', 'SAMEA3134368', 'SAMEA3134382', 'SAMEA3134364', 'SAMEA3134369', 'SAMEA3134366', 'SAMEA3134372', 'SAMEA3134385']
pig_gut_test = ['SAMEA3663026', 'SAMEA3663021', 'SAMEA3663031', 'SAMEA3663022', 'SAMEA3663029', 'SAMEA3663013', 'SAMEA3663015', 'SAMEA3663007', 'SAMEA3663010', 'SAMEA3663017']
wastewater_test = ['SAMN04262559', 'SAMN04262501', 'SAMN04262589', 'SAMN04262510', 'SAMN04262568', 'SAMN04262552', 'SAMN04262573', 'SAMN04262500', 'SAMN04262508', 'SAMN04262553']

def extra_env_benchmark():
    def get_hq_num(result):
        binning_result = pd.read_csv(result, index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(
                contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        return len(high_quality)
    Metabat2 = {'cat_gut':0, 'built_environment':0, 'human_oral':0, 'mouse_gut':0, 'pig_gut':0, 'wastewater':0}
    SemiBin = {'cat_gut':0, 'built_environment':0, 'human_oral':0, 'mouse_gut':0, 'pig_gut':0, 'wastewater':0}

    for env in ['cat_gut', 'built_environment', 'human_oral', 'mouse_gut', 'pig_gut', 'wastewater']:
        if env == 'built_environment':
            sample_list = built_environment_test
        if env == 'cat_gut':
            sample_list = cat_gut_test
        if env == 'human_oral':
            sample_list = human_oral_test
        if env == 'mouse_gut':
            sample_list = mouse_gut_test
        if env == 'pig_gut':
            sample_list = pig_gut_test
        if env == 'wastewater':
            sample_list = wastewater_test
        for sample in sample_list:
            Metabat2[env] += get_hq_num('updated_results/extra_env/{0}/Metabat2/{1}/result.csv'.format(env, sample))
            SemiBin[env] += get_hq_num('updated_results/extra_env/{0}/SemiBin/{1}/result.csv'.format(env, sample))
    print(Metabat2)
    print(SemiBin)

    for env in ['cat_gut', 'built_environment', 'human_oral', 'mouse_gut', 'pig_gut', 'wastewater']:
        print(env)
        print((SemiBin[env] - Metabat2[env]) / Metabat2[env])

    subset = pd.DataFrame(np.array([[Metabat2['cat_gut'],SemiBin['cat_gut']],[Metabat2['mouse_gut'],SemiBin['mouse_gut']],[Metabat2['pig_gut'],SemiBin['pig_gut']],[Metabat2['built_environment'],SemiBin['built_environment']],[Metabat2['human_oral'],SemiBin['human_oral']],[Metabat2['wastewater'],SemiBin['wastewater']]]),columns = ['Metabat2','SemiBin'], index=['Cat gut','Mouse gut','Pig gut','Built environment','Human oral','Wastewater'])
    print(subset)

    ax = subset.plot(kind='bar',width = 0.6,color = ['#7570b3', '#1b9e77'])
    ax.set_xticklabels(labels=['Cat gut','Mouse gut','Pig gut','Built environment','Human oral','Wastewater'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    plt.savefig('Extra_env_benchmark.pdf', dpi=300, bbox_inches='tight')
    plt.close()

def plot_cross_validation():
    def get_hq_num(result):
        binning_result = pd.read_csv(result, index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(
                contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        return len(high_quality)

    environment = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater']

    result = []
    Metabat2_result = []
    train_whole_result = []
    SemiBin_pretrain_result = []

    for train_env in environment:
        result_train_env = []
        for test_env in environment:
            if test_env == 'human_gut':
                testing_set = transfer_multi_human
                num_result = 0
                if train_env == 'human_gut':
                    for sample in testing_set:
                        num_result += get_hq_num('Results/Real/CheckM_GUNC/human/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/human/single_sample/{}/Metabat2/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)
                    
                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'dog_gut':
                testing_set = transfer_multi_dog
                num_result = 0
                if train_env == 'dog_gut':
                    for sample in testing_set:
                        num_result += get_hq_num('Results/Real/CheckM_GUNC/dog/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/dog/single_sample/{}/Metabat2/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)
                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'ocean':
                testing_set = transfer_multi_tara
                num_result = 0
                if train_env == 'ocean':
                    for sample in testing_set:
                        num_result += get_hq_num('Results/Real/CheckM_GUNC/tara/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/tara/single_sample/{}/Metabat2/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)
                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'soil':
                testing_set = transfer_multi_soil
                num_result = 0
                if train_env == 'soil':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{}/Metabat2/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'cat_gut':
                testing_set = cat_gut_test
                num_result = 0
                if train_env == 'cat_gut':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/extra_env/cat_gut/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/cat_gut/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'human_oral':
                testing_set = human_oral_test
                num_result = 0
                if train_env == 'human_oral':
                    for sample in testing_set:
                        if sample == 'SAMEA2738039':
                            continue
                        num_result += get_hq_num('updated_results/extra_env/human_oral/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/human_oral/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'mouse_gut':
                testing_set = mouse_gut_test
                num_result = 0
                if train_env == 'mouse_gut':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/extra_env/mouse_gut/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/mouse_gut/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'pig_gut':
                testing_set = pig_gut_test
                num_result = 0
                if train_env == 'pig_gut':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/extra_env/pig_gut/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/pig_gut/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'built_environment':
                testing_set = built_environment_test
                num_result = 0
                if train_env == 'built_environment':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/extra_env/built_environment/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/built_environment/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            if test_env == 'wastewater':
                testing_set = wastewater_test
                num_result = 0
                if train_env == 'wastewater':
                    for sample in testing_set:
                        num_result += get_hq_num('updated_results/extra_env/wastewater/SemiBin/{}/result.csv'.format(sample))
                    result_train_env.append(num_result)
                    SemiBin_pretrain_result.append(num_result)

                    num_result_metabat2 = 0
                    for sample in testing_set:
                        num_result_metabat2 += get_hq_num('updated_results/extra_env/wastewater/Metabat2/{}/result.csv'.format(sample))
                    Metabat2_result.append(num_result_metabat2)

                    num_result_whole = 0
                    for sample in testing_set:
                        num_result_whole += get_hq_num('updated_results/training_whole/{0}/{1}/result.csv'.format(train_env, sample))
                    train_whole_result.append(num_result_whole)
                    continue

            num_result = 0
            for sample in testing_set:
                if sample == 'SAMEA2738039':
                    continue
                num_result += get_hq_num('updated_results/cross_validation/{0}/{1}/{2}/result.csv'.format(train_env, test_env, sample))
            result_train_env.append(num_result)
        result.append(result_train_env)

    result.insert(0, train_whole_result)
    result.insert(0, Metabat2_result)
    result = np.array(result)

    human_names = {
            'human_gut': 'Human gut',
            'dog_gut': 'Dog gut',
            'ocean': 'Ocean',
            'soil': 'Soil',
            'cat_gut': 'Cat gut',
            'human_oral': 'Human oral',
            'mouse_gut': 'Mouse gut',
            'pig_gut': 'Pig gut',
            'built_environment': 'Built environment',
            'wastewater': 'Wastewater',
            }

    df = pd.DataFrame(result, index=['Metabat2', 'SemiBin(whole)'] + environment, columns=environment)
    df.rename(index=human_names, columns=human_names, inplace=True)
    order1 = [
           'Human gut',
           'Dog gut',
           'Cat gut',
           'Pig gut',
           'Human oral',
           'Mouse gut',
           'Ocean',
           'Built environment',
           'Soil',
           'Wastewater',
           ]
    order0 = ['Metabat2', 'SemiBin(whole)']+order1
    same_env = pd.Series({e:df.loc[e,e] for e in order1})
    df = df.loc[order0, order1]

    df_norm = (df+10).divide(same_env+10)

    fig,(ax_top,ax_main) = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios':[2,10]})
    im0 = ax_top.imshow(
            df_norm.iloc[:2]*100.,
            vmin=100*df_norm.min().min(),
            vmax=100*df_norm.max().max(),
            aspect='auto', cmap=cm.YlOrBr)
    im1 = ax_main.imshow(
            df_norm.iloc[2:]*100.,
            vmin=100*df_norm.min().min(),
            vmax=100*df_norm.max().max(),
            aspect='auto', cmap=cm.YlOrBr)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            color = ('#333333' if df_norm.iloc[i,j] < .90 else '#eeeeee')
            (ax_main if i > 1 else ax_top).text(s='{}'.format(df.iloc[i,j]),
                    y=(i if i < 2 else i - 2),
                    x=j,
                    color=color,
                    horizontalalignment='center',
                    verticalalignment='center')
    ax_main.set_xticks(np.arange(len(order1)))
    ax_main.set_xticklabels(order1, rotation=90)
    ax_main.set_yticks(np.arange(len(order0[2:])))
    ax_main.set_yticklabels(order0[2:])
    ax_top.set_yticks(np.arange(len(order0[:2])))
    ax_top.set_yticklabels(order0[:2])
    fig.colorbar(im0, orientation='horizontal', )
    fig.tight_layout()
    fig.savefig('model-transfer.svg')


def plot_training_whole():
    def get_hq_num(result):
        binning_result = pd.read_csv(result, index_col=0)
        high_quality = binning_result[(binning_result['Completeness'].astype(float) > float(90)) & (
                    binning_result['Contamination'].astype(float) < float(
                contamination * 100)) & (binning_result['pass.GUNC'] == True)]
        return len(high_quality)

    environment = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater']
    Metabat2_result = []
    SemiBin_pretrain_result = []
    SemiBin_whole_rsult = []

    for test_env in environment:
        if test_env == 'human_gut':
            testing_set = transfer_multi_human
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/human/single_sample/{}/Metabat2/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('Results/Real/CheckM_GUNC/human/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/human_gut/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'dog_gut':
            testing_set = transfer_multi_dog
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/dog/single_sample/{}/Metabat2/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('Results/Real/CheckM_GUNC/dog/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/dog_gut/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'ocean':
            testing_set = transfer_multi_tara
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('Results/Real/CheckM_GUNC/tara/single_sample/{}/Metabat2/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('Results/Real/CheckM_GUNC/tara/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/ocean/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'soil':
            testing_set = transfer_multi_soil
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{}/Metabat2/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/Soil_benchmark/CheckM_GUNC/single_sample/{}/SemiBin_pretrain/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/soil/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'cat_gut':
            testing_set = cat_gut_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/cat_gut/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/cat_gut/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/cat_gut/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'human_oral':
            testing_set = human_oral_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/human_oral/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/human_oral/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/human_oral/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'mouse_gut':
            testing_set = mouse_gut_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/mouse_gut/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/mouse_gut/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/mouse_gut/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'pig_gut':
            testing_set = pig_gut_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/pig_gut/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/pig_gut/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/pig_gut/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'built_environment':
            testing_set = built_environment_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/built_environment/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/built_environment/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/built_environment/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

        if test_env == 'wastewater':
            testing_set = wastewater_test
            num_result_metabat2 = 0
            num_result_pertrain = 0
            num_result_whole = 0
            for sample in testing_set:
                num_result_metabat2 += get_hq_num('updated_results/extra_env/wastewater/Metabat2/{}/result.csv'.format(sample))
                num_result_pertrain += get_hq_num('updated_results/extra_env/wastewater/SemiBin/{}/result.csv'.format(sample))
                num_result_whole += get_hq_num('updated_results/training_whole/wastewater/{}/result.csv'.format(sample))
            Metabat2_result.append(num_result_metabat2)
            SemiBin_pretrain_result.append(num_result_pertrain)
            SemiBin_whole_rsult.append(num_result_whole)

    result = []
    result.append(Metabat2_result)
    result.append(SemiBin_pretrain_result)
    result.append(SemiBin_whole_rsult)
    result = np.array(result).T
    subset = pd.DataFrame(result,columns = ['Metabat2','SemiBin(pertrain)','SemiBin(whole)'], index=['Human gut', 'Dog gut', 'ocean', 'Soil', 'Cat gut', 'Human oral', 'Mouse gut', 'Pig gut', 'Built environment', 'Wastewater'])
    print(subset)

    ax = subset.plot(kind='bar',width = 0.6,color = ['#7570b3', '#41AB5D', '#005A32'])
    # ax.set_yticks(ticks=[0,800,1600,2400,3200])
    # ax.set_yticklabels(labels=[0,800,1600,2400,3200],fontsize=12,color = 'black')
    ax.set_xticklabels(labels=['Human gut', 'Dog gut', 'ocean', 'Soil', 'Cat gut', 'Human oral', 'Mouse gut', 'Pig gut', 'Built environment', 'Wastewater'], fontsize=15,color = 'black',rotation = 50)
    ax.set_ylabel('High-quality bins', fontsize=15,color = 'black')
    plt.savefig('Whole_training.pdf', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # print('human')
    # plot_bar_per_sample_com('human',[-10,0,10,20,30],[0,500,1000,1500])
    # print('dog')
    # plot_bar_per_sample_com('dog',[-10,0,10,20,30],[0,1000,2000,3000])
    # print('tara')
    # plot_bar_per_sample_com('tara',[-10,0,10,20,30], [0,200,400,600])
    # print('soil')
    # plot_bar_per_sample_com('soil', [-10, 0, 10, 20, 30], [0, 60, 120, 180])

    # tranfer_multi()

    ### bar plot high quality genomes comparison
    # plot_checkm_high_quality_comparison()

    ### venn plot multi annotation comparison
    # plot_multi_venn_comparison()

    #alternative_real_compare_test()
    # plot_sankey_overlap(dataset='human',output='human_sankey.pdf')
    # plot_sankey_overlap(output='dog_sankey.pdf')
    # plot_sankey_overlap(dataset='tara',output='tara_sankey.pdf')
    # plot_sankey_overlap(dataset='soil', output='soil_sankey.pdf')

    ### recall, precision, F1-score box plot
    #
    # plot_overlap_F1('human')
    # plot_overlap_F1()
    # plot_overlap_F1('tara')
    # plot_overlap_F1('soil')

    # ### bar plot the overlap of annotation in all taxa
    # plot_taxa_overlap()

    # ## comparison of known and unknown species
    # plot_comparison_known_unknown(y_label=[0,500,1000,1500,2000,2500], output='dog_taxi_known_unknown.pdf')
    # plot_comparison_known_unknown(dataset='human', y_label=[0,300,600,900,1200,1500], output='human_taxi_known_unknown.pdf')
    # plot_comparison_known_unknown(dataset='tara', y_label=[0,100,200,300,400], output='tara_taxi_known_unknown.pdf')
    # plot_comparison_known_unknown(dataset='soil',y_label=[0, 20, 40, 60, 80, 100],output='soil_taxi_known_unknown.pdf')

    # plot_transfer()

    # plot_extra_bar()

    # plot_extra_per_sample()

    # CAT_mmseqs()

    # seq_depth_effect()

    # max_edge_effect()

    # extra_env_benchmark()

    plot_cross_validation()

    # plot_training_whole()

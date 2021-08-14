import os

human_list = ['CCMD41521570ST', 'CCMD75147712ST', 'CCMD18579000ST', 'CCMD53508245ST', 'CCMD19168690ST', 'CCMD52117727ST', 'CCMD42956136ST', 'CCMD79349503ST', 'CCMD89306485ST', 'CCMD76409700ST', 'CCMD31134579ST', 'CCMD71242853ST', 'CCMD89107682ST', 'CCMD76222476ST', 'CCMD10032470ST', 'CCMD17410933ST', 'CCMD38158721ST', 'CCMD35081859ST', 'CCMD54057834ST', 'CCMD28738636ST', 'CCMD98702133ST', 'CCMD30626189ST', 'CCMD32965613ST', 'CCMD53522274ST', 'CCMD37575804ST', 'CCMD68973846ST', 'CCMD25475945ST', 'CCMD65406197ST', 'CCMD21703880ST', 'CCMD50300306ST', 'CCMD51228890ST', 'CCMD59540613ST', 'CCMD49942357ST', 'CCMD95431029ST', 'CCMD41202658ST', 'CCMD15562448ST', 'CCMD21593359ST', 'CCMD92404903ST', 'CCMD50538120ST', 'CCMD49461418ST', 'CCMD72690923ST', 'CCMD85481373ST', 'CCMD39882286ST', 'CCMD18829815ST', 'CCMD51154251ST', 'CCMD85661207ST', 'CCMD71915439ST', 'CCMD39157124ST', 'CCMD22852639ST', 'CCMD35801800ST', 'CCMD27463710ST', 'CCMD59583015ST', 'CCMD89967135ST', 'CCMD52145360ST', 'CCMD95676152ST', 'CCMD45004878ST', 'CCMD67373733ST', 'CCMD99929634ST', 'CCMD89643949ST', 'CCMD26625622ST', 'CCMD23541216ST', 'CCMD31009081ST', 'CCMD99440714ST', 'CCMD66848156ST', 'CCMD65222621ST', 'CCMD98531134ST', 'CCMD45812507ST', 'CCMD46727384ST', 'CCMD73128545ST', 'CCMD30627121ST', 'CCMD50529145ST', 'CCMD98198513ST', 'CCMD93755960ST', 'CCMD35633353ST', 'CCMD56948710ST', 'CCMD27867141ST', 'CCMD32288175ST', 'CCMD29706695ST', 'CCMD72666896ST', 'CCMD10191450ST', 'CCMD49025643ST', 'CCMD74592084ST']

dog_list = ['SAMN06172456', 'SAMN06172425', 'SAMN06172487', 'SAMN06172450', 'SAMN06172459', 'SAMN06172479', 'SAMN06172435', 'SAMN06172414', 'SAMN06172409', 'SAMEA103957796', 'SAMN06172442', 'SAMN06172500', 'SAMN06172437', 'SAMN06172413', 'SAMN06172514', 'SAMN06172403', 'SAMN06172471', 'SAMN06172490', 'SAMN06172448', 'SAMN06172504', 'SAMN06172457', 'SAMN06172441', 'SAMN06172422', 'SAMN06172408', 'SAMN06172429', 'SAMN06172420', 'SAMN06172503', 'SAMN06172410', 'SAMN06172458', 'SAMN06172493', 'SAMEA103957794', 'SAMN06172402', 'SAMN06172515', 'SAMN06172462', 'SAMN06172421', 'SAMN06172411', 'SAMN06172511', 'SAMN06172516', 'SAMN06172465', 'SAMN06172419', 'SAMN06172517', 'SAMN06172510', 'SAMN06172418', 'SAMN06172424', 'SAMN06172427', 'SAMN06172453', 'SAMN06172491', 'SAMN06172496', 'SAMN06172513', 'SAMN06172461', 'SAMN06172449', 'SAMN06172426', 'SAMN06172452', 'SAMN06172522', 'SAMN06172400', 'SAMN06172405', 'SAMN06172521', 'SAMN06172407', 'SAMN06172455', 'SAMN06172446', 'SAMN06172467', 'SAMN06172499', 'SAMN06172474', 'SAMN06172412', 'SAMN06172468', 'SAMN06172478', 'SAMN06172423', 'SAMN06172447', 'SAMN06172415', 'SAMN06172523', 'SAMN06172417', 'SAMN06172497', 'SAMN06172498', 'SAMN06172489', 'SAMN06172436', 'SAMN06172432', 'SAMN06172406', 'SAMN06172488', 'SAMN06172502', 'SAMN06172401', 'SAMN06172434', 'SAMN06172416', 'SAMN06172445', 'SAMN06172431', 'SAMN06172438', 'SAMN06172473', 'SAMN06172486', 'SAMN06172472', 'SAMN06172428', 'SAMEA103957793', 'SAMEA103957795', 'SAMN06172443', 'SAMN06172475', 'SAMN06172520', 'SAMN06172495', 'SAMN06172440', 'SAMN06172430', 'SAMN06172481', 'SAMN06172524', 'SAMN06172519', 'SAMN06172454', 'SAMN06172404', 'SAMN06172460', 'SAMN06172433', 'SAMN06172469', 'SAMN06172451', 'SAMN06172476', 'SAMN06172492', 'SAMN06172484', 'SAMN06172509', 'SAMN06172506', 'SAMN06172518', 'SAMN06172477', 'SAMN06172470', 'SAMN06172482', 'SAMN06172512', 'SAMN06172494', 'SAMN06172485', 'SAMN06172508', 'SAMN06172466', 'SAMN06172507', 'SAMN06172444', 'SAMN06172505', 'SAMN06172464', 'SAMN06172439', 'SAMN06172501', 'SAMN06172483', 'SAMN06172463', 'SAMN06172480']

ocean_list = ['TARA_041_SRF_0.1-0.22', 'TARA_038_SRF_0.22-1.6', 'TARA_076_SRF_0.22-3', 'TARA_023_SRF_0.22-1.6', 'TARA_042_SRF_0.22-1.6', 'TARA_124_SRF_0.22-3', 'TARA_124_SRF_0.22-0.45', 'TARA_066_SRF_lt-0.22', 'TARA_057_SRF_0.22-3', 'TARA_124_SRF_0.45-0.8', 'TARA_004_SRF_0.22-1.6', 'TARA_018_SRF_0.22-1.6', 'TARA_070_SRF_0.22-0.45', 'TARA_034_SRF_lt-0.22', 'TARA_064_SRF_0.22-3', 'TARA_125_SRF_0.22-0.45', 'TARA_111_SRF_0.22-3', 'TARA_122_SRF_0.22-0.45', 'TARA_145_SRF_0.22-3', 'TARA_099_SRF_0.22-3', 'TARA_038_SRF_lt-0.22', 'TARA_082_SRF_0.22-3', 'TARA_041_SRF_lt-0.22', 'TARA_146_SRF_0.22-3', 'TARA_151_SRF_0.22-3', 'TARA_123_SRF_0.22-3', 'TARA_110_SRF_0.22-3', 'TARA_150_SRF_0.22-3', 'TARA_072_SRF_lt-0.22', 'TARA_085_SRF_0.22-3', 'TARA_098_SRF_0.22-3', 'TARA_078_SRF_0.22-3', 'TARA_149_SRF_0.22-3', 'TARA_094_SRF_0.22-3', 'TARA_068_SRF_0.22-3', 'TARA_148_SRF_0.22-3', 'TARA_067_SRF_0.45-0.8', 'TARA_018_SRF_lt-0.22', 'TARA_138_SRF_0.22-3', 'TARA_093_SRF_0.22-3', 'TARA_041_SRF_0.22-1.6', 'TARA_122_SRF_0.22-3', 'TARA_078_SRF_0.45-0.8', 'TARA_070_SRF_lt-0.22', 'TARA_065_SRF_lt-0.22', 'TARA_122_SRF_0.1-0.22', 'TARA_036_SRF_0.22-1.6', 'TARA_031_SRF_0.22-1.6', 'TARA_142_SRF_0.22-3', 'TARA_124_SRF_0.1-0.22', 'TARA_036_SRF_lt-0.22', 'TARA_065_SRF_0.22-3', 'TARA_067_SRF_lt-0.22', 'TARA_112_SRF_0.22-3', 'TARA_109_SRF_0.22-3', 'TARA_068_SRF_lt-0.22', 'TARA_109_SRF_lt-0.22', 'TARA_064_SRF_lt-0.22', 'TARA_048_SRF_0.22-1.6', 'TARA_034_SRF_0.22-1.6', 'TARA_070_SRF_0.45-0.8', 'TARA_025_SRF_lt-0.22', 'TARA_133_SRF_0.22-3', 'TARA_096_SRF_0.22-3', 'TARA_038_SRF_0.1-0.22', 'TARA_007_SRF_0.22-1.6', 'TARA_048_SRF_0.1-0.22', 'TARA_140_SRF_0.22-3', 'TARA_034_SRF_0.1-0.22', 'TARA_067_SRF_0.22-3', 'TARA_125_SRF_0.45-0.8', 'TARA_030_SRF_0.22-1.6', 'TARA_031_SRF_lt-0.22', 'TARA_032_SRF_0.22-1.6', 'TARA_070_SRF_0.22-3', 'TARA_132_SRF_0.22-3', 'TARA_076_SRF_lt-0.22', 'TARA_125_SRF_0.1-0.22', 'TARA_123_SRF_0.45-0.8', 'TARA_078_SRF_lt-0.22', 'TARA_068_SRF_0.45-0.8', 'TARA_068_SRF_0.22-0.45', 'TARA_067_SRF_0.22-0.45', 'TARA_100_SRF_0.22-3', 'TARA_122_SRF_0.45-0.8', 'TARA_137_SRF_0.22-3', 'TARA_076_SRF_0.22-0.45', 'TARA_125_SRF_0.22-3', 'TARA_078_SRF_0.22-0.45', 'TARA_076_SRF_0.45-0.8', 'TARA_084_SRF_0.22-3', 'TARA_032_SRF_lt-0.22', 'TARA_025_SRF_0.22-1.6', 'TARA_062_SRF_0.22-3', 'TARA_066_SRF_0.22-3', 'TARA_036_SRF_0.1-0.22', 'TARA_056_SRF_0.22-3', 'TARA_072_SRF_0.22-3', 'TARA_128_SRF_0.22-3', 'TARA_052_SRF_0.22-1.6', 'TARA_033_SRF_0.22-1.6', 'TARA_123_SRF_0.22-0.45', 'TARA_102_SRF_0.22-3', 'TARA_065_SRF_0.1-0.22', 'TARA_009_SRF_0.22-1.6', 'TARA_141_SRF_0.22-3', 'TARA_045_SRF_0.22-1.6', 'TARA_042_SRF_lt-0.22', 'TARA_152_SRF_0.22-3']

def run_single(dataset = None):
    if dataset == 'human':
        sample_list = human_list
        env = 'human_gut'
    elif dataset == 'dog':
        sample_list = dog_list
        env = 'dog_gut'
    elif dataset == 'ocean':
        sample_list = ocean_list
        env = 'ocean'

    for sample in sample_list:
        bam_file = '{}.mapped.sorted.bam'.format(sample)
        if dataset != 'ocean':
            contig_file = '{}-assembled.fa'.format(sample)
        else:
            contig_file = '{}_G-assembled.fa'.format(sample)
        # SemiBin
        os.system('SemiBin single_easy_bin -i {0} -o {1}/{2}/SemiBin -b {3} --recluster'.format(contig_file, dataset, sample, bam_file))

        # SemiBin(pretrain)
        os.system('SemiBin single_easy_bin -i {0} -o {1}/{2}/SemiBin_pretrain -b {3} --recluster --environment {4}'.format(contig_file, dataset, sample, bam_file, env))

        # Metabat2
        os.system('runMetaBat.sh {0} {1}'.format(contig_file, bam_file))

        # VAMB
        os.system('vamb --outdir {1}/{2}/VAMB --fasta {0} --bamfiles {3} --minfasta 200000 -m 2000 --cuda'.format(contig_file, dataset, sample, bam_file))

        # Maxbin2
        # For running Maxbin2, first we generated the abundance file from the bam file and then we running Maxbin2
        os.makedirs('{0}/{1}/Maxbin2', exist_ok=True)
        os.system('perl run_MaxBin.pl -contig {0} -out {1}/{2}/Maxbin2 -abund {1}/{2}/coverage.tsv'.format(contig_file, dataset, sample))

def run_multi(dataset = None):
    if dataset == 'human':
        sample_list = human_list
    elif dataset == 'dog':
        sample_list = dog_list
    elif dataset == 'ocean':
        sample_list = ocean_list

    concatenated_fasta = 'concatenated,fasta'
    bam_files = '*.mapped.sorted,bam'

    # VAMB
    if dataset == 'human':
        os.system('vamb --outdir {1}/VAMB_multi --fasta {0} --bamfiles {2} --minfasta 200000 -m 2000 --cuda -o TC'.format(concatenated_fasta, dataset, bam_files))
    else:
        os.system('vamb --outdir {1}/VAMB_multi --fasta {0} --bamfiles {2} --minfasta 200000 -m 2000 --cuda -o C'.format(concatenated_fasta, dataset, bam_files))

    # SemiBin(We can use the cannot_link generated in single-sample binning)

    if dataset == 'human':
        os.system('SemiBin generate_data_multi -i {0} -b {1} -s TC -o {2}/SemiBin_multi'.format(concatenated_fasta, bam_files, dataset))
    else:
        os.system('SemiBin generate_data_multi -i {0} -b {1} -s C -o {2}/SemiBin_multi'.format(concatenated_fasta, bam_files, dataset))

    for sample in sample_list:
        if dataset != 'ocean':
            contig_file = '{}-assembled.fa'.format(sample)
        else:
            contig_file = '{}_G-assembled.fa'.format(sample)

        # train model
        os.system('SemiBin train --data {0}/SemiBin_multi/samples/{1}/data.csv --data-split {0}/SemiBin_multi/samples/{1}/data_split.csv -c {0}/{1}/SemiBin/cannot/cannot.txt --mode single -i {2} -o {0}/SemiBin_multi/samples/{1} -b {3}'.format(dataset, sample, contig_file, bam_files))
        # bin
        os.system('SemiBin bin --data {0}/SemiBin_multi/samples/{1}/data.csv --recluster --model {0}/SemiBin_multi/samples/{1}/model.h5 -i {2} -o {0}/SemiBin_multi/samples/{1} -b {3}'.format(dataset, sample, contig_file, bam_files))

if __name__ == '__main__':
    run_single('human')
    run_single('dog')
    run_single('ocean')

    run_multi('human')
    run_multi('dog')
    run_multi('ocean')


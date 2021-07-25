import os


def run_CAMI_I(dataset = None):
    if dataset == 'low':
        fasta_file = 'CAMI_low_RL_S001__insert_270_GoldStandardAssembly.fasta'
        bam_files = 'RL_S001__insert_270.mapped.sorted.bam'
        output = 'Low'
        n_sample = 1
    elif dataset == 'medium':
        fasta_file = 'CAMI_medium_GoldStandardAssembly.fasta'
        bam_files = 'RM2_S00*__insert_270.mapped.sorted.bam'
        output = 'Medium'
        n_sample = 2
    elif dataset == 'high':
        fasta_file = 'CAMI_high_GoldStandardAssembly.fasta'
        bam_files = 'RH_S00*__insert_270.mapped.sorted.bam'
        output = 'High'
        n_sample = 5

    # SemiBin
    os.system('SemiBin single_easy_bin -i {0} -o {1}/SemiBin_output -b {2} --recluster'.format(fasta_file, output, bam_files))

    # modified version of SemiBin

    # Generate data.csv for these versions
    os.system('SemiBin generate_data_single -i {0} -o {1}/generalization -b {2}'.format(fasta_file, output, bam_files))

    # Generate must_link and cannot_link constraints
    os.system('mmseqs createdb {0} {1}/mmseqs_annotation/contig_DB'.format(fasta_file, output))
    os.system('mmseqs taxonomy {0}/mmseqs_annotation/contig_DB ~/.cache/SemiBin/mmseqs2-GTDB/GTDB {0}/mmseqs_annotation/mmseqs_contig_annotation {0}/mmseqs_annotation_tmp --tax-lineage 1'.format(output))
    os.system('mmseqs createtsv {0}/mmseqs_annotation/contig_DB {0}/mmseqs_annotation/mmseqs_contig_annotation {0}/mmseqs_annotation/taxonomyResult.tsv'.format(output))
    os.system('python generate_constraints.py -i {0}/mmseqs_annotation/taxonomyResult.tsv -c {1} -o {0}/mmseqs_annotation/ --mmseqs'.format(output, fasta_file))

    # NoSemi
    os.system('python SemiBin_generalization.py -i {0} -o {1}/generalization/NoSemi --data {1}/generalization/data.csv -n {2}'.format(fasta_file, output, n_sample))
    # SemiBin_m
    os.system('python SemiBin_generalization.py -i {0} -o {1}/generalization/SemiBin_m --data {1}/generalization/data.csv -m {1}/mmseqs_annotation/must.txt -n {2}'.format(fasta_file, output, n_sample))
    # SemiBin_c
    os.system('python SemiBin_generalization.py -i {0} -o {1}/generalization/SemiBin_c --data {1}/generalization/data.csv -c {1}/mmseqs_annotation/cannot.txt -n {2}'.format(fasta_file, output, n_sample))
    # SemiBin_mc
    os.system('python SemiBin_generalization.py -i {0} -o {1}/generalization/SemiBin_mc --data {1}/generalization/data.csv -m {1}/mmseqs_annotation/must.txt -c {1}/mmseqs_annotation/cannot.txt -n {2}'.format(fasta_file, output, n_sample))

    # Metabat2
    os.system('runMetaBat.sh {0} {1}'.format(fasta_file, bam_files))

    # Maxbin2
    if dataset == 'low':
        os.system('perl run_MaxBin.pl -contig {0} -out {1}/Maxbin2_output -reads RL_S001__insert_270.fq.gz'.format(fasta_file, output))
    elif dataset == 'medium':
        os.system('perl run_MaxBin.pl -contig {0} -out {1}/Maxbin2_output -reads RM2_S001__insert_270.fq.gz -reads1 RM2_S002__insert_270.fq.gz'.format(fasta_file, output))
    elif dataset == 'high':
        os.system('perl run_MaxBin.pl -contig {0} -out {1}/Maxbin2_output -reads RH_S001__insert_270.fq.gz -reads1 RH_S002__insert_270.fq.gz -reads2 RH_S003__insert_270.fq.gz -reads3 RH_S004__insert_270.fq.gz -reads4 RH_S005__insert_270.fq.gz'.format(fasta_file, output))

    # We need to generate composition_profiles and coverage_profiles for SolidBin and COCACOLA first
    # We used the script from SolidBin
    # composition_profiles: https://github.com/sufforest/SolidBin/blob/master/scripts/gen_kmer.py
    # coverage_profiles:
    # https://github.com/sufforest/SolidBin/blob/master/scripts/gen_cov.sh

    os.system('python gen_kmer.py {0} 1000 4 {1}/kmer'.format(fasta_file, output))
    os.system('bedtools genomecov -ibam {0} > {1}/coverage/RL_S001__insert_270_cov.txt'.format(bam_files, output))
    os.system('bash gen_cov.sh {0}/coverage'.format(output))
    os.system('perl Collate.pl {0}/coverage > {0}/coverage/coverage.csv'.format(output))
    os.system(
        'perl -pe "s/,/\t/g;" {0}/coverage/coverage.csv > {0}/coverage/coverage.tsv'.format(output))

    # COCACOLA
    os.system('python cocacola.py --contig_file {0} --abundance_profiles {1}/coverage/coverage.tsv --composition_profiles {1}/kmer/kmer.csv --output {1}/COCACOLA_output/result.csv'.format(fasta_file, output))

    # SolidBin-naive
    os.system('python SolidBin.py --contig_file {0} --coverage_profiles {1}/coverage/coverage.tsv --composition_profiles {1}/kmer/kmer.csv --output {1}/SolidBin_naive_output/result.tsv'.format(fasta_file, output))

    # Generate must_link and cannot_link for SolidBin
    os.system('python generate_constraints.py -i {0}/mmseqs_annotation/taxonomyResult.tsv -c {1} -o {0}/mmseqs_constraints_solidbin/ --mmseqs --solidbin'.format(output, fasta_file))

    # SolidBin_coalign
    os.system('python SolidBin.py --contig_file {0} --coverage_profiles {1}/coverage/coverage.tsv --composition_profiles {1}/kmer/kmer.csv --output {1}/SolidBin_coalign_output/result.tsv  --priori_ml_list {1}/mmseqs_constraints_solidbin/must.txt'.format(fasta_file, output))

    # SolidBin_CL
    os.system('python SolidBin.py --contig_file {0} --coverage_profiles {1}/coverage/coverage.tsv --composition_profiles {1}/kmer/kmer.csv --output {1}/SolidBin_CL_output/result.tsv --priori_cl_list {1}/mmseqs_constraints_solidbin/cannot.txt'.format(fasta_file, output))

    # SolidBin_SFS_CL
    os.system('python SolidBin.py --contig_file {0} --coverage_profiles {1}/coverage/coverage.tsv --composition_profiles {1}/kmer/kmer.csv --output {1}/SolidBin_SFS_CL_output/result.tsv --priori_cl_list {1}/mmseqs_constraints_solidbin/cannot.txt --use_sfs'.format(fasta_file, output))


if __name__ == '__main__':
    run_CAMI_I('low')
    run_CAMI_I('medium')
    run_CAMI_I('high')





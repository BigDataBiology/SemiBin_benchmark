import os

skin_fasta = 'skin_contig_whole.fna'  # contig from 10 samples
skin_bamfiles = 'skin_*.mapped.sorted.bam'  # 10 bam files

oral_fasta = 'oral_contig_whole.fna'  # contig from 10 samples
oral_bamfiles = 'oral_*.mapped.sorted.bam'  # 10 bam files

n_sample = 10

def run_multi_sample():
    # VAMB multi skin
    os.system('vamb --outdir VAMB_skin --fasta {0} --bamfiles {1} -o C --minfasta 200000 --cuda -m 2000'.format(skin_fasta, skin_bamfiles))

    # VAMB multi oral
    os.system('vamb --outdir VAMB_oral --fasta {0} --bamfiles {1} -o C --minfasta 200000 --cuda -m 2000'.format(oral_fasta, oral_bamfiles))

    # SemiBin multi skin
    os.system('SemiBin multi_easy_bin -i {0} -o SemiBin_skin -b {1} --recluster -s C'.format(skin_fasta, skin_bamfiles))

    # SemiBin multi oral
    os.system('SemiBin multi_easy_bin -i {0} -o SemiBin_oral -b {1} --recluster -s C'.format(oral_fasta, oral_bamfiles))

def run_modified_SemiBin():
    # Skin
    os.system('SemiBin generate_data_multi -i {0} -b {1} -s C -o SemiBin_skin_modified'.format(skin_fasta, skin_bamfiles))

    for index in [1,13,14,15,16,17,18,19,20,28]:
        # generate must_link and cannot_link constraints
        # here we can use the taxonomyResult.tsv from the multi_sample binning that running before
        os.system('python generate_constraints.py -i SemiBin_skin/samples/S{0}/mmseqs_annotation/taxonomyResult.tsv -c contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation --mmseqs'.format(index))

        # NoSemi
        os.system('python SemiBin_generalization.py -i skin_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/NoSemi --data SemiBin_skin_modified/samples/S{0}/data.csv -n {1}'.format(index, n_sample))
        # SemiBin_m
        os.system('python SemiBin_generalization.py -i skin_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/SemiBin_m --data SemiBin_skin_modified/samples/S{0}/data.csv -m SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation/must.txt -n {1}'.format(index, n_sample))
        # SemiBin_c
        os.system('python SemiBin_generalization.py -i skin_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/SemiBin_c --data SemiBin_skin_modified/samples/S{0}/data.csv -c SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation/cannot.txt -n {1}'.format(index, n_sample))
        # SemiBin_mc
        os.system('python SemiBin_generalization.py -i skin_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/SemiBin_mc --data SemiBin_skin_modified/samples/S{0}/data.csv -m SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation/must.txt -c SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation/cannot.txt -n {1}'.format(index, n_sample))


    # Oral
    os.system('SemiBin generate_data_multi -i {0} -b {1} -s C -o SemiBin_oral_modified'.format(oral_fasta, oral_bamfiles))

    for index in [6,7,8,13,14,15,16,17,18,19]:
        os.system('python generate_constraints.py -i SemiBin_oral/samples/S{0}/mmseqs_annotation/taxonomyResult.tsv -c contigs/S{0}/anonymous_gsa.fasta -o SemiBin_skin_modified/samples/S{0}/mmseqs_generalziation --mmseqs'.format(index))

        # NoSemi
        os.system('python SemiBin_generalization.py -i oral_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_oral_modified/samples/S{0}/NoSemi --data SemiBin_oral_modified/samples/S{0}/data.csv -n {1}'.format(index, n_sample))
        # SemiBin_m
        os.system('python SemiBin_generalization.py -i oral_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_oral_modified/samples/S{0}/SemiBin_m --data SemiBin_oral_modified/samples/S{0}/data.csv -m SemiBin_oral_modified/samples/S{0}/mmseqs_generalziation/must.txt -n {1}'.format(index, n_sample))
        # SemiBin_c
        os.system('python SemiBin_generalization.py -i oral_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_oral_modified/samples/S{0}/SemiBin_c --data SemiBin_oral_modified/samples/S{0}/data.csv -c SemiBin_oral_modified/samples/S{0}/mmseqs_generalziation/cannot.txt -n {1}'.format(index, n_sample))
        # SemiBin_mc
        os.system('python SemiBin_generalization.py -i oral_contigs/S{0}/anonymous_gsa.fasta -o SemiBin_oral_modified/samples/S{0}/SemiBin_mc --data SemiBin_oral_modified/samples/S{0}/data.csv -m SemiBin_oral_modified/samples/S{0}/mmseqs_generalziation/must.txt -c SemiBin_oral_modified/samples/S{0}/mmseqs_generalziation/cannot.txt -n {1}'.format(index, n_sample))

if __name__ == '__main__':
    run_multi_sample()
    run_modified_SemiBin()
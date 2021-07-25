import argparse
import os
import sys
import numpy as np
import logging
from  Bio import SeqIO
import pandas as pd
import subprocess
from sklearn.cluster import KMeans
from atomicwrites import atomic_write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import shutil
from sklearn.neighbors import kneighbors_graph
from igraph import Graph
import warnings



def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Semi-supervised siamese neural network for metagenomic binning')
    parser.add_argument('-n',required=True, dest='n_sample',type=int)
    parser.add_argument('--data', required=True, dest='data', type=str)
    parser.add_argument('-i','--input-fasta',
                        required=True,
                        help='Path to the input contig fasta file.',
                        dest='contig_fasta',
                        default=None)
    parser.add_argument('-c','--cannot-link',
                        required=False,
                        help='Path to the input can not link file generated from other additional biological information,one row for one can not link                               constraint.The file format:contig_1\tcontig_2.',
                        dest='cannot_link',
                        default=None)
    parser.add_argument('-m','--must-link',
                        required=False,
                        help='Path to the input can not link file generated from other additional biological information,one row for one can not link                               constraint.The file format:contig_1\tcontig_2.',
                        dest='must_link',
                        default=None)
    parser.add_argument('-o','--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default=None)
    parser.add_argument('--epoches',
                        required=False,
                        type=int,
                        help='Epoches used in the training process.',
                        dest='epoches',
                        default=20
    )
    parser.add_argument('--batch-size',
                        required=False,
                        type=int,
                        help='Batch size used in the training process.',
                        dest='batchsize',
                        default=2048,
                        )
    parser.add_argument('--max-edges',
                        required=False,
                        type=int,
                        help='The maximun number of edges that can be connected to one contig.',
                        dest='max_edges',
                        default=200)
    parser.add_argument('--max-node',
                        required=False,
                        type=float,
                        dest='max_node',
                        default=1,
                        help='Percentage of contigs that considered to be binned.')

    return parser.parse_args()

def validate_args(args):

    def except_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)
    except_file(args.contig_fasta)


def get_threshold(contig_len):
    """
    calculate the threshold length for must link breaking up
    """
    basepair_sum = 0
    threshold = 0
    whole_len = np.sum(contig_len)
    contig_len.sort(reverse = True)
    index = 0
    while(basepair_sum / whole_len < 0.98):
        basepair_sum += contig_len[index]
        threshold = contig_len[index]
        index += 1
    threshold = max(threshold, 4000)
    return threshold

def write_bins(namelist,contig_labels,output, contig_dict , recluster = False,origin_label=0):
    from collections import defaultdict
    res = defaultdict(list)
    for label, name in zip(contig_labels, namelist):
        if label != -1:
            res[label].append(name)

    os.makedirs(output, exist_ok=True)

    for label in res:
        bin = []
        whole_bin_bp = 0
        for contig in res[label]:
            rec = SeqRecord(Seq(str(contig_dict[contig])), id=contig, description='')
            bin.append(rec)
            whole_bin_bp += len(str(contig_dict[contig]))
        if not recluster:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'bin.{}.fa'.format(label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')
        else:
            if whole_bin_bp >= 200000:
                with atomic_write(os.path.join(output, 'recluster_{0}.bin.{1}.fa'.format(origin_label,label)), overwrite=True) as ofile:
                    SeqIO.write(bin, ofile, 'fasta')


def cal_kl(m1,m2,v1,v2):
        m1 = max(m1,1e-6)
        m2 = max(m2,1e-6)
        v1 = 1 if v1 < 1 else v1
        v2 = 1 if v2 < 1 else v2
        value = np.log(np.sqrt(v2 / v1)) + np.divide(np.add(v1,np.square(m1 - m2)),2 * v2) - 0.5
        return min(max(value,1e-6),1-1e-6)

def cal_num_bins(fasta_path,contig_output,hmm_output,seed_output,binned_short):

    if not os.path.exists(contig_output + '.faa'):

        frag_out_log = open(contig_output + '.out','w')
        subprocess.check_call(
            ['run_FragGeneScan.pl',
             '-genome={}'.format(fasta_path),
             '-out={}'.format(contig_output),
             '-complete=0',
             '-train=complete',
             '-thread=48',
             ],
            stdout = frag_out_log,
            stderr = subprocess.DEVNULL,
        )


    if not os.path.exists(hmm_output):

        hmm_out_log = open(hmm_output+'.out','w')
        subprocess.check_call(
            ['hmmsearch',
             '--domtblout',
             hmm_output,
             '--cut_tc',
             '--cpu', str(48),
             'marker.hmm',
             contig_output+'.faa',
            ],
            stdout=hmm_out_log,
            stderr=subprocess.DEVNULL,
        )

    if not os.path.exists(seed_output):
            if binned_short:
                getmarker = 'test_getmarker.pl'
                subprocess.check_call(
                    ['perl', getmarker,
                     hmm_output,
                     fasta_path,
                     str(1001),
                     seed_output,
                     ],
                )

            else:
                getmarker =  'test_getmarker.pl'
                subprocess.check_call(
                    ['perl', getmarker,
                     hmm_output,
                     fasta_path,
                     str(2501),
                     seed_output,
                     ],
                )

def main(args=None):
    if args is None:
        args = sys.argv

    args = parse_args(args)
    validate_args(args)


    logger = logging.getLogger('SemiBin')
    logger.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    out = args.output
    if not os.path.exists(out):
        os.mkdir(out)



    whole_contig_bp = 0
    contig_bp_2500 = 0
    contig_length_list = []
    contig_length_dict = {}
    contig_dict = {}
    for seq_record in SeqIO.parse(args.contig_fasta , "fasta"):
        if len(seq_record) >= 1000 and len(seq_record) <= 2500:
            contig_bp_2500 += len(seq_record)
        contig_length_list.append(len(seq_record))
        whole_contig_bp += len(seq_record)
        contig_length_dict[str(seq_record.id).strip('')] = len((seq_record.seq))
        contig_dict[str(seq_record.id).strip('')] = str(seq_record.seq)

    # threshold for generating must link pairs
    threshold = get_threshold(contig_length_list)

    if contig_bp_2500 / whole_contig_bp >= 0.05:
        binned_short = False
    else:
        binned_short = True


    # generating coverage for every contig and for must link pair

    data = pd.read_csv(os.path.join(out, 'data.csv'), index_col=0)

    kmer = data.values[:,0:136]
    depth = data.values[:,136:len(data.values[0])]


    namelist = data.index.tolist()
    mapObj = dict(zip(namelist, range(len(namelist))))
    row_index = data._stat_axis.values.tolist()
    train_data = data.values

    n_sample = args.n_sample
    print(n_sample)
    is_combined = True if n_sample >= 5 else False

    if not is_combined:
        train_data_input = train_data[:,0:136]
    else:
        train_data_input = train_data

    #x = torch.from_numpy(train_data_input)
    print(train_data_input.shape)
    print(args.max_edges)
    embedding = train_data_input
    embedding_matrix = kneighbors_graph(embedding, n_neighbors=args.max_edges, mode='distance', p=2, n_jobs=-1).toarray()

    embedding_matrix[embedding_matrix >= 1] = 1
    embedding_matrix[embedding_matrix == 0] = 1
    embedding_matrix = 1 - embedding_matrix



    #cannot_list = pd.read_csv('/home1/pansj/binning/CAMI_medium/mmseqs_solidbin/medium_cannot.txt', sep=',',header=None).values
    #cannot_list = pd.read_csv('/share/inspurStorage/home1/pansj/binning/CAMI_low/mmseqs_solidbin/low_cannot.txt', sep=',',header=None).values

    if args.cannot_link is not None:
        print('can not link')
        cannot_list = pd.read_csv(args.cannot_link,sep=',',header=None).values
        print(len(cannot_list))
        for temp in cannot_list:
            embedding_matrix[mapObj[temp[0]]][mapObj[temp[1]]] = 0

    if args.must_link is not None:
        print('must link')
        must_list = pd.read_csv(args.must_link,sep=',',header=None).values
        print(len(must_list))
        for temp in must_list:
            embedding_matrix[mapObj[temp[0]]][mapObj[temp[1]]] = 1
    print(embedding_matrix.shape)
    threshold = 0.95

    while (threshold >= 0):
        num = len(list(set(np.where(embedding_matrix > threshold)[0])))
        if round(num / len(embedding_matrix), 2) >= args.max_node:
            break
        else:
            threshold -= 0.05

    embedding_matrix[embedding_matrix <= threshold] = 0
    if not is_combined:
        logger.info('Calculating depth matrix.')
        depth_matrix = np.zeros(shape=embedding_matrix.shape)
        for i in range(len(embedding_matrix)):
            for j in range(i + 1, len(embedding_matrix)):
                if embedding_matrix[i][j] > 0:
                        temp_depth = 0
                        for k in range(n_sample):
                            temp_depth  += 1 - cal_kl(depth[i][2*k], depth[j][2*k], depth[i][2*k+1], depth[j][2*k+1])
                        depth_matrix[i][j] = temp_depth / n_sample

        matrix = embedding_matrix * depth_matrix
    else:
        matrix = embedding_matrix

    edges = []
    edges_weight = []

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if matrix[i][j] > 1e-6:
                edges.append((i, j))
                edges_weight.append(matrix[i][j])

    logger.info('Edges:{}'.format(len(edges)))

    g = Graph()
    vertex = list(range(len(matrix)))
    g.add_vertices(vertex)
    g.add_edges(edges)
    length_weight = np.array([contig_length_dict[name] for name in namelist])
    result = g.community_infomap(edge_weights=edges_weight,vertex_weights=length_weight)
    contig_labels = np.zeros(shape=(len(matrix)), dtype=np.int)

    for i in range(len(result)):
        temp = result[i]
        for infomap_index in temp:
            contig_labels[infomap_index] = i

    output_bin_path = os.path.join(out,'output_bins')
    if not os.path.exists(output_bin_path):
        os.mkdir(output_bin_path)

    write_bins(namelist, contig_labels, output_bin_path, contig_dict)
    if not is_combined:
        mean_index = [2 * temp for temp in range(n_sample)]
        depth_mean = depth[:, mean_index] / 100
        scaling = np.mean(np.abs(embedding)) / np.mean(depth_mean)
        base = 10
        weight = 2 * base * math.ceil(scaling / base)
        embedding_new = np.concatenate((embedding, depth_mean * weight), axis=1)
    else:
        embedding_new = embedding

    bin_files = os.listdir(output_bin_path)
    logger.info('Reclustering.')

    for bin in bin_files:
        if os.path.exists(os.path.join(output_bin_path, bin)):
            contig_list = []
            for seq_record in SeqIO.parse(os.path.join(output_bin_path, bin), "fasta"):
                contig_list.append(seq_record.id)
            contig_output = os.path.join(output_bin_path, bin) + '.frag'
            hmm_output = os.path.join(output_bin_path, bin) + '.hmmout'
            seed_output = os.path.join(output_bin_path, bin) + '.seed'
            try:
                cal_num_bins(os.path.join(output_bin_path, bin),contig_output,hmm_output,seed_output,binned_short)
            except:
                pass
            contig_index = [mapObj[temp] for temp in contig_list]
            re_bin_features = embedding_new[contig_index]
            if not os.path.exists(os.path.join(out, 'output_recluster_bins')):
                os.mkdir(os.path.join(out, 'output_recluster_bins'))

            if os.path.exists(seed_output):
                seed = open(seed_output).read().split('\n')
                seed = [contig for contig in seed if contig != '']
                init_seed = seed
                num_bin = len(seed)
                seed_index = []
                for temp in init_seed:
                    seed_index.append(row_index.index(temp))
                length_weight = np.array([contig_length_dict[name] for name in contig_list])
                seeds_embedding = embedding_new[seed_index]
                kmeans = KMeans(n_clusters=num_bin, init=seeds_embedding,n_init=1)
                kmeans.fit(re_bin_features, sample_weight=length_weight)
                labels = kmeans.labels_
                write_bins(contig_list, labels, os.path.join(out, 'output_recluster_bins'), contig_dict,
                           recluster=True, origin_label=int(bin.split('.')[-2]))
            else:
                shutil.copy(os.path.join(output_bin_path, bin), os.path.join(out, 'output_recluster_bins'))

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main(sys.argv)

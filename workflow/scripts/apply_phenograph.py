import argparse
import logging
#from unsupervised.phenograph import Phenograph
import sys
import pip
import numpy as np
import os
import phenograph
from abc import ABCMeta, abstractmethod
import h5py
import pandas as pd

logging.basicConfig(level=logging.DEBUG)

logging.debug("system version : " + str(sys.version))

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", dest="input_file", required=True, help="input matrix hdf5 file")
parser.add_argument("--output_file", dest="output_file", required=True, help="path to the output csv file")
parser.add_argument("--distance_matrix", dest="distance_matrix", required=True, help="the distance matrix (output file) used by the louvain clustering")
parser.add_argument("--modularity_score", dest="modularity_score", required=True, help="text file containing the louvain modularity score")
parser.add_argument("--n_neighbours", dest="n_neighbours", required=True, help="the number of neighbours", type=int)
parser.add_argument("--min_size", dest="min_size", required=True, help="minimum cluster size", type=int)
parser.add_argument("--n_threads", dest="n_threads", required=True, help="the number of threads", type=int, default=1)
parser.add_argument("-l", "--log_normalize", dest="log_normalize", required=True, help="Boolean switch for the log normalization", type=str2bool)
args = parser.parse_args()


class UnsupervisedMethod(metaclass=ABCMeta):

    def __init__(self):
        self.matrix = None
        self.barcodes = None
        self.results = None

    def load_from_csv(self, input_file):
        df = pd.read_csv(input_file, header=None)
        self.barcodes = df[0].values
        df = df.drop(axis=1, columns=[0])
        self.matrix = df.values

    def load_from_hdf5(self, input_file):
        h5f = h5py.File(input_file, 'r')
        self.matrix = h5f['cor_counts'][()]
        barcodes = h5f['cell_attrs']['cell_names'][()]
        # decoder is needed since the input is a binary string
        decoder = np.vectorize(lambda t: t.decode('UTF-8'))
        self.barcodes = decoder(barcodes)
        h5f.close()


    def write_csv(self, output_file, dim_red_res=False):
        df = pd.DataFrame(self.barcodes)
        df = pd.concat([df, pd.DataFrame(self.results)], axis=1)
        logging.debug('set of communities to be written: ' + str(set(self.results)))
        if dim_red_res:
            df.to_csv(output_file,header=False, index=False)
        else:
            # final clustering results
            df.columns = ['cell_barcode','cluster']
            df.to_csv(output_file,header=True, index=False)

    @abstractmethod
    def apply(self):
        pass

    def log_normalize(self):
        self.matrix = np.log1p(self.matrix)




class Phenograph(UnsupervisedMethod):

    def __init__(self, n_neighbours, min_size, threads):
        super().__init__()
        self.n_neighbours = n_neighbours
        self.threads = threads
        self.distance_matrix = None
        self.modularity = None
        self.min_size = min_size

    def apply(self):
        communities, graph, Q = phenograph.cluster(data=self.matrix,k=self.n_neighbours, min_cluster_size=self.min_size, n_jobs=self.threads)
        # add 1 to the cluster labels to shift -1 values to zero.
        communities = communities + 1

        self.results = communities

        arr = graph.toarray()
        arr_full = arr+arr.T
        np.fill_diagonal(arr_full, 1)
        dist = (arr_full- arr_full.max())*(-1)
        np.fill_diagonal(dist, 0)


        self.distance_matrix = dist
        self.modularity = Q

        set_c = set(communities)
        logging.debug('set of communities found by phenograph: ' + str(set_c))

        # make sure at least 2 communities are found
        try:
            assert(len(set_c) > 1)
        except AssertionError as error:
            # Output expected AssertionErrors.
            logging.debug('Less than 2 communities is found')
            sys.exit(1)




def apply_pheno(input_file, output_file, distance_matrix, modularity_score, n_neighbours, min_size, threads, log_normalize):

    logging.debug('input file: ' + input_file)
    logging.debug('n_neighbours: ' +  str(n_neighbours))
    logging.debug('min_size:' + str(min_size))
    logging.debug('log_normalization value: ' + str(log_normalize))
    logging.debug('n_threads: ' +  str(threads))
    logging.debug('type of log_normalization value: ' + str(type(log_normalize)))

    pheno = Phenograph(n_neighbours, min_size, threads)
    pheno.load_from_hdf5(input_file)
    if log_normalize:
        pheno.log_normalize()
    else:
        pass
    pheno.apply()


    logging.debug(pheno.barcodes[1:5])
    logging.debug('distance matrix')
    logging.debug(pheno.distance_matrix)

    pheno.write_csv(output_file)

    np.savetxt(fname = distance_matrix, X=pheno.distance_matrix, delimiter='\t')

    # write the score file
    f = open(modularity_score, 'w')
    f.write(str(pheno.modularity))
    f.close()



apply_pheno(args.input_file, args.output_file, args.distance_matrix, args.modularity_score, args.n_neighbours, args.min_size, args.n_threads, args.log_normalize)

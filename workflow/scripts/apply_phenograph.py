import argparse
import logging
from unsupervised.phenograph import Phenograph
from utils import str2bool
import sys
import pip
import numpy as np

logging.basicConfig(level=logging.DEBUG)

logging.debug("system version : " + str(sys.version))
# for pip version < 10.0.0
#pip_res = pip.get_installed_distributions()
#phenograph_version = [i for i in pip_res if str(i).startswith('Pheno')]
#logging.debug("phenograph version: " + str(phenograph_version))

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

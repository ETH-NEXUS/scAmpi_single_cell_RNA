import phenograph
import logging
from .unsupervised_method import UnsupervisedMethod
import sys
import numpy as np

logging.basicConfig(level=logging.DEBUG)

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

from abc import ABCMeta, abstractmethod
import h5py
import numpy as np
import pandas as pd
import logging

logging.basicConfig(level=logging.DEBUG)

logging.debug('pandas version: ' + str(pd.__version__))
logging.debug('numpy version: ' + str(np.__version__))
logging.debug('h5py version: ' + str(h5py.__version__))

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

from ZIFA import ZIFA, block_ZIFA
from .unsupervised_method import UnsupervisedMethod

class Zifa(UnsupervisedMethod):

    def __init__(self, n_components, n_blocks):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.n_blocks = n_blocks

    def apply(self):
        Zhat, params = block_ZIFA.fitModel(self.matrix, self.n_components, n_blocks = self.n_blocks)
        self.results = Zhat

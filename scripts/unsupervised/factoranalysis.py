from sklearn.decomposition import FactorAnalysis as sk_fa
from .unsupervised_method import UnsupervisedMethod

class FactorAnalysis(UnsupervisedMethod):

    def __init__(self, n_components):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components

    def apply(self):
        fa = sk_fa(n_components = self.n_components)
        self.results = fa.fit_transform(self.matrix)

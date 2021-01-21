from sklearn.decomposition import PCA
from .unsupervised_method import UnsupervisedMethod

class Pca(UnsupervisedMethod):

    def __init__(self, n_components):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components

    def apply(self):
        pca = PCA(n_components = self.n_components)
        self.results = pca.fit_transform(self.matrix)

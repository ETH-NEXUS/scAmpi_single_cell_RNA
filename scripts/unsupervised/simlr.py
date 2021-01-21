import SIMLR
from .unsupervised_method import UnsupervisedMethod

class Simlr(UnsupervisedMethod):

    def __init__(self,n_components, pca_components, n_neighbours, max_iter, threads):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.pca_components = pca_components
        self.n_neighbours = n_neighbours
        self.max_iter = max_iter
        self.threads = threads

    def set_params(self, n_clusters):
        self.n_components = n_clusters
        return self

    def fit_predict(self, matrix):
        X = matrix

        # Selecting 500 features with PCA
        if X.shape[1]>500:
        # fast_pca assumes the number of cells > 500 therefore try-catch
            try:
                X = SIMLR.helper.fast_pca(X,self.pca_components)
            except:
                pass
        # Running Simlr
        simlr = SIMLR.SIMLR_LARGE(num_of_rank=self.n_components, num_of_neighbor=self.n_neighbours, max_iter=self.max_iter)
        S, F, val, ind = simlr.fit(X)
        return F

    def apply(self):
        self.results = fit_predict(self.matrix)

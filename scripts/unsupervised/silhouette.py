from .unsupervised_method import UnsupervisedMethod
from .simlr import Simlr
from sklearn.metrics.cluster import silhouette_score
from sklearn.cluster import KMeans, AgglomerativeClustering
import numpy as np

class Silhouette(UnsupervisedMethod):

    def __init__(self, model, k_min, k_max, metric):
        UnsupervisedMethod.__init__(self)
        self.model = model
        self.k_min = k_min
        self.k_max = k_max
        self.metric = metric


    @classmethod
    def init_with_kmeans(cls, n_init, n_jobs, k_min, k_max, metric):

        # K-means falls in local minima. That’s why it can be useful to restart it several times using n_init
        # The classical EM-style algorithm is "full" and it is suggested for sparse data

        model = KMeans(algorithm='full', n_init=n_init, n_jobs=n_jobs)
        return cls(model, k_min, k_max, metric)

    @classmethod
    def init_with_hierarchical(cls, h_affinity, h_linkage, k_min, k_max, metric):
        model = AgglomerativeClustering(affinity=h_affinity, linkage=h_linkage, compute_full_tree=False)
        return cls(model, k_min, k_max, metric)

    @classmethod
    def init_with_simlr(cls, n_components, pca_components, n_neighbours, max_iter, threads, k_min, k_max, metric):
        model = Simlr(n_components, pca_components, n_neighbours, max_iter, threads)
        return cls(model, k_min, k_max, metric)

    def apply(self):
        k_range = range(self.k_min, self.k_max)
        if isinstance(self.model, Simlr):
        # uses generative object and np.fromiter for memory efficiency
            go = (silhouette_score(X=self.matrix, labels=KMeans(k).fit_predict(self.model.set_params(n_clusters=k).fit_predict(self.matrix))) for k in k_range)
            silhouette_scores = np.fromiter(go, dtype=float, count=len(k_range))
        else:
            predicted_labels = [self.model.set_params(n_clusters=k).fit_predict(self.matrix) for k in k_range]
            silhouette_scores = [silhouette_score(X=self.matrix, labels=obj, metric=self.metric) for obj in predicted_labels]

        max_index = np.argmax(silhouette_scores)
        self.results = predicted_labels[max_index]

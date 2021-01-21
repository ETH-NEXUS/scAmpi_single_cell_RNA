from sklearn.manifold import TSNE
from .unsupervised_method import UnsupervisedMethod

class Tsne(UnsupervisedMethod):

    def __init__(self, n_components, init):
        UnsupervisedMethod.__init__(self)
        self.n_components = n_components
        self.init = init

    def apply(self):
        tsne = TSNE(init=self.init, n_components=self.n_components)
        self.results = tsne.fit_transform(self.matrix)

import numpy as np
import networkx as nx
from sklearn.neighbors import NearestNeighbors


class SNNGraph:
    def __init__(self, data_matrix, k=9, round_weights=False, decimal_places=1):

        self.data_matrix = data_matrix
        self.k = k
        self.round_weights = round_weights
        self.decimal_places = decimal_places
        self.graph = None

    def compute_nearest_neighbors(self):
        nn = NearestNeighbors(n_neighbors=self.k + 1, metric='euclidean')
        nn.fit(self.data_matrix)
        distances, indices = nn.kneighbors(self.data_matrix)
        return indices

    def construct_graph(self, indices):
        G = nx.Graph()
        num_nodes = self.data_matrix.shape[0]

        G.add_nodes_from(range(num_nodes))

        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):
                shared_neighbors = np.intersect1d(indices[i], indices[j])
                if len(shared_neighbors) > 0:
                    rank_i = np.mean([np.where(indices[i] == neighbor)[0][0] for neighbor in shared_neighbors])
                    rank_j = np.mean([np.where(indices[j] == neighbor)[0][0] for neighbor in shared_neighbors])
                    weight = self.k - 0.5 * (rank_i + rank_j)
                    if self.round_weights:
                        weight = round(weight, self.decimal_places)
                    G.add_edge(i, j, weight=weight)

        self.graph = G

    def get_graph(self):
        if self.graph is None:
            indices = self.compute_nearest_neighbors()
            self.construct_graph(indices)
        return self.graph

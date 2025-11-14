import scipy.sparse as sp
import scanpy as sc
import networkx as nx
from networkx import Graph
import sys
class SNNClusterRunner:
    def __init__(self, graph: Graph, adata, weight_key="weight", decimal_places=None):
        self.graph = graph
        self.adata = adata
        self.weight_key = weight_key
        self.decimal_places = decimal_places
        self.adj_matrix = None

    def _to_sparse_matrix(self):
        adj = nx.to_scipy_sparse_array(self.graph, weight=self.weight_key, format='coo')
        if self.decimal_places is not None:
            adj.data = adj.data.round(self.decimal_places)
        self.adj_matrix = adj.tocsr()

    def inject_graph(self):

        if self.adj_matrix is None:
            self._to_sparse_matrix()
        self.adata.uns['neighbors'] = {
            'connectivities': self.adj_matrix,
            'params': {'method': 'custom_snn', 'n_neighbors': 'custom'}
        }
        self.adata.obsp['connectivities'] = self.adj_matrix

    def run_leiden(self, resolution=1.0):
        self.inject_graph()
        sc.tl.leiden(self.adata, resolution=resolution)
    def run_umap(self):
        sc.tl.umap(self.adata)

    def plot(self, color='leiden'):
        sc.pl.umap(self.adata, color=color)
import numpy as np
import networkx as nx
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

def create_cluster_graph(data, labels):
    labels = np.array(labels)

    unique_labels = np.unique(labels)

    cluster_centers = {}
    for label in unique_labels:
        cluster_data = data[labels == label]
        cluster_center = np.mean(cluster_data, axis=0)
        cluster_centers[label] = cluster_center

    centers_list = [cluster_centers[label] for label in unique_labels]
    distances = pairwise_distances(centers_list)

    G = nx.Graph()

    for label in unique_labels:
        G.add_node(label, name=f"Cluster {label}", pos=cluster_centers[label])

    for i, label1 in enumerate(unique_labels):
        for j, label2 in enumerate(unique_labels):
            if i < j:
                G.add_edge(label1, label2, weight=distances[i, j])

    return G

if __name__ == "__main__":
    np.random.seed(0)
    data = np.random.rand(10, 2)
    labels = [0, 0, 0, 1, 1, 1, 2, 2, 3, 3]

    G = create_cluster_graph(data, labels)

    print("节点:", G.nodes(data=True))
    print("边:", G.edges(data=True))

    import matplotlib.pyplot as plt

    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', font_size=12)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    plt.show()

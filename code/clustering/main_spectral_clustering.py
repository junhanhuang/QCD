from preprocess import process_data
from snn_graph import SNNGraph
from ClusterGraph import create_cluster_graph
from draw import plot_umap_clusters_with_mst
import networkx as nx
import numpy as np
import os
import pandas as pd
from sklearn.cluster import SpectralClustering

def run_pipeline_spectral(file_path, output_dir, params):
    print(f"开始处理数据集: {file_path}")

    processed_df = process_data(file_path, top_n=params['top_n'])
    os.makedirs(output_dir, exist_ok=True)
    processed_df.to_csv(os.path.join(output_dir, "norm_matrix.csv"))

    snn = SNNGraph(data_matrix=processed_df, k=params['k'], decimal_places=params['decimal_places'])
    G_nx = snn.get_graph()

    nodes = list(G_nx.nodes)
    adj_matrix = nx.to_scipy_sparse_array(G_nx, nodelist=nodes, format="csr")
    adj_matrix.indptr = adj_matrix.indptr.astype(np.int32)
    adj_matrix.indices = adj_matrix.indices.astype(np.int32)
    n_clusters = params.get("n_clusters", 5)
    clustering = SpectralClustering(
        n_clusters=n_clusters,
        affinity='precomputed',
        assign_labels='kmeans',
        random_state=0
    )
    spectral_labels = clustering.fit_predict(adj_matrix)

    cluster_df = pd.DataFrame({
        'cell_id': processed_df.index,
        'group_id': [f'Path{label}' for label in spectral_labels]
    })
    cluster_df.to_csv(os.path.join(output_dir, "spectral_cluster.csv"), index=False)

    G_cluster = create_cluster_graph(processed_df, spectral_labels)
    for u, v in G_cluster.edges():
        weight = np.linalg.norm(np.array(G_cluster.nodes[u]['pos']) - np.array(G_cluster.nodes[v]['pos']))
        G_cluster[u][v]['weight'] = weight
    mst_nx = nx.minimum_spanning_tree(G_cluster, weight='weight')
    mst_edges = list(mst_nx.edges())

    embedding = plot_umap_clusters_with_mst(
        data=processed_df,
        labels=spectral_labels,
        mst_edges=mst_edges,
        folder_name=output_dir,
        file_name="spectral_result.png"
    )
    np.savetxt(os.path.join(output_dir, "umap_coords.csv"), embedding, delimiter=",")
    np.savetxt(os.path.join(output_dir, "spectral_labels.csv"), np.array(spectral_labels), fmt='%d', delimiter=",")

    print(f" 数据集处理完成: {file_path}")

if __name__ == "__main__":
    import yaml
    root_dir = os.path.dirname(os.path.dirname(__file__))
    config_path = os.path.join(root_dir, 'src', 'config_spectral.yaml')

    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    default_params = config['default']

    for data_type, datasets in config['datasets'].items():
        for file_path, specific_params in datasets.items():
            params = default_params.copy()
            params.update(specific_params)
            full_file_path = os.path.join(root_dir, 'data', file_path)
            output_dir = os.path.join(root_dir, 'results', os.path.dirname(file_path))
            run_pipeline_spectral(full_file_path, output_dir, params)

    print(" 所有数据集处理完成！")
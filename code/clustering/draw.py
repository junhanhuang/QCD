import numpy as np
import matplotlib.pyplot as plt
import umap
import seaborn as sns
import os

def plot_umap_clusters_with_mst(data, labels, mst_edges, folder_name, file_name,
                                random_state=42, figsize=(12, 8), point_size=100, alpha=0.8):


    reducer = umap.UMAP(random_state=random_state)
    embedding = reducer.fit_transform(data)

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    label_to_index = {label: idx for idx, label in enumerate(unique_labels)}

    centroids = np.array([
        embedding[np.array(labels) == label].mean(axis=0) for label in unique_labels
    ])

    fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(
        x=embedding[:, 0],
        y=embedding[:, 1],
        hue=labels,
        palette='tab10',
        s=point_size,
        alpha=alpha,
        legend='full',
        ax=ax
    )

    for i, j in mst_edges:
        idx_i, idx_j = label_to_index[i], label_to_index[j]
        xi, yi = centroids[idx_i]
        xj, yj = centroids[idx_j]

        ax.plot([xi, xj], [yi, yj], 'k-', alpha=0.7, linewidth=3)
        ax.text(xi, yi, str(i), fontsize=14, ha='left', va='bottom', weight='bold', color='black')
        ax.text(xj, yj, str(j), fontsize=14, ha='left', va='bottom', weight='bold', color='black')

    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.legend(title='State', title_fontsize='13', fontsize='11')
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    folder_path = os.path.join('result', folder_name)
    os.makedirs(folder_path, exist_ok=True)
    filepath = os.path.join(folder_path, file_name)

    plt.savefig(filepath, format='png', bbox_inches='tight', dpi=300)
    plt.show()

    return embedding

import pandas as pd
import numpy as np
from scipy.stats import zscore
import os

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(ROOT_DIR, 'result')


def save_results(df, filename):
    if not os.path.exists(RESULT_DIR):
        os.makedirs(RESULT_DIR)


def load_data(file_path):
    return pd.read_csv(file_path, index_col=0)


def filter_genes(df: pd.DataFrame):
    non_zero_counts = (df != 0).sum(axis=0)
    threshold = 0.05 * df.shape[0]
    filtered_genes = df.columns[non_zero_counts >= threshold]
    return df[filtered_genes], df.shape[1], len(filtered_genes), df.shape[1] - len(filtered_genes)


def log_transform(df):
    return np.log2(df + 1)


def standardize_data(df):
    return df.apply(zscore)


def sort_genes_by_variance(df):
    return df[df.var(axis=0).sort_values(ascending=False).index]


def select_top_n_genes(df, n):
    return df.iloc[:, :n]


def process_data(file_path, top_n):
    df = load_data(file_path)
    # df = df.T
    original_gene_count = df.shape[1]

    df, original_count, remaining_count, removed_count = filter_genes(df)
    print(f"原始共有基因数量: {original_count}")
    print(f"去除后的基因数量: {remaining_count}")
    print(f"去除的基因数量: {removed_count}")

    df = sort_genes_by_variance(df)

    df = log_transform(df)

    df = standardize_data(df)

    df = select_top_n_genes(df, top_n)

    save_results(df, f'processed_data_top_{top_n}_genes.csv')

    return df

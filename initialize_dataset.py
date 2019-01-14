import os
from scRNA_seq import Gene_Expression_Dataset

seed_matrices_file_paths = []

seed_matrices_file_paths.append(os.path.join(os.getcwd(), "transcripts", "gene_counts.csv"))

dataset_path = os.path.join(os.getcwd(), "workspace")

dataset = Gene_Expression_Dataset(
    dataset_path,
    seed_matrix_file_path=seed_matrices_file_paths
)

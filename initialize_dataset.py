import os
from scRNA_seq import Gene_Expression_Dataset

seed_matrices_file_paths = []

seed_matrices_file_paths.append(os.path.expanduser(
    os.path.join(
        "~", "Nf-1_Transcriptomics", "TENX-70-NUCSEQ-SMALL", "mm10_premrna",
        "matrix.csv"
    )
))

# seed_matrices_file_paths.append(os.path.expanduser(
#     os.path.join(
#         "~", "Virus_Transcriptomics", "C2-small", "matrixdata",
#         "C2-small_matrix.csv"
#     )
# ))

dataset_path = os.path.expanduser(
    os.path.join(
        "~", "Nf-1_Transcriptomics", "TENX-70-NUCSEQ-SMALL", "workspace"
    )
)

dataset = Gene_Expression_Dataset(
    dataset_path,
    seed_matrix_file_path=seed_matrices_file_paths
)

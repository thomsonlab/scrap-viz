import os
from scRNA_seq import Gene_Expression_Dataset

seed_matrices_file_paths = []

seed_matrices_file_paths.append("/media/dibidave/mobinoodle/Virus_Transcriptomics/DRG/transcripts/gene_counts.csv")

# seed_matrices_file_paths.append(os.path.expanduser(
#     os.path.join(
#         "~", "Virus_Transcriptomics", "C2-small", "matrixdata",
#         "C2-small_matrix.csv"
#     )
# ))

dataset_path = "/media/dibidave/mobinoodle/Virus_Transcriptomics/DRG/workspace"

dataset = Gene_Expression_Dataset(
    dataset_path,
    seed_matrix_file_path=seed_matrices_file_paths
)

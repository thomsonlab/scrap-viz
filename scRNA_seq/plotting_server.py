import os
from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot

dataset_path = os.path.expanduser(
    os.path.join(
        "~", "Virus_Transcriptomics", "C2-small", "workspace"
    )
)

gene_expression_dataset = Gene_Expression_Dataset(dataset_path)


gene_expression_dataset.filter_low_gene_counts(3)
gene_expression_dataset.filter_low_transcript_cells(1000)
gene_expression_dataset.normalize(
    Gene_Expression_Dataset.Normalization_Method.STD)
# gene_expression_dataset.transform(
#     Gene_Expression_Dataset.Transformation_Method.PCA, num_dimensions=10)
# gene_expression_dataset.transform(
#     Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2)

plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

plot.start()

import os
from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot

dataset_path = os.path.expanduser(
    os.path.join(
        "~", "Virus_Transcriptomics", "C2-small", "workspace"
    )
)

# dataset_path = os.path.expanduser(
#     os.path.join(
#         "~", "Nf-1_Transcriptomics", "workspace"
#     )
# )

print("Loading dataset...")
gene_expression_dataset = Gene_Expression_Dataset(dataset_path)

print("Filtering...")
gene_expression_dataset.filter_low_gene_counts(5)
gene_expression_dataset.filter_low_transcript_cells(1000)
print("Normalizing...")
gene_expression_dataset.normalize_cells(
    Gene_Expression_Dataset.Data_Mode.READS_PER_MILLION_TRANSCRIPTS)
gene_expression_dataset.normalize_genes(
    Gene_Expression_Dataset.Normalization_Method.STD)
print("Transforming...")
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.PCA, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2,
    use_normalized=True)
print("Ready to plot!")

plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

plot.start()


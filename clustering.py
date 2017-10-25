import os
from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot

Nf1_dataset_path = os.path.expanduser(
    os.path.join(
        "~", "NF-1_transcriptomics", "analysis"
    )
)

gene_expression_dataset = Gene_Expression_Dataset(Nf1_dataset_path)
gene_expression_dataset.filter_low_gene_counts(5)
gene_expression_dataset.filter_low_transcript_cells(1500)
gene_expression_dataset.normalize(
    Gene_Expression_Dataset.Normalization_Method.STD)
# gene_expression_dataset.transform(
#     Gene_Expression_Dataset.Transformation_Method.PCA, num_dimensions=10)
# gene_expression_dataset.transform(
#     Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2)

plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

plot.plot_tSNE()


# plotting.plot_tSNE(transformed_tSNE)
#
# sample_means = cdf_dataset.get_sample_means()
# sample_means = gene_expression_dataset.get_sample_means()
#
# plotting.plot_differential_expression(sample_means)

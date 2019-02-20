import os
from scrap_viz import Gene_Expression_Dataset
import sys


dataset_path = os.getcwd()

pipeline_name = "5_RPM_SQRT"

gene_expression_dataset = Gene_Expression_Dataset(dataset_path)

print("Filtering...")
gene_expression_dataset.filter_low_gene_counts(5)

print("Normalizing...")
gene_expression_dataset.normalize_cells(
    Gene_Expression_Dataset.Data_Mode.READS_PER_MILLION_TRANSCRIPTS)
gene_expression_dataset.normalize_genes(
    Gene_Expression_Dataset.Normalization_Method.SQUARE_ROOT)

print("Transforming...")
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.PCA, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.NMF, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.SVD, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2,
    use_normalized=True)

gene_expression_dataset.save(pipeline_name)

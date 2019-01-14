import os
from scRNA_seq import Gene_Expression_Dataset
import sys


dataset_path = os.getcwd()

pipeline_name = "3_800_RPM_SD"

gene_expression_dataset = Gene_Expression_Dataset(dataset_path)

print("Filtering...")
gene_expression_dataset.filter_low_gene_counts(3)
gene_expression_dataset.filter_low_transcript_cells(800)

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
    Gene_Expression_Dataset.Transformation_Method.NMF, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.SVD, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2,
    use_normalized=True)

gene_expression_dataset.save(pipeline_name)

import os
import sys

import scrap_viz
from scrap_viz import Gene_Expression_Dataset
from scrap_viz import Data_Mode
from scrap_viz import Transformation_Method
from scrap_viz import Normalization_Method


dataset_path = os.getcwd()

pipeline_name = "5_RPM_SQRT"

gene_expression_dataset = Gene_Expression_Dataset(dataset_path)

print("Filtering...")
gene_expression_dataset.filter_low_gene_counts(5)

print("Normalizing...")
gene_expression_dataset.normalize_cells(Data_Mode.READS_PER_MILLION_TRANSCRIPTS)
gene_expression_dataset.normalize_genes(Normalization_Method.SQUARE_ROOT,
    use_normalized=True)

print("Transforming...")
gene_expression_dataset.transform(
    Transformation_Method.PCA, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Transformation_Method.NMF, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Transformation_Method.SVD, num_dimensions=30,
    use_normalized=True)
gene_expression_dataset.transform(
    Transformation_Method.TSNE, num_dimensions=2,
    use_normalized=True)

gene_expression_dataset.save(pipeline_name)

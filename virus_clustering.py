import os
from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot
import pandas

virus_dataset_path = os.path.expanduser(
    os.path.join(
        "~", "Virus_Transcriptomics", "C2-small", "workspace"
    )
)

virus_count_file_path = os.path.expanduser(
    os.path.join(
        "~", "GradinaruLab", "AAVs", "workspace", "drop_seq_1", "export",
        "C2-small_virus_matrix.csv"
    )
)

gene_expression_dataset = Gene_Expression_Dataset(virus_dataset_path)

virus_count_data_frame = pandas.read_csv(
    virus_count_file_path, sep=",", header=0, index_col=0)

for virus_name, cell_virus_counts in virus_count_data_frame.iterrows():
    virus_cells = set(cell_virus_counts[cell_virus_counts > 1].index)
    gene_expression_dataset.label_cells(virus_name, virus_cells)

gene_expression_dataset.filter_low_gene_counts(3)
gene_expression_dataset.filter_low_transcript_cells(1000)
gene_expression_dataset.normalize(
    Gene_Expression_Dataset.Normalization_Method.STD)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.PCA, num_dimensions=10)
gene_expression_dataset.transform(
    Gene_Expression_Dataset.Transformation_Method.TSNE, num_dimensions=2)

plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

plot.plot_tSNE()

diffs = gene_expression_dataset.compare_gene_expression(["Neurons", "PHP.eB"],
                                                        ["Neurons", "PHP.S"])

neurons_PHP_eB = gene_expression_dataset._label_cells["Neurons"]\
    .intersection(gene_expression_dataset._label_cells["PHP.eB"])
neurons_PHP_S = gene_expression_dataset._label_cells["Neurons"]\
    .intersection(gene_expression_dataset._label_cells["PHP.S"])

uninfected_cells = set(gene_expression_dataset._data_frame.columns)\
    .difference(gene_expression_dataset._label_cells["PHP.eB"])\
    .difference(gene_expression_dataset._label_cells["PHP.S"])

gene_expression_dataset.label_cells("Uninfected", uninfected_cells)

diffs = gene_expression_dataset.compare_gene_expression(["Neurons", "PHP.eB"],
                                                        ["Neurons", "PHP.S"])

cells2 = gene_expression_dataset._data_frame[
    list(gene_expression_dataset._label_cells["PHP.eB"].
         intersection(gene_expression_dataset._label_cells["Neurons"]))]

i = 0
for key in sorted(diffs.keys(), key=lambda x: diffs[x][1], reverse=False):
    print(key)
    print(diffs[key])
    i += 1
    if i > 20:
        break

# plotting.plot_tSNE(transformed_tSNE)
#
# sample_means = cdf_dataset.get_sample_means()
# sample_means = gene_expression_dataset.get_sample_means()
#
# plotting.plot_differential_expression(sample_means)

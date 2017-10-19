import os
import plotly
from plotly.graph_objs import Scatter
from plotly.graph_objs import Layout
from plotly.graph_objs import Figure
from scRNA_Seq import Gene_Expression_Dataset

gene_expression_matrix_path = os.path.expanduser(
    os.path.join(
        "~", "NF1 scRNA-Seq", "WT_KO_matrix_gene_threshold_5.csv"
    )
)

gene_expression_matrix_eCDF_path = os.path.expanduser(
    os.path.join(
        "~", "NF1 scRNA-Seq", "WT_KO_matrix_gene_threshold_5_eCDF.csv"
    )
)

gene_expression_dataset = Gene_Expression_Dataset(gene_expression_matrix_path)
cdf_dataset = Gene_Expression_Dataset(gene_expression_matrix_eCDF_path)

sample_means = cdf_dataset.get_sample_means()
sample_means = gene_expression_dataset.get_sample_means()

x_values = []
y_values = []
gene_names = []

for gene, gene_means in sample_means.iteritems():

    x_values.append(gene_means.iloc[0])
    y_values.append(gene_means.iloc[1])
    gene_names.append(gene)

min_value = min(min(x_values), min(y_values))
max_value = max(max(x_values), max(y_values))

gene_counts_scatter = Scatter(
    x=x_values,
    y=y_values,
    mode="markers",
    text=gene_names,
    name="Gene Expression")

x_y_line = Scatter(
    x=[min_value, max_value],
    y=[min_value, max_value],
    name="X=Y"
)

data = [gene_counts_scatter, x_y_line]
layout = Layout(
    xaxis=dict(
        range=[min_value, max_value],
        title="%s Gene Expression" % sample_means.iloc[0].name
    ),
    yaxis=dict(
        range=[min_value, max_value],
        title="%s Gene Expression" % sample_means.iloc[1].name
    ),
    hovermode="closest"
)

figure = Figure(data=data, layout=layout)

filename = "differential_expression_%i_mgc_%i_mtc.html" % \
           (5, 0)

plotly.offline.plot(figure, filename=filename)

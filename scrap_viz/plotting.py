import plotly
from plotly.graph_objs import Scatter
from plotly.graph_objs import Layout
from plotly.graph_objs import Figure
import pandas


def plot_eCDF(eCDF_values):

    eCDF_values = pandas.Series(eCDF_values)

    counts = eCDF_values.value_counts().sort_index()
    cum_sum = counts.cumsum()

    eCDF_values = cum_sum / sum(counts.values)

    eCDF_scatter = Scatter(
        x=eCDF_values.index,
        y=eCDF_values.values)

    data = [eCDF_scatter]

    layout = Layout(
        xaxis=dict(
            title="Value"
        ),
        yaxis=dict(
            title="Cumulative Probability"
        ),
        hovermode="closest"
    )

    figure = Figure(data=data, layout=layout)

    plotly.offline.plot(figure)


def plot_scatter(x, y):

    scatter = Scatter(x=x, y=y)

    data = [scatter]

    layout = Layout(hovermode="closest")

    figure = Figure(data=data, layout=layout)

    plotly.offline.plot(figure)


def plot_differential_expression(sample_means):

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

    filename = "differential_expression_%s_vs_%s_mtc.html" % \
               (sample_means.iloc[0].name, sample_means.iloc[1].name)

    plotly.offline.plot(figure, filename=filename)


def plot_tSNE(tSNE_transformed):

    x_values = tSNE_transformed.T.iloc[0]
    y_values = tSNE_transformed.T.iloc[1]

    tSNE_scatter = Scatter(
        x=x_values,
        y=y_values,
        mode="markers",
        name="tSNE")

    plot_data_sets = [tSNE_scatter]

    layout = Layout(
        xaxis=dict(
            title="tSNE 1"
        ),
        yaxis=dict(
            title="tSNE 2"
        ),
        hovermode="closest"
    )

    figure = Figure(data=plot_data_sets, layout=layout)

    filename = "tSNE.html"

    plotly.offline.plot(figure, filename=filename)

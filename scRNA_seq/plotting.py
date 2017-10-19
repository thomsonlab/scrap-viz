import plotly
from plotly.graph_objs import Scatter
from plotly.graph_objs import Layout
from plotly.graph_objs import Figure


def plot_eCDF(eCDF_values):

    eCDF_scatter = Scatter(
        x=eCDF_values.index,
        y=eCDF_values.values)

    data = [eCDF_scatter]

    layout = Layout(
        xaxis=dict(
            title="Transcript Count"
        ),
        yaxis=dict(
            title="Cumulative Probability"
        ),
        hovermode="closest"
    )

    figure = Figure(data=data, layout=layout)

    plotly.offline.plot(figure)
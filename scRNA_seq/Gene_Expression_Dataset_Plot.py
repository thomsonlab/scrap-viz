import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from .Gene_Expression_Dataset import Gene_Expression_Dataset

class Gene_Expression_Dataset_Plot:

    def __init__(self, gene_expression_dataset):

        self._app = None
        self._gene_expression_dataset = gene_expression_dataset
        self._dfs = {}

    def get_tSNE_figure(self):

        figure = {
            'data': [
                go.Scatter(
                    x=df['tSNE_1'],
                    y=df['tSNE_2'],
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size': 5,
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                    name=name,
                    text=df.index
                ) for name, df in self._dfs
                ],
            'layout': go.Layout(
                xaxis={'title': 'tSNE 1'},
                yaxis={'title': 'tSNE 2'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }

        return figure

    def plot_tSNE(self):

        self._dfs = self._gene_expression_dataset.get_cell_gene_expression_by_label(
            Gene_Expression_Dataset.Transformation_Method.TSNE)

        self._dfs = [(k, v) for k, v in sorted(self._dfs.items())]

        self._app = dash.Dash()

        self._tSNE_figure = self.get_tSNE_figure()

        self._tab = [

                dcc.Graph(
                    id='tSNE',
                    figure=self._tSNE_figure
                ),

                dcc.Input(id="label_name", type="text", value=""),

                html.Button("Add Label", id="add_label_button"),

                html.Div(className="row", children=[

                    html.Div([
                        html.Pre(id='selected-data')
                    ]),
                    html.Div([
                        html.Pre(id='trash')
                    ])
                ])
            ]

        self._app.layout = html.Div([

            html.Button("Clustering", id="clustering_button"),
            html.Button("Differential Expression",
                        id="differential_expression_button"),

            html.Div(id="tab", children=self._tab)
        ])

        @self._app.callback(
            dash.dependencies.Output("tab", "children"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tab_clicked(clustering_button_clicks, differential_expression_clicks):
            print(clustering_button_clicks)
            print(differential_expression_clicks)
            return self._tab

        @self._app.callback(
            dash.dependencies.Output("tSNE", "figure"),
            [dash.dependencies.Input("add_label_button", "n_clicks")],
            [dash.dependencies.State("label_name", "value"),
             dash.dependencies.State("tSNE", "selectedData")])
        def data_edited(_, label_name, selected_data):

            if selected_data is None:
                self._selected_points = None
                return self.get_tSNE_figure()
            else:
                self._selected_points = selected_data["points"]

            curve_number = 0

            for name, data_frame in self._dfs:
                point_indices = [
                    point["pointNumber"] for point in
                    self._selected_points if point["curveNumber"] ==
                    curve_number]
                cells = data_frame.index[point_indices]

                self._gene_expression_dataset.label_cells(label_name, cells)

                curve_number += 1

            return self.get_tSNE_figure()

        self._app.run_server()

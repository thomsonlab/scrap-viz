import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from .Gene_Expression_Dataset import Gene_Expression_Dataset
import plotly.figure_factory as ff
from scRNA_seq import Gene
import pandas


class Gene_Expression_Dataset_Plot:

    def __init__(self, gene_expression_dataset):

        self._app = None
        self._gene_expression_dataset = gene_expression_dataset
        self._dfs = {}
        self._differential_expression_n_clicks = 0
        self._clustering_n_clicks = 0

    def get_tSNE_figure(self):

        self._dfs = self._gene_expression_dataset.get_cell_gene_expression_by_label(
            Gene_Expression_Dataset.Transformation_Method.TSNE)

        self._dfs = [(k, v) for k, v in sorted(self._dfs.items())]

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

    @staticmethod
    def split_string(string):

        new_string = ""
        last_line_break_index = 0

        for character_index in range(len(string)):
            # TODO: This should not be a magic number
            if character_index - last_line_break_index > 50:
                if string[character_index] == " ":
                    new_string += string[last_line_break_index:character_index]
                    new_string += "<BR>"
                    last_line_break_index = character_index + 1

        new_string += string[last_line_break_index:]

        return new_string

    @staticmethod
    def generate_table(data_frame):

        genes = data_frame.index

        gene_name_descriptions = Gene.get_gene_summaries(genes)

        num_rows = data_frame.shape[0]

        num_columns = data_frame.shape[1]

        hover_texts = [[""] * (num_columns + 1)]

        for row in range(num_rows):

            gene_description = gene_name_descriptions[row][1]
            gene_description = Gene_Expression_Dataset_Plot.split_string(
                gene_description)
            hover_text = "<B>%s</B><BR>%s" %\
                (gene_name_descriptions[row][0], gene_description)

            hover_texts.append([hover_text] + [""] * num_columns)

        return dcc.Graph(
            id="de_figure", figure=ff.create_table(
                data_frame, index=True, index_title="Gene", hoverinfo="text",
                text=hover_texts),
            config={'displayModeBar': False})

    def _get_label_dropdowns(self):

        labels = self._gene_expression_dataset.get_labels()

        label_options = [{"label": x, "value": x} for x in labels]

        label_dropdowns = [
            html.Label('Subgroups 1'),
            dcc.Dropdown(
                id="subgroup_1_dropdown",
                options=label_options,
                value=[],
                multi=True
            ),
            html.Label('Subgroups 2'),
            dcc.Dropdown(
                id="subgroup_2_dropdown",
                options=label_options,
                value=[],
                multi=True
            )
        ]

        return label_dropdowns

    def _get_differential_expression_tab(self):

        label_dropdowns = self._get_label_dropdowns()

        differential_expression_tab = \
            html.Div(id="de_tab", children=[
                html.Div(
                    id="label_dropdowns",
                    children=self._get_label_dropdowns()),
                html.Button("Go!", id="de_button"),
                html.H4(children='Differential Expression'),
                html.Div(id="de_table", children=[])
            ])


        return differential_expression_tab

    def _get_tSNE_tab(self):

        tSNE_tab = \
            html.Div(id="tSNE_tab", children=[

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
            ])])

        return tSNE_tab

    def start(self):

        self._app = dash.Dash()

        self._tSNE_figure = self.get_tSNE_figure()

        self._tabs = [
            self._get_tSNE_tab(),
            self._get_differential_expression_tab()
        ]

        self._app.layout = html.Div([

            html.Button("Clustering", id="clustering_button"),
            html.Button("Differential Expression",
                        id="differential_expression_button"),

            html.Div(id="tabs", children=self._tabs)
        ])

        @self._app.callback(
            dash.dependencies.Output("tSNE_tab", "style"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tabs_clicked_update_tSNE_tab(clustering_button_clicks, _):

            if clustering_button_clicks is not None and \
                clustering_button_clicks > self._clustering_n_clicks:

                self._clustering_n_clicks += 1
                return {"display": "block"}

            else:
                return {"display": "none"}

        @self._app.callback(
            dash.dependencies.Output("label_dropdowns", "children"),
            [dash.dependencies.Input("add_label_button", "n_clicks")])
        def label_added_update_label_dropdowns(n_clicks):

            if n_clicks is not None and n_clicks == 0:
                return []

            return self._get_label_dropdowns()

        @self._app.callback(
            dash.dependencies.Output("de_tab", "style"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tabs_clicked_update_de_tab(_, differential_expression_clicks):

            if differential_expression_clicks is not None and \
                differential_expression_clicks > \
                self._differential_expression_n_clicks:

                self._differential_expression_n_clicks += 1
                return {"display": "block"}
            else:
                return {"display": "none"}

        @self._app.callback(
            dash.dependencies.Output("de_table", "children"),
            [dash.dependencies.Input("de_button", "n_clicks")],
            [dash.dependencies.State("subgroup_1_dropdown", "value"),
             dash.dependencies.State("subgroup_2_dropdown", "value")])
        def subgroup_changed(n_clicks, subgroup_1_labels, subgroup_2_labels):

            if n_clicks is None or n_clicks == 0:
                return []

            if len(subgroup_1_labels) > 0:

                if len(subgroup_2_labels) == 0:
                    subgroup_2_labels = None

                de = self._gene_expression_dataset.compare_gene_expression(
                    subgroup_1_labels, subgroup_2_labels)

                de = de.sort_values("p-value").iloc[:20]

                return self.generate_table(de)
            else:
                return []

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

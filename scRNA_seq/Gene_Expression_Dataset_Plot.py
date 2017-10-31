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
        self._tSNE_figure = None
        self._tabs = []
        self._data_containers = []
        self._cells = []
        self._de_stats = None
        self._de_stats_normalized = None
        self._differential_expression_n_clicks = 0
        self._clustering_n_clicks = 0
        self._de_start_range = 0
        self._subgroup_1_labels = []
        self._subgroup_2_labels = []
        self._n_clicks_de = 0
        self._n_clicks_de_previous = 0
        self._n_clicks_de_next = 0
        self._n_clicks_add_label = None
        self._n_clicks_delete_label = None

    def get_tSNE_figure(self, highlighted_clusters=None):

        gene_expression = \
            self._gene_expression_dataset.get_cell_gene_expression(
                Gene_Expression_Dataset.Transformation_Method.TSNE)

        highlighted_cells = \
            self._gene_expression_dataset.get_cells(highlighted_clusters)

        cell_read_counts = self._gene_expression_dataset._cell_transcript_counts

        x_values = []
        y_values = []
        colors = []
        hover_text = []
        self._cells = []

        for cell, gene_expression in gene_expression.iterrows():

            x_values.append(gene_expression['tSNE_1'])
            y_values.append(gene_expression['tSNE_2'])

            read_count = cell_read_counts.loc[cell].values[0]

            if cell in highlighted_cells:
                colors.append(read_count)
            else:
                colors.append('rgb(150, 150, 150')
            hover_text.append("%s<BR>%s" % (cell, read_count))
            self._cells.append(cell)

        figure = {
            'data': [
                go.Scatter(
                    x=x_values,
                    y=y_values,
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size': 5,
                        # 'line': {'width': 0.5, 'color': 'white'},
                        'color': colors,
                        'colorscale': 'Viridis',
                        'showscale': True
                    },
                    text=hover_text
                )
                ],
            'layout': go.Layout(
                xaxis={'title': 'tSNE 1'},
                yaxis={'title': 'tSNE 2'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': 0, 'y': 1},
                hovermode='closest',
                showlegend=False
            )
        }

        return figure

    @staticmethod
    def split_string(string):

        new_string = ""
        last_line_break_index = 0

        for character_index in range(len(string)):
            if character_index - last_line_break_index > 50:
                if string[character_index] == " ":
                    new_string += string[last_line_break_index:character_index]
                    new_string += "<BR>"
                    last_line_break_index = character_index + 1

        new_string += string[last_line_break_index:]

        return new_string

    @staticmethod
    def generate_differential_gene_expression_table(data_frame):

        if data_frame is None:
            return [
                dcc.Graph(
                    id="de_table_graph",
                    figure=ff.create_table(pandas.DataFrame()))]

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

        data_frame = data_frame.apply(lambda x: x.apply(lambda y: "%.3e" % y))

        return [
            dcc.Graph(
                id="de_table_graph", figure=ff.create_table(
                    data_frame, index=True, index_title="Gene",
                    hoverinfo="text", text=hover_texts),
                config={'displayModeBar': False})]

    @staticmethod
    def generate_label_counts_table(label_counts):

        label_counts = label_counts.sort_index()
        label_counts["# Cells"] =\
            label_counts["# Cells"].apply(lambda x: "%i" % x)
        label_counts["Ratio"] =\
            label_counts["Ratio"].apply(lambda x: "%.2f%%" % (x*100))

        return dcc.Graph(
            id="label_counts_figure", figure=ff.create_table(
                label_counts, index=True, index_title="Label"),
            config={'displayModeBar': False})

    def _get_gene_histograms(self, gene_index):

        if self._de_stats is None:
            return []

        gene = self._de_stats.index[gene_index + self._de_start_range]

        graph_title = "%s eCDF" % gene

        subgroup_1_title = " ".join(self._subgroup_1_labels)

        if self._subgroup_2_labels is None:
            subgroup_2_title = "Not %s" % " ".join(self._subgroup_1_labels)
        else:
            subgroup_2_title = "%s" % " ".join(self._subgroup_2_labels)

        group_1_counts = self._gene_expression_dataset.get_gene_counts(
            gene, self._subgroup_1_labels
        )
        group_1_counts = group_1_counts.value_counts()

        group_1_eCDF = group_1_counts.sort_index().cumsum() /\
            sum(group_1_counts.values)

        group_1_scatter = go.Scatter(
            x=group_1_eCDF.index,
            y=group_1_eCDF.values,
            name=subgroup_1_title
        )

        group_2_counts = self._gene_expression_dataset.get_gene_counts(
            gene, self._subgroup_2_labels
        )
        group_2_counts = group_2_counts.value_counts()

        group_2_eCDF = group_2_counts.sort_index().cumsum() /\
            sum(group_2_counts.values)

        group_2_scatter = go.Scatter(
            x=group_2_eCDF.index,
            y=group_2_eCDF.values,
            name=subgroup_2_title
        )

        layout = go.Layout(
            title=graph_title,
            xaxis=dict(
                title="Transcript Count"
            ),
            yaxis=dict(
                title="Cumulative Probability",
                range=[0, 1]
            ),
            hovermode="closest"
        )

        data = [group_1_scatter, group_2_scatter]

        graph = dcc.Graph(
            id="gene_eCDF",
            figure={
                "data": data,
                "layout": layout
            }
        )

        return [graph]

    def _get_de_plot(self):

        if self._de_stats_normalized is None:
            return {}

        x_values = []
        y_values = []
        gene_names = []

        for gene, gene_stats in self._de_stats_normalized.iterrows():
            x_values.append(gene_stats["Group 1 Mean"])
            y_values.append(gene_stats["Group 2 Mean"])
            gene_names.append(gene)

        min_value = min(min(x_values), min(y_values))
        max_value = max(max(x_values), max(y_values))

        gene_counts_scatter = go.Scatter(
            x=x_values,
            y=y_values,
            mode="markers",
            text=gene_names,
            name="Gene Expression")

        x_y_line = go.Scatter(
            x=[min_value, max_value],
            y=[min_value, max_value],
            name="X=Y"
        )

        subgroup_1_labels = " ".join(self._subgroup_1_labels)

        if self._subgroup_2_labels is None:
            subgroup_2_labels = "Not %s" % subgroup_1_labels
        else:
            subgroup_2_labels = " ".join(self._subgroup_2_labels)

        data = [gene_counts_scatter, x_y_line]
        layout = go.Layout(
            xaxis=dict(
                range=[min_value, max_value],
                title="%s Gene Expression" % subgroup_1_labels
            ),
            yaxis=dict(
                range=[min_value, max_value],
                title="%s Gene Expression" % subgroup_2_labels
            ),
            hovermode="closest",
            showlegend=False
        )

        return {
            "data": data,
            "layout": layout
        }

    def _get_label_options(self):

        labels = self._gene_expression_dataset.get_labels()

        label_options = [{"label": x, "value": x} for x in labels]

        return label_options

    def _get_label_dropdowns(self):

        label_options = self._get_label_options()

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

        differential_expression_tab = \
            html.Div(id="de_tab", children=[
                html.Div(
                    id="label_dropdowns",
                    children=self._get_label_dropdowns(),
                    style={'width': "50%"}),
                html.Div(id="navigation_pane", children=[
                    html.Button("Go!", id="de_button"),
                    html.Button("Previous", id="de_previous_button"),
                    html.Button("Next", id="de_next_button")
                    ]
                ),
                html.Div(id="de_pane", children=[
                    html.Div(
                        id="de_table", children=
                        Gene_Expression_Dataset_Plot.
                        generate_differential_gene_expression_table(None)
                    ),
                    html.Div(
                        id="de_plot_holder",
                        children=[
                            dcc.Graph(id="de_plot")
                        ]
                    ),
                    html.Div(
                        id="de_gene_range_label",
                        children=[],
                        style={"display": "none"}),
                    html.Div(id="de_data", children=[])
                ], style={"width": "60%", "display": "inline-block"}),
                html.Div(
                    id="gene_pane",
                    children=[],
                    style={
                        "width": "39%",
                        "display": "inline-block",
                        "vertical-align": "top"}),
            ])

        return differential_expression_tab

    def _get_tSNE_tab(self):

        tSNE_tab = html.Div(id="tSNE_tab", children=[

            html.Div(
                id="cluster_filter_container",
                children=[
                    dcc.Dropdown(
                        id="cluster_filter_dropdown",
                        options=self._get_label_options(),
                        value=[],
                        multi=True
                    )
                ],
                style={"width": "50%"}
            ),

            html.Div(
                id="tSNE_div",
                children=[

                    dcc.Graph(
                        id='tSNE',
                        figure=self._tSNE_figure
                    )
                ],
                style={
                    "width": "50%",
                    "display": "inline-block"}
            ),

            html.Div(
                id="label_table",
                children=[
                    Gene_Expression_Dataset_Plot.generate_label_counts_table(
                        self._gene_expression_dataset.get_label_counts()
                    )],
                style={
                    "width": "30%",
                    "display": "inline-block",
                    "vertical-align": "top"
                }
            ),

            html.Div(
                id="cluster_labeling",
                children=[

                    dcc.Input(id="label_name", type="text", value=""),
                    html.Button("Label Cluster", id="add_label_button"),
                ]
            ),
            html.Div(className="row", children=[

                html.Div([
                    html.Pre(id='selected-data')
                ])
            ]),
            html.Div(
                id="label_management",
                children=[
                    dcc.Dropdown(
                        id="manage_label_dropdown",
                        options=self._get_label_options(),
                        value=[]
                    ),
                    html.Button("Delete", id="delete_label_button"),
                    html.Label("", id="label_cell_count_text")
                ],
                style={"width": "25%"}
            )
            ]
        )

        return tSNE_tab

    def start(self):

        self._app = dash.Dash()

        self._tSNE_figure = self.get_tSNE_figure()

        self._tabs = [
            self._get_tSNE_tab(),
            self._get_differential_expression_tab()
        ]

        self._data_containers = [
            html.Div(id="labels", children=[], style={"display": "none"})
        ]

        self._app.layout = html.Div([

            html.Button("Clustering", id="clustering_button"),
            html.Button("Differential Expression",
                        id="differential_expression_button"),

            html.Div(id="tabs", children=self._tabs),
            html.Div(id="data", children=self._data_containers)
        ])

        @self._app.callback(
            dash.dependencies.Output("labels", "children"),
            [dash.dependencies.Input("delete_label_button", "n_clicks"),
             dash.dependencies.Input("add_label_button", "n_clicks")],
            [dash.dependencies.State("label_name", "value"),
             dash.dependencies.State("tSNE", "selectedData"),
             dash.dependencies.State("manage_label_dropdown", "value")]
        )
        def label_added_or_deleted(delete_label_n_clicks, add_label_n_clicks,
                                   label_name_to_add, selected_data,
                                   label_to_delete):

            current_labels = self._gene_expression_dataset.get_labels()

            if (delete_label_n_clicks is None or delete_label_n_clicks == 0)\
                    and (add_label_n_clicks is None or add_label_n_clicks == 0):
                return current_labels

            if delete_label_n_clicks == self._n_clicks_delete_label:
                self._n_clicks_add_label = add_label_n_clicks

                if selected_data is None:
                    self._selected_points = None
                    return current_labels

                if label_name_to_add == "":
                    return current_labels

                self._selected_points = selected_data["points"]

                point_indices = [
                    point["pointNumber"] for point in
                    self._selected_points]

                cells = [self._cells[i] for i in point_indices]

                self._gene_expression_dataset.label_cells(label_name_to_add,
                                                          cells)

                return self._get_label_dropdowns()
            else:
                self._n_clicks_delete_label = delete_label_n_clicks
                self._gene_expression_dataset.delete_label(label_to_delete)

                return self._get_label_dropdowns()

        @self._app.callback(
            dash.dependencies.Output("manage_label_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")]
        )
        def update_label_management_dropdown(label_dropdowns):
            return self._get_label_options()

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
            [dash.dependencies.Input("labels", "children")])
        def label_added_update_label_dropdowns(n_clicks):

            return self._get_label_dropdowns()

        @self._app.callback(
            dash.dependencies.Output("label_name", "value"),
            [dash.dependencies.Input("add_label_button", "n_clicks")])
        def label_added_clear_label_field(_):
            return ""

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
            dash.dependencies.Output("de_gene_range_label", "children"),
            [dash.dependencies.Input("de_previous_button", "n_clicks"),
             dash.dependencies.Input("de_next_button", "n_clicks")])
        def update_ge_range(n_clicks_de_previous, n_clicks_de_next):

            if n_clicks_de_next is None or n_clicks_de_previous is None:
                return

            if n_clicks_de_next > self._n_clicks_de_next:
                self._n_clicks_de_next = n_clicks_de_next
                self._de_start_range += 20

            if n_clicks_de_previous > self._n_clicks_de_previous:
                self._n_clicks_de_previous = n_clicks_de_previous

                if self._de_start_range != 0:
                    self._de_start_range -= 20
            #
            # label_text = "Genes %i-%i" % (
            #     self._de_start_range, self._de_start_range + 20)

            return

        @self._app.callback(
            dash.dependencies.Output("de_next_button", "n_clicks"),
            [dash.dependencies.Input("de_button", "n_clicks")])
        def clear_next_button_n_clicks(n_clicks):
            if n_clicks is None or n_clicks == 0:
                return
            else:
                self._n_clicks_de_next = 0
                self._n_clicks_de_previous = 0
                return 0

        @self._app.callback(
            dash.dependencies.Output("de_previous_button", "n_clicks"),
            [dash.dependencies.Input("de_button", "n_clicks")])
        def clear_previous_button_n_clicks(n_clicks):
            if n_clicks is None or n_clicks == 0:
                return
            else:
                return 0

        @self._app.callback(
            dash.dependencies.Output("tSNE", "figure"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value")])
        def update_plot_from_cluster_filter(selected_clusters):

            if not selected_clusters:
                selected_clusters = None

            return self.get_tSNE_figure(highlighted_clusters=selected_clusters)

        @self._app.callback(
            dash.dependencies.Output("label_table", "children"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value"),
             dash.dependencies.Input("labels", "children")])
        def cluster_filter_updated(cluster_filter_values, clusters):

            label_counts = self._gene_expression_dataset.get_label_counts(
                cluster_filter_values)

            return self.generate_label_counts_table(label_counts)

        @self._app.callback(
            dash.dependencies.Output("de_plot", "figure"),
            [dash.dependencies.Input("de_data", "children")])
        def de_clicked_update_plot(_):
            return self._get_de_plot()

        @self._app.callback(
            dash.dependencies.Output("de_data", "children"),
            [dash.dependencies.Input("de_button", "n_clicks")],
            [dash.dependencies.State("subgroup_1_dropdown", "value"),
             dash.dependencies.State("subgroup_2_dropdown", "value")])
        def new_de_clicked(n_clicks_de, subgroup_1_labels, subgroup_2_labels):

            if n_clicks_de is None or n_clicks_de == 0:
                return []

            if n_clicks_de > self._n_clicks_de:
                self._de_start_range = 0
                self._n_clicks_de = n_clicks_de

            if len(subgroup_1_labels) > 0:

                if self._subgroup_1_labels != subgroup_1_labels or \
                        self._subgroup_2_labels != subgroup_2_labels:

                    if len(subgroup_2_labels) == 0:
                        subgroup_2_labels = None

                    self._de_stats = \
                        self._gene_expression_dataset.compare_gene_expression(
                            subgroup_1_labels, subgroup_2_labels,
                            use_normalized=False)
                    self._de_stats_normalized = \
                        self._gene_expression_dataset.compare_gene_expression(
                            subgroup_1_labels, subgroup_2_labels,
                            use_normalized=True)
                    self._de_stats["Log2 Change Abs"] = \
                        abs(self._de_stats["Log2 Change"])
                    self._de_stats = self._de_stats.sort_values(
                        ["p-value", "Log2 Change Abs"],
                        ascending=[True, False]
                    )
                    self._de_stats = self._de_stats.drop("Log2 Change Abs",
                                                         axis=1)

                    self._subgroup_1_labels = subgroup_1_labels
                    self._subgroup_2_labels = subgroup_2_labels
            return []

        @self._app.callback(
            dash.dependencies.Output("de_table", "children"),
            [dash.dependencies.Input("de_data", "children"),
             dash.dependencies.Input("de_gene_range_label", "children")])
        def get_differential_expression_clicked(_1, _2):

            if self._de_stats is not None:

                de = self._de_stats.iloc[
                     self._de_start_range:self._de_start_range+20]

                return self.generate_differential_gene_expression_table(de)
            else:
                return []

        @self._app.callback(
            dash.dependencies.Output("gene_pane", "children"),
            [dash.dependencies.Input("de_table_graph", "clickData")]
        )
        def de_figure_clicked(click_data):

            if click_data is None:
                return []

            if click_data["points"][0]["y"] == 0:
                return []

            return self._get_gene_histograms(click_data["points"][0]["y"] - 1)

        @self._app.callback(
            dash.dependencies.Output("cluster_filter_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")])
        def data_edited(labels):
            return self._get_label_options()

        self._app.run_server()

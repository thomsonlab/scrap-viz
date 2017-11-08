import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from .Gene_Expression_Dataset import Gene_Expression_Dataset
import plotly.figure_factory as ff
from scRNA_seq import Gene
import pandas
import json
import math

verbose = True


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
        self._n_clicks_filter_by = 0
        self._n_clicks_add_label = None
        self._n_clicks_delete_label = None
        self._column_to_sort = 1
        self._regenerate_de = False
        self._de_figure = None

    def get_tSNE_figure(self, highlighted_cells=None, cell_color_values=None):

        gene_expression = \
            self._gene_expression_dataset.get_cell_gene_expression(
                Gene_Expression_Dataset.Transformation_Method.TSNE)

        x_values = []
        y_values = []
        colors = []
        hover_text = []
        self._cells = []

        for cell, gene_expression in gene_expression.iterrows():

            x_values.append(gene_expression['tSNE_1'])
            y_values.append(gene_expression['tSNE_2'])

            if cell_color_values is not None:
                color_value = cell_color_values[cell]
            else:
                color_value = "rgb(0, 0, 255)"

            if highlighted_cells is not None and cell in highlighted_cells:
                colors.append(color_value)
            else:
                colors.append("rgba(150, 150, 150, 0.25")
            hover_text.append("%s<BR>%s" % (cell, color_value))
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
                    text=hover_text,
                    hoverinfo="text"
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
            return ff.create_table(pandas.DataFrame())

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

        return ff.create_table(
            data_frame, index=True, index_title="Gene", hoverinfo="text",
            text=hover_texts)

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

    def _get_gene_eCDF(self, gene_index):

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
        group_1_counts = group_1_counts.value_counts().sort_index()
        group_1_cum_sum = group_1_counts.cumsum()

        group_1_eCDF = group_1_cum_sum / sum(group_1_counts.values)

        group_1_hover_text = []
        for gene_count_index in range(len(group_1_eCDF)):
            gene_count = group_1_eCDF.index[gene_count_index]
            eCDF_value = group_1_eCDF.values[gene_count_index]
            num_cells = group_1_cum_sum.values[gene_count_index]
            group_1_hover_text.append(
                "%i cell(s) (%.2f%%)<BR>with count <= %f" % (
                    num_cells, eCDF_value*100, gene_count))

        group_1_scatter = go.Scatter(
            x=group_1_eCDF.index,
            y=group_1_eCDF.values,
            name=subgroup_1_title,
            text=group_1_hover_text,
            hoverinfo="text"
        )

        group_2_counts = self._gene_expression_dataset.get_gene_counts(
            gene, self._subgroup_2_labels
        )
        group_2_counts = group_2_counts.value_counts().sort_index()
        group_2_cum_sum = group_2_counts.cumsum()

        group_2_eCDF = group_2_cum_sum / sum(group_2_counts.values)

        group_2_hover_text = []
        for gene_count_index in range(len(group_2_eCDF)):
            gene_count = group_2_eCDF.index[gene_count_index]
            eCDF_value = group_2_eCDF.values[gene_count_index]
            num_cells = group_2_cum_sum.values[gene_count_index]
            group_2_hover_text.append(
                "%i cell(s) (%.2f%%)<BR>with count <= %f" % (
                    num_cells, eCDF_value*100, gene_count))

        group_2_scatter = go.Scatter(
            x=group_2_eCDF.index,
            y=group_2_eCDF.values,
            name=subgroup_2_title,
            text=group_2_hover_text,
            hoverinfo="text"
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

        labels = sorted(self._gene_expression_dataset.get_labels())

        label_options = [{"label": x, "value": x} for x in labels]

        return label_options

    def _get_gene_options(self):

        genes = sorted(self._gene_expression_dataset.get_genes())

        gene_options = [{"label": x, "value": x} for x in genes]

        return gene_options

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

        differential_expression_tab = html.Div(
            id="de_tab",
            children=[
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
                        id="de_table", children=[
                            dcc.Graph(
                                id="de_table_graph",
                                figure=go.Figure(),
                                config={'displayModeBar': False})
                        ]
                    ),
                    html.Div(
                        id="de_plot_holder",
                        children=[
                            dcc.Graph(id="de_plot")
                        ]
                    )
                ], style={"width": "60%", "display": "inline-block"}),
                html.Div(
                    id="gene_pane",
                    children=[],
                    style={
                        "width": "39%",
                        "display": "inline-block",
                        "vertical-align": "top"})
            ]
        )

        return differential_expression_tab

    @staticmethod
    def get_cell_value_range_slider(min_value, max_value):

        value_range = max_value - min_value
        log_range = math.ceil(math.log10(value_range))
        value_range = math.pow(10, log_range)
        interval = value_range/10

        min_value = math.floor(min_value/interval)*interval
        max_value = math.ceil(max_value/interval)*interval

        marks = {}
        mark_value = min_value

        while mark_value <= max_value:

            if round(mark_value) == mark_value:
                mark_value = int(mark_value)

            marks[mark_value] = {'label': mark_value}
            mark_value += interval

        return dcc.RangeSlider(
            id="cell_value_range_slider",
            min=min_value,
            max=max_value,
            step=interval/5,
            value=[min_value, max_value],
            allowCross=False,
            marks=marks
        )

    def _get_tSNE_tab(self):

        tSNE_tab = html.Div(
            id="tSNE_tab",
            children=[
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
                        "display": "inline-block",
                        "marginTop": 25}
                ),
                html.Div(
                    id="gene_filter_management",
                    children=[
                        html.Div(
                            id="gene_filter_dropdown_div",
                            children=[
                                dcc.Dropdown(
                                    id="gene_filter_dropdown",
                                    options=self._get_gene_options(),
                                    value=[],
                                    multi=False
                                )
                            ],
                            style={
                                "width": "50%"
                            }
                        ),
                        html.Div(
                            id="cell_value_range_slider_div",
                            children=[
                                self.get_cell_value_range_slider(-2, 2)
                            ],
                            style={
                                "marginTop": 10,
                                "marginBottom": 20,
                                "width": "75%"
                            }
                        )
                    ],
                    style={
                        "width": "45%",
                        "display": "inline-block",
                        "vertical-align": "top",
                        "marginLeft": 25
                    }
                ),
                html.Div(
                    id="selection_labeling",
                    children=[
                        dcc.Input(
                            id="label_name",
                            type="text", value=""
                        ),
                        html.Button(
                            "Label Selected",
                            id="add_label_button"
                        )
                    ],
                    style={
                        "width": "25%"
                    }
                ),
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
                    style={
                        "width": "25%",
                        "marginTop": 10,
                        "marginBottom": 10
                    }
                ),
                html.Div(
                    id="label_table",
                    children=[
                        Gene_Expression_Dataset_Plot.
                        generate_label_counts_table(
                            self._gene_expression_dataset.get_label_counts()
                        )
                    ],
                    style={
                        "width": "50%",
                        "display": "inline-block",
                        "vertical-align": "top"
                    }
                ),
                html.Div(
                    id="gene_count_table",
                    children=[],
                    style={
                        "width": "45%",
                        "display": "inline-block",
                        "marginLeft": 25
                    }
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
            html.Div(id="labels", children=[], style={"display": "none"}),
            html.Div(id="de_data", children=[], style={"display": "none"}),
            html.Div(id="gene_range", children=[], style={"display": "none"}),
            html.Div(id="cell_color_values", children=[],
                     style={"display": "none"}),
            html.Div(id="unfiltered_cells", children=[],
                     style={"display": "none"})
        ]

        self._app.layout = html.Div([

            html.Button("Clustering", id="clustering_button"),
            html.Button("Differential Expression",
                        id="differential_expression_button"),

            html.Div(id="tabs", children=self._tabs, style={"marginTop": 10}),
            html.Div(id="data", children=self._data_containers)
        ])

        @self._app.callback(
            dash.dependencies.Output("manage_label_dropdown", "value"),
            [dash.dependencies.Input("delete_label_button", "n_clicks")]
        )
        def delete_button_clicked_clear_dropdown(n_clicks):

            if verbose:
                print("delete_button_clicked_clear_dropdown")
            return ""

        @self._app.callback(
            dash.dependencies.Output("labels", "children"),
            [dash.dependencies.Input("delete_label_button", "n_clicks"),
             dash.dependencies.Input("add_label_button", "n_clicks")],
            [dash.dependencies.State("label_name", "value"),
             dash.dependencies.State("tSNE", "selectedData"),
             dash.dependencies.State("manage_label_dropdown", "value"),
             dash.dependencies.State("unfiltered_cells", "children")]
        )
        def label_added_or_deleted(delete_label_n_clicks, add_label_n_clicks,
                                   label_name_to_add, selected_data,
                                   label_to_delete, unfiltered_cells):

            if verbose:
                print("label_added_or_deleted")

            current_labels = self._gene_expression_dataset.get_labels()

            if (delete_label_n_clicks is None or delete_label_n_clicks == 0)\
                    and (add_label_n_clicks is None or add_label_n_clicks == 0):
                return current_labels

            # If n_clicks of delete button is the same, this is an add label
            if delete_label_n_clicks == self._n_clicks_delete_label:
                self._n_clicks_add_label = add_label_n_clicks

                if label_name_to_add == "":
                    return current_labels

                if unfiltered_cells:
                    unfiltered_cells = json.loads(unfiltered_cells)
                    if not unfiltered_cells:
                        unfiltered_cells = None
                else:
                    unfiltered_cells = None

                if selected_data:
                    selected_points = selected_data["points"]

                    point_indices = [
                        point["pointNumber"] for point in selected_points]

                    selected_cells = [self._cells[i] for i in point_indices]
                else:
                    selected_cells = self._cells

                selected_cells = set(selected_cells)
                unfiltered_cells = set(unfiltered_cells)

                cells_to_label = selected_cells.intersection(unfiltered_cells)

                self._gene_expression_dataset.label_cells(label_name_to_add,
                                                          cells_to_label)
                self._gene_expression_dataset.save()

                return self._get_label_dropdowns()
            else:
                self._n_clicks_delete_label = delete_label_n_clicks
                self._gene_expression_dataset.delete_label(label_to_delete)
                self._gene_expression_dataset.save()

                return self._get_label_dropdowns()

        @self._app.callback(
            dash.dependencies.Output("manage_label_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")]
        )
        def update_label_management_dropdown(label_dropdowns):
            if verbose:
                print("update_label_management_dropdown")
            return self._get_label_options()

        @self._app.callback(
            dash.dependencies.Output("tSNE_tab", "style"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tabs_clicked_update_tSNE_tab(clustering_button_clicks, _):
            if verbose:
                print("tabs_clicked_update_tSNE_tab")

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
            if verbose:
                print("label_added_update_label_dropdowns")

            return self._get_label_dropdowns()

        @self._app.callback(
            dash.dependencies.Output("label_name", "value"),
            [dash.dependencies.Input("add_label_button", "n_clicks")])
        def label_added_clear_label_field(_):
            if verbose:
                print("label_added_clear_label_field")
            return ""

        @self._app.callback(
            dash.dependencies.Output("de_tab", "style"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tabs_clicked_update_de_tab(_, differential_expression_clicks):
            if verbose:
                print("tabs_clicked_update_de_tab")

            if differential_expression_clicks is not None and \
                    differential_expression_clicks > \
                    self._differential_expression_n_clicks:

                self._differential_expression_n_clicks += 1
                return {"display": "block"}
            else:
                return {"display": "none"}

        @self._app.callback(
            dash.dependencies.Output("gene_range", "children"),
            [dash.dependencies.Input("de_previous_button", "n_clicks"),
             dash.dependencies.Input("de_next_button", "n_clicks")])
        def update_ge_range(n_clicks_de_previous, n_clicks_de_next):

            if verbose:
                print("update_ge_range")

            if n_clicks_de_next is None or n_clicks_de_previous is None:
                pass
            elif n_clicks_de_next > self._n_clicks_de_next:
                self._n_clicks_de_next = n_clicks_de_next
                self._de_start_range += 20
                self._regenerate_de = True
            elif n_clicks_de_previous > self._n_clicks_de_previous:
                self._n_clicks_de_previous = n_clicks_de_previous

                if self._de_start_range != 0:
                    self._de_start_range -= 20
                    self._regenerate_de = True

            gene_range = {"start": self._de_start_range}

            return json.dumps(gene_range)

        @self._app.callback(
            dash.dependencies.Output("de_next_button", "n_clicks"),
            [dash.dependencies.Input("de_button", "n_clicks")])
        def clear_next_button_n_clicks(n_clicks):
            if verbose:
                print("clear_next_button_n_clicks")
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
            if verbose:
                print("clear_previous_button_n_clicks")
            if n_clicks is None or n_clicks == 0:
                return
            else:
                return 0

        @self._app.callback(
            dash.dependencies.Output("cell_value_range_slider_div", "children"),
            [dash.dependencies.Input("cell_color_values", "children")])
        def update_cell_value_range_slider(cell_color_values):

            if verbose:
                print("update_cell_value_range_slider")
                print(cell_color_values)

            if not cell_color_values:
                return []

            de_data = json.loads(cell_color_values)

            if not de_data:
                return []

            min_value = min(de_data.values())
            max_value = max(de_data.values())

            return [self.get_cell_value_range_slider(min_value, max_value)]

        @self._app.callback(
            dash.dependencies.Output("cell_color_values", "children"),
            [dash.dependencies.Input("gene_filter_dropdown", "value")])
        def update_cell_color_values(gene):

            if not gene:
                cell_read_counts = \
                    self._gene_expression_dataset.get_cell_transcript_counts()
                cell_read_counts = cell_read_counts["TOTAL_TRANSCRIPT_COUNT"]
                cell_read_count_ints = \
                    {gene: int(count)
                        for gene, count in cell_read_counts.items()}
                return json.dumps(cell_read_count_ints)

            de = self._gene_expression_dataset.get_cell_gene_differential(gene)

            return json.dumps(de.to_dict())

        @self._app.callback(
            dash.dependencies.Output("unfiltered_cells", "children"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value"),
             dash.dependencies.Input("cell_value_range_slider", "value")],
            [dash.dependencies.State("cell_color_values", "children")])
        def set_unfiltered_cells(selected_clusters, cell_value_range,
                                 cell_color_values):

            if verbose:
                print("set_unfiltered_cells")

            if not selected_clusters:
                selected_clusters = None

            if cell_color_values:
                cells_gene_de = json.loads(cell_color_values)
                if not cells_gene_de:
                    cells_gene_de = None
            else:
                cells_gene_de = None

            if not cell_value_range:
                cell_value_range = None

            cells_in_cluster = \
                self._gene_expression_dataset.get_cells(selected_clusters)

            cells_to_remove = set()
            if cell_value_range is not None and cells_gene_de is not None:
                for cell in cells_in_cluster:
                    cell_color_value = cells_gene_de[cell]
                    if cell_color_value < cell_value_range[0] \
                            or cell_color_value > cell_value_range[1]:
                        cells_to_remove.add(cell)
                unfiltered_cells = cells_in_cluster.difference(cells_to_remove)
            else:
                unfiltered_cells = cells_in_cluster

            return json.dumps(list(unfiltered_cells))

        @self._app.callback(
            dash.dependencies.Output("tSNE", "figure"),
            [dash.dependencies.Input("unfiltered_cells", "children")],
            [dash.dependencies.State("cell_color_values", "children")])
        def update_plot_from_filters(
                unfiltered_cells, cell_color_values):

            if verbose:
                print("update_plot_from_cluster_filter")

            if cell_color_values:
                cell_color_values = json.loads(cell_color_values)
                if not cell_color_values:
                    cell_color_values = None
            else:
                cell_color_values = None

            if unfiltered_cells:
                unfiltered_cells = json.loads(unfiltered_cells)
            else:
                unfiltered_cells = None

            return self.get_tSNE_figure(
                highlighted_cells=unfiltered_cells,
                cell_color_values=cell_color_values)

        @self._app.callback(
            dash.dependencies.Output("label_table", "children"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value"),
             dash.dependencies.Input("labels", "children")])
        def cluster_filter_updated(cluster_filter_values, clusters):

            if verbose:
                print("cluster_filter_updated")

            label_counts = self._gene_expression_dataset.get_label_counts(
                cluster_filter_values)

            return self.generate_label_counts_table(label_counts)

        @self._app.callback(
            dash.dependencies.Output("de_plot", "figure"),
            [dash.dependencies.Input("de_data", "children")])
        def de_clicked_update_plot(de_data):
            if verbose:
                print("de_clicked_update_plot")
            return self._get_de_plot()

        @self._app.callback(
            dash.dependencies.Output("de_data", "children"),
            [dash.dependencies.Input("de_button", "n_clicks")],
            [dash.dependencies.State("subgroup_1_dropdown", "value"),
             dash.dependencies.State("subgroup_2_dropdown", "value")])
        def new_de_clicked(n_clicks_de, subgroup_1_labels, subgroup_2_labels):
            if verbose:
                print("new_de_clicked")

            de_data = {"ready": False}

            if n_clicks_de is None or n_clicks_de == 0:
                print("Setting de_data to:")
                print(de_data)
                return json.dumps(de_data)

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

                    self._regenerate_de = True
                    de_data["ready"] = True

            print("Setting de_data to:")
            print(de_data)
            return json.dumps(de_data)

        @self._app.callback(
            dash.dependencies.Output("de_table_graph", "figure"),
            [dash.dependencies.Input("de_data", "children"),
             dash.dependencies.Input("gene_range", "children"),
             dash.dependencies.Input("de_table_graph", "clickData")])
        def get_differential_expression_clicked(
                de_data, gene_range, click_data):
            if verbose:
                print("get_differential_expression_clicked")

            if click_data is not None:
                x = click_data["points"][0]["x"]
                y = click_data["points"][0]["y"]
                if y == 0 and x > 0:
                    if x - 1 != self._column_to_sort:
                        self._column_to_sort = x - 1
                        self._de_start_range = 0

                        if self._column_to_sort == 0:
                            self._de_stats["Log2 Change Abs"] = \
                                abs(self._de_stats["Log2 Change"])
                            self._de_stats = self._de_stats.sort_values(
                                "Log2 Change Abs", ascending=False)
                            self._de_stats = self._de_stats.drop(
                                "Log2 Change Abs",
                                axis=1)
                        elif self._column_to_sort == 1:
                            self._de_stats["Log2 Change Abs"] = \
                                abs(self._de_stats["Log2 Change"])
                            self._de_stats = self._de_stats.sort_values(
                                ["p-value", "Log2 Change Abs"],
                                ascending=[True, False]
                            )
                            self._de_stats = self._de_stats.drop(
                                "Log2 Change Abs",
                                axis=1)
                        else:
                            column_to_sort_by = \
                                self._de_stats.columns[self._column_to_sort]
                            self._de_stats = self._de_stats.sort_values(
                                column_to_sort_by, ascending=False)
                        self._regenerate_de = True

            if self._de_stats is not None:
                if self._regenerate_de:

                    self._regenerate_de = False

                    de = self._de_stats.iloc[
                         self._de_start_range:self._de_start_range+20]

                    self._de_figure = \
                        self.generate_differential_gene_expression_table(de)

                    return self._de_figure
                else:
                    return self._de_figure
            else:
                self._de_figure = go.Figure()
                return self._de_figure

        @self._app.callback(
            dash.dependencies.Output("gene_pane", "children"),
            [dash.dependencies.Input("de_table_graph", "clickData")]
        )
        def de_figure_clicked(click_data):
            if verbose:
                print("de_figure_clicked")

            if click_data is None:
                return []

            if click_data["points"][0]["y"] == 0:
                return []

            return self._get_gene_eCDF(click_data["points"][0]["y"] - 1)

        @self._app.callback(
            dash.dependencies.Output("cluster_filter_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")])
        def data_edited(labels):
            if verbose:
                print("data_edited")
            return self._get_label_options()

        @self._app.callback(
            dash.dependencies.Output("gene_count_table", "children"),
            [dash.dependencies.Input("tSNE", "clickData")])
        def tSNE_clicked(click_data):

            if click_data is None:
                return []

            click_text = click_data["points"][0]["text"]
            cell_name = click_text[:click_text.find("<BR>")]

            gene_counts = pandas.DataFrame(
                self._gene_expression_dataset.
                get_gene_counts_for_cell(cell_name))

            gene_counts.columns = ["Count"]

            gene_counts = gene_counts[gene_counts != 0].\
                sort_values("Count", ascending=False)[0:50]

            gene_counts["Count"] = gene_counts["Count"].apply(
                lambda x: "%.3e" % x)

            gene_counts.columns = ["%s Count" % cell_name]

            return dcc.Graph(
                id="gene_counts_figure", figure=ff.create_table(
                    gene_counts, index=True, index_title="Gene"),
                config={'displayModeBar': False})

        self._app.run_server()

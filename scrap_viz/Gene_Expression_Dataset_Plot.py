import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.figure_factory as ff
import pandas
import json
import math
import os
from flask import send_from_directory
from matplotlib import pyplot

from scrap.dataset import Transformation_Method
from scrap.dataset import Clustering_Method
from scrap.utils import fileio

verbose = True


class Gene_Expression_Dataset_Plot:

    def __init__(self, gene_expression_dataset, gene_metadata, port):

        self._app = None
        self._gene_expression_dataset = gene_expression_dataset
        self._gene_metadata = gene_metadata
        self._projection_figure = None
        self._tabs = []
        self._data_containers = []
        self._cells = []
        self._de_stats = None
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
        self._n_clicks_auto_cluster = None
        self._n_clicks_slider = None
        self._column_to_sort = 1
        self._regenerate_de = False
        self._de_figure = None
        self._port = port
        self._eCDF_gene_index = None

    def get_projection_figure(
            self, transformation_method=None,
            highlighted_cells=None, cell_color_values=None,
            color_by_gene_count=True,
            x_axis_dimension=0, y_axis_dimension=1):

        if transformation_method is None:
            transformation_method = self._get_default_transformation_method()

        print("Getting gene expression...")

        gene_expression = \
            self._gene_expression_dataset.get_cell_gene_expression(
                transformation_method)

        if highlighted_cells is None or len(highlighted_cells) == 1 or \
                color_by_gene_count is True:

            x_values = []
            y_values = []
            colors = []
            hover_text = []
            self._cells = []

            if highlighted_cells is not None:

                print("Getting highlighted cells...")

                all_highlighted_cells = highlighted_cells[
                    sorted(highlighted_cells.keys())[0]]
                for i in range(1, len(highlighted_cells)):
                    all_highlighted_cells.extend(
                        highlighted_cells[sorted(highlighted_cells.keys())[i]])

            print("Getting x/y values, color, and hover text for all cells...")
            for cell, cell_gene_expression in gene_expression.iterrows():

                x_values.append(
                    cell_gene_expression[
                        gene_expression.columns[x_axis_dimension]])
                y_values.append(
                    cell_gene_expression[
                        gene_expression.columns[y_axis_dimension]])

                if cell_color_values is not None:
                    color_value = cell_color_values[cell]
                else:
                    color_value = "rgb(0, 0, 255)"

                if highlighted_cells is not None and \
                        cell in all_highlighted_cells:
                    colors.append(color_value)
                else:
                    colors.append("rgba(150, 150, 150, 0.25)")
                hover_text.append(
                    "%s<BR>%s" % (cell, color_value))

                self._cells.append(cell)

            print("Creating figure object...")

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
                    )],
                'layout': go.Layout(
                    xaxis={
                        'title': gene_expression.columns[x_axis_dimension]},
                    yaxis={
                        'title': gene_expression.columns[y_axis_dimension]},
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                    legend={'x': 0, 'y': 1},
                    hovermode='closest',
                    showlegend=False
                )
            }

        else:
            x_values = {}
            y_values = {}
            label_colors = {}
            hover_text = {}
            self._cells = []

            colors = self.get_colors(len(highlighted_cells))

            if highlighted_cells is not None:

                print("Getting highlighted cells...")
                for label_index, label in enumerate(sorted(highlighted_cells)):
                    x_values[label] = []
                    y_values[label] = []
                    label_colors[label] = colors[label_index]
                    hover_text[label] = []
            else:
                highlighted_cells = {}

            x_values['Other'] = []
            y_values['Other'] = []
            label_colors['Other'] = "rgba(150, 150, 150, 0.25)"
            hover_text['Other'] = []

            print("Getting x/y values, color, and hover text for all cells...")
            for cell, cell_gene_expression in gene_expression.iterrows():

                cell_label = 'Other'

                for label in highlighted_cells:
                    if cell in highlighted_cells[label]:
                        cell_label = label
                        x_values[label].append(
                            cell_gene_expression[gene_expression.columns[x_axis_dimension]])
                        y_values[label].append(
                            cell_gene_expression[gene_expression.columns[y_axis_dimension]])

                        if cell_color_values is not None:
                            color_value = cell_color_values[cell]
                        else:
                            color_value = "rgb(0, 0, 255)"

                        hover_text[label].append("%s<BR>%s" % (cell, color_value))

                if cell_label == 'Other':
                    x_values[cell_label].append(
                        cell_gene_expression[
                            gene_expression.columns[x_axis_dimension]])
                    y_values[cell_label].append(
                        cell_gene_expression[
                            gene_expression.columns[y_axis_dimension]])

                    if cell_color_values is not None:
                        color_value = cell_color_values[cell]
                    else:
                        color_value = "rgb(0, 0, 255)"

                    hover_text[cell_label].append("%s<BR>%s" % (cell, color_value))

                self._cells.append(cell)

            print("Creating figure object...")

            figure = {
                'data': [
                    go.Scatter(
                        x=x_values[label],
                        y=y_values[label],
                        mode='markers',
                        opacity=0.7,
                        marker={
                            'size': 5,
                            'color': label_colors[label]
                        },
                        text=hover_text[label],
                        hoverinfo="text",
                        name=label
                    ) for label in x_values],
                'layout': go.Layout(
                    xaxis={'title': gene_expression.columns[x_axis_dimension]},
                    yaxis={'title': gene_expression.columns[y_axis_dimension]},
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                    legend={'x': 0, 'y': 1},
                    hovermode='closest',
                    showlegend=True
                )
            }

        if verbose:
            print("Returning projection figure")

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
    def get_colors(num_colors):

        color_map = pyplot.get_cmap("viridis")

        num_possible_colors = len(color_map.colors)

        spacing = math.floor(num_possible_colors/num_colors)

        colors = []

        for i in range(num_colors):

            color = list(color_map.colors[spacing*i])
            for channel_index, channel_value in enumerate(color):
                color[channel_index] = int(channel_value * 255)

            colors.append("rgb(%i, %i, %i)" % (color[0], color[1], color[2]))

        return colors

    def generate_differential_gene_expression_table(self, data_frame):

        if data_frame is None:
            return ff.create_table(pandas.DataFrame())

        genes = data_frame.index

        gene_name_descriptions = \
            self._gene_metadata.get_gene_summaries(genes)

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

        for column in data_frame.columns:
            if column == "Cluster":
                continue

            data_frame[column] = data_frame[column].apply(
                lambda x: "%.3e" % x)

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

        self._eCDF_gene_index = gene_index

        gene = self._de_stats.index[gene_index]

        if "Cluster" in self._de_stats.columns:
            additional_cluster = self._de_stats.iloc[
                gene_index]["Cluster"]
            graph_title = "%s " % additional_cluster
            group_1_labels = list(self._subgroup_1_labels)
            group_1_labels.append(additional_cluster)
            if self._subgroup_2_labels is None:
                group_2_labels = list(self._subgroup_2_labels)
                group_2_labels.append(additional_cluster)
            else:
                group_2_labels = [additional_cluster]
        else:
            group_1_labels = self._subgroup_1_labels
            group_2_labels = self._subgroup_2_labels
            graph_title = ""

        graph_title += "%s eCDF" % gene

        subgroup_1_title = " ".join(self._subgroup_1_labels)

        if self._subgroup_2_labels is None:
            subgroup_2_title = "Not %s" % " ".join(self._subgroup_1_labels)
        else:
            subgroup_2_title = "%s" % " ".join(self._subgroup_2_labels)

        group_1_counts = self._gene_expression_dataset.get_gene_counts(
            gene, group_1_labels
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
            gene, group_2_labels
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
                range=[0, 1.05]
            ),
            hovermode="closest"
        )

        data = [group_1_scatter, group_2_scatter]

        graph = dcc.Graph(
            id="gene_eCDF",
            figure={
                "data": data,
                "layout": layout,
                "config": {
                    "displayModeBar": False,
                    "displaylogo": False
                }
            }
        )

        return [graph]

    def _get_de_plot(self):

        if self._de_stats is None:
            return {}

        x_values = []
        y_values = []
        gene_names = []

        for gene, gene_stats in self._de_stats.iterrows():
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

    @staticmethod
    def _get_label_type_options():

        label_types = ["L2FC", "Count", "Normalized"]

        label_type_options = [{"label": x, "value": x} for x in label_types]

        return label_type_options

    @staticmethod
    def _get_transformation_method_options():

        label_options = []

        for method_name, method in \
                Transformation_Method.__members__.items():
            label_options.append({"label": method_name, "value": method_name})

        return label_options

    @staticmethod
    def _get_default_transformation_method():

        return Transformation_Method.TSNE

    @staticmethod
    def _get_cluster_method_options():

        label_options = []

        for method_name, method in \
                Clustering_Method.__members__.items():
            label_options.append({"label": method_name, "value": method_name})

        return label_options

    @staticmethod
    def _get_default_cluster_method():

        return Clustering_Method.K_MEANS

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
            ),
            html.Label("Differential Across"),
            dcc.Dropdown(
                id="de_across_dropdown",
                options=label_options,
                value=[],
                multi=True
            )
        ]

        return label_dropdowns

    def _get_differential_expression_tab(self):

        print("Getting DE tab...")

        differential_expression_tab = html.Div(
            id="de_tab", className="tabcontent",
            children=[
                html.Div(
                    id="label_dropdowns",
                    children=self._get_label_dropdowns(),
                    style={'width': "50%"}),
                dcc.Checklist(
                    id="use_normalized_checklist",
                    options=[
                        {
                            'label': 'Use Normalized',
                            'value': 'use_normalized'
                        }
                    ],
                    values=['use_normalized']
                ),
                html.Div(
                    id="navigation_pane",
                    children=[
                        html.Button("Go!", id="de_button"),
                        html.Button("Previous", id="de_previous_button"),
                        html.Button("Next", id="de_next_button"),
                        html.Button("Export", id="export_de_button")
                    ],
                    style={
                        "marginTop": "1%",
                        "marginBottom": "1%"
                    }
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
                ], style={"width": "55%", "display": "inline-block"}),
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

        print("min_value: %.2f" % min_value)
        print("max_value: %.2f" % max_value)

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

    def _get_clustering_tab(self):

        print("Getting clustering tab...")

        clustering_tab = html.Div(
            id="clustering_tab", className="tabcontent",
            children=[

                html.Div(
                    id="projection_div",
                    children=[

                        dcc.Graph(
                            id='projection',
                            figure=self._projection_figure,
                            config={'displaylogo': False}
                        )
                    ],
                    style={
                        "width": "50%",
                        "display": "inline-block",
                        "marginTop": 25}
                ),
                html.Div(
                    id="cluster_display_options", className="bordered_container",
                    children=[
                        html.H4("Plot Options", className="container_title"),
                        html.Div("Label Filters:",
                                 style={"width": "50%", "display": "inline-block"}),
                        html.Div("Projection:",
                                 style={"width": "25%", "display": "inline-block"}),
                        html.Div("X:",
                                 style={"width": "12.5%", "display": "inline-block"}),
                        html.Div("Y:",
                                 style={"width": "12.5%", "display": "inline-block"}),
                        html.Div(children=[
                            dcc.Dropdown(
                                id="cluster_filter_dropdown",
                                options=self._get_label_options(),
                                value=[],
                                multi=True
                            )],
                            style={
                                "width": "50%",
                                "display": "inline-block"
                            }
                        ),
                        html.Div(children=[
                            dcc.Dropdown(
                                id="transformation_method_dropdown",
                                options=self._get_transformation_method_options(),
                                value=self._get_default_transformation_method().name
                            )],
                            style={
                                "width": "25%",
                                "display": "inline-block",
                                "vertical-align": "top"
                            }
                        ),
                        html.Div(children=[
                            dcc.Dropdown(
                                id="x_axis_dropdown",
                                options=[],
                                value=0
                            )],
                            style={
                                "width": "12.5%",
                                "display": "inline-block",
                                "vertical-align": "top"
                            }
                        ),
                        html.Div(children=[
                            dcc.Dropdown(
                                id="y_axis_dropdown",
                                options=[],
                                value=1
                            )],
                            style={
                                "width": "12.5%",
                                "display": "inline-block",
                                "vertical-align": "top"
                            }
                        ),
                        dcc.Checklist(
                            id="union_checklist",
                            options=[
                                {
                                    'label': 'Union',
                                    'value': 'union'
                                }
                            ],
                            values=['']
                        ),
                        html.Div(
                            id="gene_filter", className="bordered_container",
                            children=[
                                html.H4("Gene Filter", className="container_title"),
                                html.Div(children=[
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
                                ),
                                html.Div(
                                    id="cell_value_range_slider_buttons",
                                    children=[
                                        html.Div("Min:",
                                                 style={
                                                     "width": "10%",
                                                     "display": "inline-block",
                                                     "marginLeft": 10
                                                 }),
                                        dcc.Input(
                                            id="manual_cell_value_range_slider_min",
                                            type="text", value=""
                                        ),
                                        html.Div("Max:",
                                                 style={
                                                     "width": "10%",
                                                     "display": "inline-block",
                                                     "marginLeft": 10
                                                 }),
                                        dcc.Input(
                                            id="manual_cell_value_range_slider_max",
                                            type="text", value=""
                                        ),
                                        html.Button(
                                            "Set",
                                            id="manual_cell_value_range_slider_button"
                                        )
                                    ],
                                    style={
                                        "marginTop": 20,
                                        "marginBottom": 20,
                                        "width": "100%"
                                    }
                                )
                            ]
                        ),
                        dcc.Checklist(
                            id="color_by_gene_count_checklist",
                            options=[
                                {
                                    'label': 'Color By Gene Count',
                                    'value': 'color_by_gene_count'
                                }
                            ],
                            values=['color_by_gene_count']
                        ),
                        html.Div(
                            id="label_type", className="bordered_container",
                            children=[
                                html.H4("Label Type", className="container_title"),
                                html.Div(children=[
                                    dcc.Dropdown(
                                        id="label_type_dropdown",
                                        options=Gene_Expression_Dataset_Plot._get_label_type_options(),
                                        value="L2FC",
                                        multi=False
                                    )
                                    ],
                                    style={
                                        "width": "50%"
                                    }
                                )
                            ]
                        )
                    ],
                    style={
                        "width": "45%",
                        "display": "inline-block",
                        "vertical-align": "top",
                        "marginLeft": "3%"
                    }
                ),
                html.Div(
                    id="label_management", className="bordered_container",
                    children=[
                        html.H4(children="Label Management",
                                className="container_title"),
                        html.Div(children=[
                            html.Label("Label currently highlighted cells:")
                        ]),
                        dcc.Input(
                            id="label_name",
                            type="text", value=""
                        ),
                        html.Button(
                            "+",
                            id="add_label_button"
                        ),
                        html.Div(children=[
                            html.Label("Current labels:")
                        ]),
                        html.Div(children=[
                            dcc.Dropdown(
                                id="manage_label_dropdown",
                                options=self._get_label_options(),
                                value=[])
                            ],
                            style={
                                "width": "50%",
                                "display": "inline-block"
                            }
                        ),
                        html.Div(children=[
                            html.Button("Delete", id="delete_label_button")
                            ]
                        )
                    ],
                    style={
                        "width": "50%",
                    }
                ),
                html.Div(
                    id="auto_clustering", className="bordered_container",
                    children=[
                        html.H4(children="Auto Clustering",
                                className="container_title"),
                        html.Div("Clustering method:"),
                        dcc.Dropdown(
                            id="cluster_method_dropdown",
                            options=self._get_cluster_method_options(),
                            value=self._get_default_cluster_method().name
                        ),
                        html.Div("# Clusters:"),
                        dcc.Input(
                            id="num_clusters",
                            type="text", value=""
                        ),
                        html.Div("Separate clusters by:"),
                        html.Div(
                            children=[
                                dcc.Dropdown(
                                    id="label_1_auto_cluster_dropdown",
                                    options=self._get_label_options(),
                                    value=[]
                                )
                            ],
                            style={
                                "width": "50%",
                                "display": "inline-block"
                            }
                        ),
                        html.Div(
                            children=[
                                dcc.Dropdown(
                                    id="label_2_auto_cluster_dropdown",
                                    options=self._get_label_options(),
                                    value=[]
                                )
                            ],
                            style={
                                "width": "50%",
                                "display": "inline-block"
                            }
                        ),
                        html.Div(children=[
                            html.Button(
                                "Auto cluster",
                                id="auto_cluster_button"
                            ),
                            html.Button(
                                "Show Auto Clusters",
                                id="show_auto_clusters_button"
                            )
                        ])
                    ],
                    style={
                        "width": "50%"
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

        return clustering_tab

    def start(self):

        external_stylesheets = ['https://codepen.io/chriddyp/pen/brPBPO.css']

        self._app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

        self._app.title = "SCRAP-viz"
        self._app.css.config.serve_locally = True
        self._app.scripts.config.serve_locally = True

        print("Creating projection figure...")

        self._projection_figure = self.get_projection_figure(
            Transformation_Method.TSNE)

        self._tabs = [
            self._get_clustering_tab(),
            self._get_differential_expression_tab()
        ]

        print("Creating data containers...")

        self._data_containers = [
            html.Div(id="labels", children=[], style={"display": "none"}),
            html.Div(id="de_data", children=[], style={"display": "none"}),
            html.Div(id="gene_range", children=[], style={"display": "none"}),
            html.Div(id="cell_color_values", children=[],
                     style={"display": "none"}),
            html.Div(id="unfiltered_cells", children=[],
                     style={"display": "none"}),
            html.Div(id="current_tab", children=[],
                     style={"display": "none"}),
            html.Div(id="projection_dimensions", children=[],
                     style={"display": "none"}),
            html.Div(id="projection_values", children=[],
                     style={"display": "none"}),
            html.Div(id="cell_clusters", children=[],
                     style={"display": "none"}),
            html.Div(id="dummy", children=[], style={"display": "none"})
        ]

        print("Creating layout...")

        self._app.layout = html.Div([

            html.Link(
                rel="stylesheet",
                href='/static/dropdown.css'
            ),

            html.Link(
                rel="stylesheet",
                href='/static/tabs.css'
            ),

            html.Link(
                rel="stylesheet",
                href='/static/containers.css'
            ),

            html.Link(
                rel="stylesheet",
                href='/static/plots.css'
            ),

            html.Div(id="tab_buttons", className="tab", children=[
                html.Button("Clustering", id="clustering_button",
                            className="tablinks"),
                html.Button("Differential Expression", className="tablinks",
                            id="differential_expression_button")
            ]),


            html.Div(id="tabs", children=self._tabs, style={"marginTop": 10}),
            html.Div(id="data", children=self._data_containers)
        ])

        print("Setting route...")

        @self._app.server.route('/static/<path:path>')
        def static_file(path):
            static_folder = os.path.join(os.getcwd(), 'static')
            return send_from_directory(static_folder, path)

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
             dash.dependencies.Input("add_label_button", "n_clicks"),
             dash.dependencies.Input("auto_cluster_button", "n_clicks")],
            [dash.dependencies.State("label_name", "value"),
             dash.dependencies.State("projection", "selectedData"),
             dash.dependencies.State("manage_label_dropdown", "value"),
             dash.dependencies.State("unfiltered_cells", "children"),
             dash.dependencies.State("num_clusters", "value"),
             dash.dependencies.State("transformation_method_dropdown", "value"),
             dash.dependencies.State("cluster_method_dropdown", "value"),
             dash.dependencies.State("label_1_auto_cluster_dropdown", "value"),
             dash.dependencies.State("label_2_auto_cluster_dropdown", "value")]
        )
        def label_added_or_deleted(delete_label_n_clicks, add_label_n_clicks,
                                   auto_cluster_n_clicks,
                                   label_name_to_add, selected_data,
                                   label_to_delete, unfiltered_cells,
                                   num_clusters, transformation_method,
                                   cluster_method,
                                   label_1_auto_cluster_dropdown,
                                   label_2_auto_cluster_dropdown):

            if verbose:
                print("label_added_or_deleted")

            current_labels = self._gene_expression_dataset.get_labels()

            if (delete_label_n_clicks is None or delete_label_n_clicks == 0)\
                    and (add_label_n_clicks is None or add_label_n_clicks == 0)\
                    and (auto_cluster_n_clicks is None or \
                                     auto_cluster_n_clicks == 0):
                return current_labels

            # Check if auto cluster was clicked
            if auto_cluster_n_clicks is not None and \
                    auto_cluster_n_clicks != self._n_clicks_auto_cluster:

                if transformation_method is not None:
                    transformation_method = \
                        Transformation_Method[
                            transformation_method]
                if cluster_method is not None:
                    cluster_method = \
                        Clustering_Method[
                            cluster_method]

                num_clusters = int(num_clusters)

                if (label_1_auto_cluster_dropdown != [] or \
                        label_2_auto_cluster_dropdown != []) and \
                        cluster_method != \
                        Clustering_Method.MAX_FEATURE:

                    self._gene_expression_dataset.get_matched_clusters(
                        label_1=label_1_auto_cluster_dropdown,
                        label_2=label_2_auto_cluster_dropdown,
                        num_clusters=num_clusters,
                        transformation_method=transformation_method,
                        clustering_method=cluster_method
                    )
                else:
                    self._gene_expression_dataset.auto_cluster(
                        num_clusters,
                        transformation_method=transformation_method,
                        clustering_method=cluster_method)
                self._gene_expression_dataset.save_labels()
                self._n_clicks_auto_cluster = auto_cluster_n_clicks
            # If n_clicks of delete button is the same, this is an add label
            elif delete_label_n_clicks == self._n_clicks_delete_label:
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

                if unfiltered_cells is not None:
                    cells_to_label = set()
                    for label in unfiltered_cells:
                        cells_to_label = \
                            cells_to_label.union(set(unfiltered_cells[label]))

                    cells_to_label = cells_to_label.intersection(selected_cells)
                else:
                    cells_to_label = selected_cells

                self._gene_expression_dataset.label_cells(label_name_to_add,
                                                          cells_to_label)
                self._gene_expression_dataset.save_labels()
            else:
                self._n_clicks_delete_label = delete_label_n_clicks
                self._gene_expression_dataset.delete_label(label_to_delete)
                self._gene_expression_dataset.save_labels()

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
            dash.dependencies.Output("clustering_button", "className"),
            [dash.dependencies.Input("current_tab", "children")])
        def current_tab_changed_update_clustering_button(current_tab):

            if verbose:
                print("current_tab_changed_update_clustering_button")

            if current_tab == "clustering":
                return "tablinks active"
            else:
                return "tablinks"

        @self._app.callback(
            dash.dependencies.Output("differential_expression_button", "className"),
            [dash.dependencies.Input("current_tab", "children")])
        def current_tab_changed_update_clustering_button(current_tab):

            if verbose:
                print("current_tab_changed_update_clustering_button")

            if current_tab == "differential_expression":
                return "tablinks active"
            else:
                return "tablinks"

        @self._app.callback(
            dash.dependencies.Output("current_tab", "children"),
            [dash.dependencies.Input("clustering_button", "n_clicks"),
             dash.dependencies.Input("differential_expression_button",
                                     "n_clicks")])
        def tabs_clicked_update_current_tab(
                clustering_button_clicks,
                differential_button_clicks):

            if verbose:
                print("tabs_clicked_update_current_tab")
            if clustering_button_clicks is not None and \
                    clustering_button_clicks > self._clustering_n_clicks:
                self._clustering_n_clicks += 1
                return "clustering"
            elif differential_button_clicks is not None and \
                    differential_button_clicks > \
                    self._differential_expression_n_clicks:
                self._differential_expression_n_clicks += 1
                return "differential_expression"

            return []

        @self._app.callback(
            dash.dependencies.Output("clustering_tab", "style"),
            [dash.dependencies.Input("current_tab", "children")])
        def tabs_clicked_update_clustering_tab(current_tab):
            if verbose:
                print("tabs_clicked_update_clustering_tab")

            if current_tab == "clustering":
                return {"display": "block"}
            else:
                return {"display": "none"}

        @self._app.callback(
            dash.dependencies.Output("de_tab", "style"),
            [dash.dependencies.Input("current_tab", "children")])
        def tabs_clicked_update_de_tab(current_tab):
            if verbose:
                print("tabs_clicked_update_de_tab")

            if current_tab == "differential_expression":
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
            [dash.dependencies.Input("cell_color_values", "children"),
             dash.dependencies.Input("manual_cell_value_range_slider_button", "n_clicks")],
            [dash.dependencies.State("manual_cell_value_range_slider_min", "value"),
             dash.dependencies.State("manual_cell_value_range_slider_max", "value")])
        def update_cell_value_range_slider(cell_color_values,
                                           n_clicks_slider,
                                           state_value_1,
                                           state_value_2):

            if verbose:
                print("update_cell_value_range_slider")

            if not cell_color_values:
                return []

            de_data = json.loads(cell_color_values)

            if not de_data:
                return []

            min_value = min(de_data.values())
            max_value = max(de_data.values())

            print("Min value: %.2f" % min_value)
            print("Max value: %.2f" % max_value)

            slider = [self.get_cell_value_range_slider(min_value, max_value)]

            if n_clicks_slider:
                if not self._n_clicks_slider or \
                        n_clicks_slider > self._n_clicks_slider:
                    print("Manually setting slider max")
                    slider[0].value = [float(state_value_1), float(state_value_2)]
                    self._n_clicks_slider = n_clicks_slider

            return slider

        @self._app.callback(
            dash.dependencies.Output("cell_color_values", "children"),
            [dash.dependencies.Input("gene_filter_dropdown", "value"),
             dash.dependencies.Input("label_type_dropdown", "value")])
        def update_cell_color_values(gene, label_type):

            if not gene:
                cell_read_counts = \
                    self._gene_expression_dataset.\
                        get_cell_total_transcript_counts()
                cell_read_count_ints = \
                    {cell_barcode: int(cell_read_counts[cell_index])
                        for cell_index, cell_barcode in enumerate(self._gene_expression_dataset._cell_transcript_counts.row_names)
                     }
                return json.dumps(cell_read_count_ints)

            print(label_type)

            if label_type == "Count":
                de = self._gene_expression_dataset.get_gene_counts(gene)
            elif label_type == "Normalized":
                de = self._gene_expression_dataset.get_gene_counts(
                    gene, normalized=True)
            else:
                de = self._gene_expression_dataset.get_cell_gene_differential(
                    gene)

            return json.dumps(de.to_dict())

        @self._app.callback(
            dash.dependencies.Output("manual_cell_value_range_slider_min", "value"),
            [dash.dependencies.Input("cell_value_range_slider", "value")])
        def set_manual_filter_input_min(input_values):

            if not input_values or len(input_values) < 2:
                return ""

            return str(input_values[0])

        @self._app.callback(
            dash.dependencies.Output("manual_cell_value_range_slider_max", "value"),
            [dash.dependencies.Input("cell_value_range_slider", "value")])
        def set_manual_filter_input_max(input_values):

            if not input_values or len(input_values) < 2:
                return ""

            return str(input_values[1])

        @self._app.callback(
            dash.dependencies.Output("unfiltered_cells", "children"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value"),
             dash.dependencies.Input("cell_value_range_slider", "value"),
             dash.dependencies.Input("union_checklist", "values")],
            [dash.dependencies.State("cell_color_values", "children")])
        def set_unfiltered_cells(selected_clusters, cell_value_range,
                                 union_checklist_values, cell_color_values):

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

            union = "union" in union_checklist_values

            if union:

                cells_by_cluster = {}

                if selected_clusters is None:
                    cells_by_cluster['All Cells'] = \
                        list(self._gene_expression_dataset.get_cells())
                    selected_clusters = []

                for label in selected_clusters:
                    cells_in_cluster = \
                        self._gene_expression_dataset.get_cells(label)

                    cells_to_remove = set()
                    if cell_value_range is not None and cells_gene_de is not None:
                        for cell in cells_in_cluster:
                            cell_color_value = cells_gene_de[cell]
                            if cell_color_value < cell_value_range[0] \
                                    or cell_color_value >= cell_value_range[1]:
                                cells_to_remove.add(cell)
                        unfiltered_cells = cells_in_cluster.difference(
                            cells_to_remove)
                    else:
                        unfiltered_cells = cells_in_cluster

                    cells_by_cluster[label] = list(unfiltered_cells)

            else:
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

                if selected_clusters is None:
                    selected_clusters = ['All Cells']

                cells_by_cluster = {" ".join(selected_clusters): list(unfiltered_cells)}

            return json.dumps(cells_by_cluster)

        @self._app.callback(
            dash.dependencies.Output("projection_dimensions", "children"),
            [dash.dependencies.Input("transformation_method_dropdown", "value")]
        )
        def update_projection_dimensions(transformation_method):

            if transformation_method is None:
                return []

            transformation_method = \
                Transformation_Method[
                    transformation_method
                ]

            num_dimensions = \
                len(self._gene_expression_dataset.get_cell_gene_expression(
                    transformation_method).columns)\

            dimension_options = [{"label": str(x+1), "value": x}
                                 for x in range(num_dimensions)]

            return json.dumps(dimension_options)

        @self._app.callback(
            dash.dependencies.Output("x_axis_dropdown", "options"),
            [dash.dependencies.Input("projection_dimensions", "children")]
        )
        def update_x_axis_dropdown(projection_dimensions):

            if projection_dimensions is None:
                return []

            return json.loads(projection_dimensions)

        @self._app.callback(
            dash.dependencies.Output("y_axis_dropdown", "options"),
            [dash.dependencies.Input("projection_dimensions", "children")]
        )
        def update_y_axis_dropdown(projection_dimensions):

            if projection_dimensions is None:
                return []

            return json.loads(projection_dimensions)

        @self._app.callback(
            dash.dependencies.Output("projection", "figure"),
            [dash.dependencies.Input("unfiltered_cells", "children"),
             dash.dependencies.Input("transformation_method_dropdown", "value"),
             dash.dependencies.Input("x_axis_dropdown", "value"),
             dash.dependencies.Input("y_axis_dropdown", "value"),
             dash.dependencies.Input("color_by_gene_count_checklist", "values")],
            [dash.dependencies.State("cell_color_values", "children")])
        def update_plot_from_filters(
                unfiltered_cells, transformation_method, x_axis_dimension,
                y_axis_dimension, color_by_gene_count_checklist_values,
                cell_color_values):

            if x_axis_dimension is None:
                x_axis_dimension = 0
            if y_axis_dimension is None:
                y_axis_dimension = 1

            if transformation_method is not None:
                transformation_method = \
                    Transformation_Method[
                        transformation_method]

            if verbose:
                print("update_plot_from_cluster_filter")

            color_by_gene_count = 'color_by_gene_count' in \
                                  color_by_gene_count_checklist_values

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

            proj_figure = self.get_projection_figure(
                transformation_method,
                highlighted_cells=unfiltered_cells,
                cell_color_values=cell_color_values,
                color_by_gene_count=color_by_gene_count,
                x_axis_dimension=x_axis_dimension,
                y_axis_dimension=y_axis_dimension)

            if verbose:
                print("update_plot_from_cluster_filter finished")

            return proj_figure

        @self._app.callback(
            dash.dependencies.Output("label_table", "children"),
            [dash.dependencies.Input("cluster_filter_dropdown", "value"),
             dash.dependencies.Input("labels", "children"),
             dash.dependencies.Input("union_checklist", "values")])
        def cluster_filter_updated(cluster_filter_values, clusters,
                                   union_checklist):

            if verbose:
                print("cluster_filter_updated")

            union = union_checklist is not None and "union" in union_checklist

            label_counts = self._gene_expression_dataset.get_label_counts(
                cluster_filter_values, union=union)

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
             dash.dependencies.State("subgroup_2_dropdown", "value"),
             dash.dependencies.State("de_across_dropdown", "value"),
             dash.dependencies.State("use_normalized_checklist", "values")])
        def new_de_clicked(n_clicks_de, subgroup_1_labels, subgroup_2_labels,
                           de_across_dropdown, use_normalized_checklist):
            if verbose:
                print("new_de_clicked")

            de_data = {"ready": False}

            use_normalized = "use_normalized" in use_normalized_checklist

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
                            differential_clusters=de_across_dropdown,
                            use_normalized=use_normalized)
                    self._de_stats = self._de_stats.sort_values(
                        ["p-value", "difference"],
                        ascending=[True, False]
                    )

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

                    # Ignore the difference column... wow this is hacky
                    if x >= 2:
                        column_to_sort_by = x
                    else:
                        column_to_sort_by = x - 1

                    if column_to_sort_by != self._column_to_sort:

                        self._column_to_sort = column_to_sort_by
                        self._de_start_range = 0

                        print("Sorting by column %i" % self._column_to_sort)
                        print(self._de_stats.columns)

                        if self._column_to_sort == 0:
                            self._de_stats["Log2 Change Abs"] = \
                                abs(self._de_stats["Log2 Change"])
                            self._de_stats = self._de_stats.sort_values(
                                "Log2 Change Abs", ascending=False)
                            self._de_stats = self._de_stats.drop(
                                "Log2 Change Abs",
                                axis=1)
                        elif self._column_to_sort == 1:
                            self._de_stats = self._de_stats.sort_values(
                                ["p-value", "difference"],
                                ascending=[True, False]
                            )
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

                    de = de.drop("difference", axis=1)

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
            [dash.dependencies.Input("de_table_graph", "clickData"),
             dash.dependencies.Input("de_plot", "clickData")]
        )
        def de_figure_clicked(table_click_data, plot_click_data):
            if verbose:
                print("de_figure_clicked")

            if table_click_data is not None:
                table_gene_index = table_click_data["points"][0]["y"] - 1 + \
                                   self._de_start_range
            else:
                table_gene_index = None

            if plot_click_data is not None:
                plot_gene_index = plot_click_data["points"][0]["pointNumber"]
            else:
                plot_gene_index = None

            if plot_gene_index is not None and \
                    plot_gene_index != self._eCDF_gene_index:
                gene_index = plot_gene_index
            elif table_gene_index is not None and \
                    table_gene_index != self._eCDF_gene_index:
                gene_index = table_gene_index
            else:
                gene_index = self._eCDF_gene_index

            if gene_index is None:
                return []

            return self._get_gene_eCDF(gene_index)

        @self._app.callback(
            dash.dependencies.Output("cluster_filter_dropdown", "value"),
            [dash.dependencies.Input("show_auto_clusters_button", "n_clicks")],
            [dash.dependencies.State("cluster_filter_dropdown", "options")]
        )
        def show_auto_clusters_button_clicked(n_clicks, cluster_filter_options):

            if n_clicks is None:
                return []

            auto_cluster_names = []

            for label_name in cluster_filter_options:
                if label_name["value"].find("Auto Cluster ") != -1:
                    auto_cluster_names.append(label_name["value"])

            return auto_cluster_names

        @self._app.callback(
            dash.dependencies.Output("union_checklist", "values"),
            [dash.dependencies.Input("show_auto_clusters_button", "n_clicks")],
            [dash.dependencies.State("union_checklist", "values")]
        )
        def show_auto_clusters_button_clicked_force_union(
                n_clicks, union_checklist):

            if n_clicks is None:
                return union_checklist

            return ["union"]

        @self._app.callback(
            dash.dependencies.Output("color_by_gene_count_checklist", "values"),
            [dash.dependencies.Input("show_auto_clusters_button", "n_clicks")],
            [dash.dependencies.State("color_by_gene_count_checklist", "values")]
        )
        def show_auto_clusters_button_clicked_force_cluster_coloring(
                n_clicks,
                color_by_gene_count_checklist):

            if n_clicks is None:
                return color_by_gene_count_checklist

            return []

        @self._app.callback(
            dash.dependencies.Output("cluster_filter_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")])
        def data_edited(labels):
            if verbose:
                print("data_edited")
            return self._get_label_options()

        @self._app.callback(
            dash.dependencies.Output("label_1_auto_cluster_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")])
        def labels_changed_update_label_1_dropdown(labels):
            if verbose:
                print("data_edited")
            return self._get_label_options()

        @self._app.callback(
            dash.dependencies.Output("label_2_auto_cluster_dropdown", "options"),
            [dash.dependencies.Input("labels", "children")])
        def labels_changed_update_label_2_dropdown(labels):
            if verbose:
                print("data_edited")
            return self._get_label_options()

        @self._app.callback(
            dash.dependencies.Output("dummy", "children"),
            [dash.dependencies.Input("export_de_button", "n_clicks")]
        )
        def export_de_clicked(n_clicks):

            if not n_clicks:
                return []

            if not self._subgroup_1_labels:
                return []

            if self._de_stats is None:
                return []

            file_name = "_".join(self._subgroup_1_labels)

            if self._subgroup_2_labels:
                file_name += "_vs_" + "_".join(self._subgroup_2_labels)

            file_name = file_name.replace(" ", "_")

            file_name += "_DE.csv"

            fileio.write_pandas_csv(self._de_stats, file_name)

            return []

        @self._app.callback(
            dash.dependencies.Output("gene_count_table", "children"),
            [dash.dependencies.Input("projection", "clickData")])
        def projection_clicked(click_data):

            if click_data is None:
                return []

            click_text = click_data["points"][0]["text"]
            cell_name = click_text[:click_text.find("<BR>")]

            gene_expression = self._gene_expression_dataset.\
                get_gene_expression_for_cell(cell_name)

            gene_expression = gene_expression[gene_expression["Count"] != 0]

            gene_expression["Log2 Change Abs"] = \
                abs(gene_expression["Log2 Fold Change"])

            gene_expression = gene_expression.sort_values(
                "Log2 Change Abs", ascending=False)

            gene_expression = gene_expression.drop("Log2 Change Abs", axis=1)

            gene_expression = gene_expression[0:50]

            gene_expression = gene_expression.\
                apply(lambda x: x.apply(lambda y: "%.3e" % y))

            return dcc.Graph(
                id="gene_counts_figure", figure=ff.create_table(
                    gene_expression, index=True, index_title="Gene"),
                config={'displayModeBar': False})

        self._app.run_server(port=self._port)

import pandas
from enum import Enum
import random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import os
import math
from scipy import stats
import csv
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.cluster import KMeans
import numpy
import scipy
import json
from sklearn import mixture
import urllib3 as urllib

Entrez_gene_summary_marker = " Entrez Gene Summary for "
gene_summary_marker = " Summary for "
gene_name_marker = "aliasMainName"


class Gene_Expression_Dataset:

    class Data_Mode(Enum):
        READ_COUNT = 0
        READS_PER_MILLION_TRANSCRIPTS = 1
        GENE_PROBABILITIES = 2

    class Normalization_Method(Enum):
        STD = 0
        ECDF = 1

    class Transformation_Method(Enum):
        PCA = 0
        TSNE = 1

    @staticmethod
    def get_sample_name(cell_name):
        return cell_name[:cell_name.find("_")]

    @staticmethod
    def read_pandas_csv(file_path):
        return pandas.read_csv(file_path, sep=",", header=0, index_col=0)

    @staticmethod
    def write_pandas_csv(data_frame, file_path):
        pandas.DataFrame(data_frame)\
            .to_csv(file_path, sep=',', encoding='utf-8', chunksize=1000)

    @staticmethod
    def get_cell_labels_file_path(dataset_path):
        return os.path.join(dataset_path, "labels.csv")

    @staticmethod
    def get_cell_transcript_counts_file_path(dataset_path):
        return os.path.join(dataset_path, "cell_transcript_counts.csv")

    @staticmethod
    def get_gene_counts_file_path(dataset_path):
        return os.path.join(dataset_path, "gene_counts.csv")

    @staticmethod
    def initialize_dataset(dataset_path, seed_matrices_file_path):

        if isinstance(seed_matrices_file_path, str):
            seed_matrix_file_path = seed_matrices_file_path
            data_frame = Gene_Expression_Dataset.\
                read_pandas_csv(seed_matrix_file_path)
        elif isinstance(seed_matrices_file_path, list):

            seed_matrix_file_path = seed_matrices_file_path[0]

            data_frames = [Gene_Expression_Dataset.
                read_pandas_csv(seed_matrix_file_path)]

            for file_index in range(1, len(seed_matrices_file_path)):
                seed_matrix_file_path = seed_matrices_file_path[file_index]
                data_frames.append(Gene_Expression_Dataset.
                    read_pandas_csv(seed_matrix_file_path))

            data_frame = pandas.concat(data_frames, axis=1)
        else:
            raise Exception("seed_matrices_file_path must be a path or list "\
                             "of paths")

        genes_below_threshold = \
            (data_frame.max(axis=1) < 1)

        data_frame = data_frame[~genes_below_threshold]

        Gene_Expression_Dataset.write_pandas_csv(
            data_frame, Gene_Expression_Dataset.get_gene_counts_file_path(
                dataset_path))

        samples = pandas.Series([
            Gene_Expression_Dataset.get_sample_name(cell_name)
            for cell_name in list(data_frame.columns)])

        samples = set(samples)

        label_cells = {}

        for sample_name in samples:
            label_cells[sample_name] = set()

        for cell_name in data_frame.columns:
            sample_name = Gene_Expression_Dataset.get_sample_name(cell_name)
            label_cells[sample_name].add(cell_name)

        Gene_Expression_Dataset.write_label_cells_to_file(
            label_cells,
            Gene_Expression_Dataset.get_cell_labels_file_path(dataset_path))

        cell_transcript_counts = data_frame.sum()
        cell_transcript_counts.name = "TOTAL_TRANSCRIPT_COUNT"

        Gene_Expression_Dataset.write_pandas_csv(
            cell_transcript_counts,
            Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
                dataset_path))

    @staticmethod
    def get_label_cells_from_file(file_path):

        cell_labels = {}

        cell_labels_file = open(file_path, "r")

        cell_labels_reader = csv.reader(cell_labels_file)

        for row in cell_labels_reader:

            label = row[0]
            cells = row[1:]

            cell_labels[label] = set(cells)

        cell_labels_file.close()

        return cell_labels

    @staticmethod
    def write_label_cells_to_file(label_cells, file_path):

        cell_labels_file = open(file_path, "w")
        writer = csv.writer(cell_labels_file)
        for label, cells in label_cells.items():
            label_array = [label]
            label_array.extend(list(cells))
            writer.writerow(label_array)
        cell_labels_file.close()

    def __init__(self, dataset_path, name=None, seed_matrix_file_path=None):

        if not os.path.exists(dataset_path):
            os.makedirs(dataset_path)
        elif os.path.isfile(dataset_path):
            raise Exception("Requested to initialize a dataset in '%s' \
                but this is a file, not a folder!" % dataset_path)

        if seed_matrix_file_path is not None:
            Gene_Expression_Dataset.initialize_dataset(
                dataset_path, seed_matrix_file_path)

        self._dataset_path = dataset_path
        self._gene_counts = None
        self._zero_genes = None
        self._normalized_gene_counts = None
        self._pca = None
        self._transformed = {}
        self._label_cells = {}
        self._cell_transcript_counts = None
        self._gene_means = None
        self._normalized_gene_means = None
        self._gene_cache = {}

        if name is None:
            self._load_dataset_from_path()
        else:
            self.load(name)

        self._gene_count_threshold = 0

    def reload(self):

        self.__init__(self._dataset_path)

    def get_labels(self):

        return list(self._label_cells.keys())

    def get_label_counts(self):

        label_counts = []

        for label, cells in self._label_cells:
            label_counts.append((label, len(self.get_cells(label))))

        return label_counts

    def filter_low_gene_counts(self, gene_count_threshold):

        self._gene_count_threshold = gene_count_threshold

        if gene_count_threshold <= 0:
            return

        genes_below_threshold =\
            (self._gene_counts.max(axis=1) < gene_count_threshold)

        self._gene_counts = self._gene_counts[~genes_below_threshold]
        if self._normalized_gene_counts is not None:
            self._normalized_gene_counts = \
                self._normalized_gene_counts[~genes_below_threshold]
            self._normalized_gene_means = \
                self._normalized_gene_means[~genes_below_threshold]
        self._gene_means = self._gene_means[~genes_below_threshold]

    def filter_low_transcript_cells(self, transcript_count_threshold):

        cells_below_threshold = \
            self._cell_transcript_counts["TOTAL_TRANSCRIPT_COUNT"] <\
            transcript_count_threshold

        cells_below_threshold = cells_below_threshold[cells_below_threshold]

        cell_names_below_threshold = set(cells_below_threshold.index)\
            .intersection(self._gene_counts.columns)

        self._gene_counts = self._gene_counts.drop(
            cell_names_below_threshold, axis=1)

        self.filter_low_gene_counts(self._gene_count_threshold)

        self._gene_means = self._gene_counts.mean(axis=1)

        if self._normalized_gene_counts is not None:
            self._normalized_gene_means = \
                self._normalized_gene_counts.mean(axis=1)

        self._cell_transcript_counts = self._cell_transcript_counts.drop(
            cell_names_below_threshold)

    def normalize_cells(
            self, data_mode=Data_Mode.READS_PER_MILLION_TRANSCRIPTS):

        if data_mode == self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS:

            self._gene_counts = self._gene_counts.div(
                self._cell_transcript_counts.loc[
                    self._gene_counts.columns]["TOTAL_TRANSCRIPT_COUNT"])

            self._gene_counts *= 1e6

        elif data_mode == self.Data_Mode.GENE_PROBABILITIES:

            for cell, gene_counts in self._gene_counts.iteritems():
                zero_counts = sum(gene_counts == 0)
                single_counts = sum(gene_counts == 1)
                total_count = self._cell_transcript_counts[
                    "TOTAL_TRANSCRIPT_COUNT"][cell]
                zero_probability = single_counts/total_count

                gene_counts = gene_counts/total_count*(1-zero_probability)
                gene_counts[gene_counts == 0] = zero_probability/zero_counts

                self._gene_counts[cell] = gene_counts

        self._gene_means = self._gene_counts.mean(axis=1)

    def label_cells(self, label, cells):

        if label not in self._label_cells:
            self._label_cells[label] = set()

        self._label_cells[label] = self._label_cells[label].union(cells)

    def delete_label(self, label):

        if label not in self._label_cells:
            return

        del self._label_cells[label]

    def transform(self, method=Transformation_Method.PCA, num_dimensions=2,
                  use_normalized=False):

        if use_normalized and self._normalized_gene_counts is None:
            use_normalized = False

        if method == self.Transformation_Method.PCA:

            self._pca = PCA(n_components=num_dimensions)

            if use_normalized:
                transformed = self._pca.fit_transform(
                    self._normalized_gene_counts.transpose())
            else:
                transformed = self._pca.fit_transform(
                    self._gene_counts.transpose())

            self._transformed[method] = pandas.DataFrame(transformed)

            self._transformed[method].columns = \
                ["PC_%i" % i for i in range(1, num_dimensions + 1)]

        elif method == self.Transformation_Method.TSNE:

            if self.Transformation_Method.PCA in self._transformed:
                transformed = TSNE(
                    verbose=True, perplexity=30, n_components=num_dimensions).\
                    fit_transform(
                        self._transformed[self.Transformation_Method.PCA])
            else:
                transformed = TSNE(
                    verbose=True, perplexity=30, n_components=num_dimensions).\
                    fit_transform(
                        self._normalized_gene_counts.transpose())

            self._transformed[self.Transformation_Method.TSNE] = \
                pandas.DataFrame(transformed)

            self._transformed[self.Transformation_Method.TSNE].columns = \
                ["tSNE_%i" % i for i in range(1, num_dimensions + 1)]

        if use_normalized:
            self._transformed[method].index = \
                self._normalized_gene_counts.columns
        else:
            self._transformed[method].index = \
                self._gene_counts.columns

    @property
    def num_cells(self):

        return self._gene_counts.shape[1]

    def get_cells(self, labels=None):

        if labels is None:
            return self._gene_counts.columns
        else:
            if isinstance(labels, str):
                return self._label_cells[labels].intersection(
                    self._gene_counts.columns)
            else:
                cells = self._gene_counts.columns

                for label in labels:
                    cells = cells.intersection(self._label_cells[label])

                return cells

    def get_genes(self):

        return list(self._gene_counts.index)

    def get_cell_gene_expression(self, transform=None):

        if transform is None:
            return self._gene_counts
        else:
            return self._transformed[transform]

    def get_cell_gene_expression_by_label(self, transform=None):

        label_cells = {}

        for label in self._label_cells:

            cell_names = self.get_cells(label)

            if transform is None:
                label_cells[label] = self._gene_counts[cell_names]
            else:
                label_cells[label] = \
                    self._transformed[transform].loc[cell_names]

        return label_cells

    def get_label_means(self, random_shuffle=False, is_median=False):

        if random_shuffle:
            label_counts = self.get_label_counts()
            labels = [label_count[0] for label_count in label_counts]
            label_weights = [label_count[1] for label_count in label_counts]

            label_cells = {}

            for label, _ in label_counts:
                label_cells[label] = set()

            for cell_name in self._gene_counts.columns:

                sample = random.choices(labels, weights=label_weights)
                label_cells[sample[0]].add(cell_name)

        else:
            label_cells = self._label_cells

        label_means = pandas.DataFrame()

        for label in label_cells:

            cell_names = self.get_cells(label)

            label_data_frame = self._gene_counts[list(cell_names)]

            if is_median:
                label_mean = label_data_frame.median(axis=1)
            else:
                label_mean = label_data_frame.mean(axis=1)

            label_mean.name = label
            label_means = label_means.append(label_mean)

        return label_means

    def get_cell_gene_differential(self, gene):

        cells_gene_count = self._gene_counts.loc[gene]

        non_zero_min = cells_gene_count[cells_gene_count > 0].min()

        cells_gene_count[cells_gene_count == 0] = non_zero_min / 2

        cells_gene_differential = cells_gene_count.apply(
            lambda x: math.log2(x / self._gene_means[gene]))

        return cells_gene_differential

    def compare_gene_expression(self, label_1, label_2=None,
                                use_normalized=True):

        if use_normalized and self._normalized_gene_counts is None:
            use_normalized = False

        gene_value_scores = {}

        label_1_cells = self.get_cells(label_1)

        if label_2 is not None:
            label_2_cells = self.get_cells(label_2)
        else:
            label_2_cells =\
                set(self._gene_counts.columns).difference(label_1_cells)

        if use_normalized:
            cell_gene_counts = self._normalized_gene_counts.copy()
        else:
            cell_gene_counts = self._gene_counts.copy()

        all_cells = label_1_cells.union(label_2_cells)
        zero_genes = self._zero_genes[list(all_cells)]
        genes_to_remove = set()

        for gene, gc in cell_gene_counts[list(all_cells)].iterrows():
            try:
                num_non_zero_genes = sum(~zero_genes.loc[gene])
            except TypeError as e:
                print(e)

            if num_non_zero_genes == 0:
                genes_to_remove.add(gene)

        cell_gene_counts = cell_gene_counts.drop(list(genes_to_remove))

        min_value = cell_gene_counts[cell_gene_counts > 0].min().min()

        for gene, gene_counts in cell_gene_counts.iterrows():

            sample_1_values = gene_counts[label_1_cells]
            sample_2_values = gene_counts[label_2_cells]

            sample_1_mean = sample_1_values.mean()
            sample_1_SD = sample_1_values.std()
            sample_2_mean = sample_2_values.mean()
            sample_2_SD = sample_2_values.std()

            if sample_1_mean == 0:
                sample_1_mean_for_log = min_value / 2
            else:
                sample_1_mean_for_log = sample_1_mean
            if sample_2_mean == 0:
                sample_2_mean_for_log = min_value / 2
            else:
                sample_2_mean_for_log = sample_2_mean

            log_2_fold_change = math.log2(sample_1_mean_for_log /
                                          sample_2_mean_for_log)

            _, p_value = stats.ks_2samp(sample_1_values, sample_2_values)

            gene_value_scores[gene] = (log_2_fold_change, p_value,
                                       sample_1_mean, sample_1_SD,
                                       sample_2_mean, sample_2_SD)

        df = pandas.DataFrame.from_dict(gene_value_scores, orient="index")

        df.columns = ["Log2 Change", "p-value", "Group 1 Mean",
                      "Group 1 SD", "Group 2 Mean", "Group 2 SD"]

        p_values = df["p-value"]

        _, p_values, _, _ = multipletests(p_values, method="bonferroni")

        df["p-value"] = p_values

        return df

    def get_gene_expression_for_cell(self, cell):

        cell_gene_counts = self._gene_counts[cell]
        cell_gene_de = cell_gene_counts.copy()

        non_zero_min = self._gene_counts[self._gene_counts > 0].min().min()
        cell_gene_de[cell_gene_de == 0] = non_zero_min / 2

        for gene, value in cell_gene_counts.items():
            cell_gene_de[gene] = \
                math.log2(cell_gene_de[gene]/self._gene_means[gene])

        cell_gene_expression = pandas.DataFrame(
            {
                "Count": cell_gene_counts,
                "Log2 Fold Change": cell_gene_de
            }
        )

        return cell_gene_expression

    def get_gene_counts(self, gene, filter_labels=None):

        cells = self.get_cells(filter_labels)

        return self._gene_counts.loc[gene][cells]

    def get_label_counts(self, filter_labels=None):

        label_counts = {}

        cells = self.get_cells(filter_labels)

        total_cells = len(cells)

        if total_cells == 0:
            return pandas.DataFrame()

        for label in self._label_cells.keys():

            num_cells = len(cells.intersection(self.get_cells(label)))
            cell_ratio = num_cells / total_cells

            if num_cells > 0:
                label_counts[label] = (num_cells, cell_ratio)

        df = pandas.DataFrame.from_dict(label_counts, orient="index")
        df.columns = ["# Cells", "Ratio"]

        return df

    def normalize_genes(self, method=Normalization_Method.STD):

        if method == self.Normalization_Method.STD:

            gene_sds = self._gene_counts.std(axis=1)

            self._normalized_gene_counts = self._gene_counts.copy().div(
                gene_sds, axis=0)

        elif method == self.Normalization_Method.ECDF:

            self._normalized_gene_counts = self._gene_counts.copy()
            gene_index = 0
            for gene, gene_counts in self._normalized_gene_counts.iterrows():

                value_counts = gene_counts.value_counts()
                eCDF = value_counts.sort_index().cumsum() * 1. / self.num_cells

                value_count_map = {}
                for i, j in eCDF.iteritems():
                    value_count_map[i] = j

                for cell, gene_count in gene_counts.iteritems():
                    gene_counts[cell] = value_count_map[gene_count]

                gene_index += 1

        self._normalized_gene_means = \
            self._normalized_gene_counts.mean(axis=1)

    def save_labels(self):

        Gene_Expression_Dataset.write_label_cells_to_file(
            self._label_cells, self._get_cell_labels_file_path())

    def load(self, name):

        gene_counts_path = os.path.join(self._dataset_path,
                                        "gene_counts_%s.csv" % name)

        self._gene_counts = Gene_Expression_Dataset.read_pandas_csv(
            gene_counts_path
        )

        normalized_path = os.path.join(self._dataset_path,
                                       "normalized_%s.csv" % name)

        if os.path.isfile(normalized_path):

            self._normalized_gene_counts = \
                Gene_Expression_Dataset.read_pandas_csv(normalized_path)

        for method_name, method in \
                self.Transformation_Method.__members__.items():

            file_name = "transformed_%s_%s.csv" % (method_name, name)

            file_path = os.path.join(self._dataset_path, file_name)

            if not os.path.isfile(file_path):
                continue

            self._transformed[method] = Gene_Expression_Dataset.read_pandas_csv(
                file_path
            )

        self._initialize_cache()

    def save(self, name):

        self.save_labels()

        Gene_Expression_Dataset.write_pandas_csv(
            self._gene_counts,
            os.path.join(self._dataset_path, "gene_counts_%s.csv" % name))

        Gene_Expression_Dataset.write_pandas_csv(
            self._normalized_gene_counts,
            os.path.join(self._dataset_path, "normalized_%s.csv" % name))

        for method_name, method in \
                self.Transformation_Method.__members__.items():

            if method not in self._transformed:
                continue

            file_name = "transformed_%s_%s.csv" % (method_name, name)

            file_path = os.path.join(self._dataset_path, file_name)

            Gene_Expression_Dataset.write_pandas_csv(
                self._transformed[method],
                file_path)

    def _save_gene_cache(self):

        gene_cache_path = os.path.join(self._dataset_path, "gene_cache.json")

        with open(gene_cache_path, "w") as gene_cache_file:
            json.dump(self._gene_cache, gene_cache_file, indent=4)

    def get_cell_transcript_counts(self):
        return self._cell_transcript_counts

    def get_matched_clusters(self, label_1, label_2=None, num_clusters=20,
                             transformation_method=Transformation_Method.PCA):

        label_1_cells = list(self.get_cells(label_1))
        label_2_cells = list(self.get_cells(label_2))

        label_1_transformed = \
            self._transformed[transformation_method].loc[label_1_cells]
        label_2_transformed = \
            self._transformed[transformation_method].loc[label_2_cells]

        label_1_k_means = KMeans(n_clusters=num_clusters, random_state=0)
        label_2_k_means = KMeans(n_clusters=num_clusters, random_state=0)

        label_1_fitted = label_1_k_means.fit(label_1_transformed)
        label_2_fitted = label_2_k_means.fit(label_2_transformed)

        label_1_cluster_centers = label_1_fitted.cluster_centers_
        label_2_cluster_centers = label_2_fitted.cluster_centers_

        # label_1_mixture = mixture.GaussianMixture(n_components=num_clusters)
        # label_2_mixture = mixture.GaussianMixture(n_components=num_clusters)
        #
        # label_1_fitted = label_1_mixture.fit(label_1_transformed)
        # label_2_fitted = label_2_mixture.fit(label_2_transformed)
        #
        # label_1_cluster_centers = label_1_fitted.means_
        # label_2_cluster_centers = label_2_fitted.means_

        cluster_distances = numpy.empty((num_clusters, num_clusters))

        for label_1_cluster_index in range(num_clusters):
            label_1_cluster_center = \
                label_1_cluster_centers[label_1_cluster_index]
            for label_2_cluster_index in range(0, num_clusters):
                label_2_cluster_center = \
                    label_2_cluster_centers[label_2_cluster_index]
                distance = numpy.linalg.norm(
                    label_1_cluster_center - label_2_cluster_center)
                cluster_distances[label_1_cluster_index][label_2_cluster_index]\
                    = distance

        cluster_assignments = scipy.optimize.linear_sum_assignment(
            cluster_distances)

        label_cells = {}

        label_1_clusters = label_1_fitted.labels_
        label_2_clusters = label_2_fitted.labels_

        # label_1_clusters = label_1_fitted.predict(label_1_transformed)
        # label_2_clusters = label_2_fitted.predict(label_2_transformed)

        for cluster_index in range(num_clusters):
            label_2_cluster_index = cluster_assignments[1][cluster_index]

            label_1_cluster_cells = label_1_transformed[
                label_1_clusters == cluster_index
            ]

            label_2_cluster_cells = label_2_transformed[
                label_2_clusters == label_2_cluster_index
            ]

            cluster_cells = list(label_1_cluster_cells.index)
            cluster_cells.extend(list(label_2_cluster_cells.index))

            label_cells["Cluster %i" % cluster_index] = cluster_cells

        return label_cells

    def get_gene_summaries(self, genes):

        gene_name_descriptions = []

        http = urllib.PoolManager()

        added_gene = False

        for gene in genes:

            if gene in self._gene_cache and self._gene_cache[gene] != \
                    (gene, "No response"):
                gene_name_descriptions.append(self._gene_cache[gene])
                continue

            url = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" % gene
            response = http.request("GET", url)

            if response.status != 200:
                gene_name_descriptions.append((gene, "No response"))
                continue

            added_gene = True

            string = response.data.decode("UTF-8")

            if string.find("was not found") != -1:
                self._gene_cache[gene] = (gene, "Gene not found")
                gene_name_descriptions.append(self._gene_cache[gene])
                continue

            gene_name_marker_location = string.find(gene_name_marker)

            if gene_name_marker_location == -1:
                gene_name = gene
            else:
                name_start_index = string.find(">", gene_name_marker_location)
                name_end_index = string.find("<", name_start_index)

                gene_name = string[name_start_index + 1:name_end_index]

            Entrez_gene_summary_location = string.find(
                Entrez_gene_summary_marker)

            if Entrez_gene_summary_location == -1:
                gene_summary_marker_location = string.find(gene_summary_marker)

                if gene_summary_marker_location == -1:
                    self._gene_cache[gene] = \
                        (gene_name, "Gene summary not found")
                    gene_name_descriptions.append(self._gene_cache[gene])
                    continue
            else:
                gene_summary_marker_location = Entrez_gene_summary_location

            summary_start_index = string.find("<p>",
                                              gene_summary_marker_location)

            summary_end_index = string.find("</p>", summary_start_index)

            gene_description = string[summary_start_index + 3:summary_end_index]
            gene_description = " ".join(gene_description.split())

            self._gene_cache[gene] = (gene_name, gene_description)
            gene_name_descriptions.append(self._gene_cache[gene])

        if added_gene:
            self._save_gene_cache()

        return gene_name_descriptions

    def _load_dataset_from_path(self):

        gene_count_file_path = os.path.join(dataset_path, "gene_counts.csv")

        self._gene_counts = Gene_Expression_Dataset.read_pandas_csv(
            gene_count_file_path)

        self._initialize_cache()

    def _initialize_cache(self):

        self._zero_genes = self._gene_counts.apply(lambda x: x == 0)

        self._gene_means = self._gene_counts.mean(axis=1)

        if self._normalized_gene_counts is not None:

            self._normalized_gene_means = self._normalized_gene_counts.mean(
                axis=1
            )

        self._cell_transcript_counts = Gene_Expression_Dataset.read_pandas_csv(
            self.get_cell_transcript_counts_file_path(self._dataset_path))

        self._label_cells = Gene_Expression_Dataset.get_label_cells_from_file(
            self.get_cell_labels_file_path(self._dataset_path))

        gene_cache_path = os.path.join(self._dataset_path, "gene_cache.json")

        if os.path.isfile(gene_cache_path):
            with open(gene_cache_path, "r") as gene_cache_file:
                self._gene_cache = json.load(gene_cache_file)

    def _get_cell_labels_file_path(self):
        return Gene_Expression_Dataset.get_cell_labels_file_path(
            self._dataset_path)

    def _get_cell_transcript_counts_file_path(self):
        return Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
            self._dataset_path)

    def _get_gene_counts_file_path(self):
        return Gene_Expression_Dataset.get_gene_counts_file_path(
            self._dataset_path)

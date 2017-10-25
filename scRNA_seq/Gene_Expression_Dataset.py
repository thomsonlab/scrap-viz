import pandas
from enum import Enum
from . import plotting
import random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import os
import math
from scipy import stats
import csv


class Gene_Expression_Dataset:

    class Data_Mode(Enum):
        READ_COUNT = 0
        READS_PER_MILLION_TRANSCRIPTS = 1

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
    def initialize_dataset(dataset_path, seed_matrix_file_path):

        data_frame = Gene_Expression_Dataset.\
            read_pandas_csv(seed_matrix_file_path)

        Gene_Expression_Dataset.write_pandas_csv(
            data_frame, Gene_Expression_Dataset.get_gene_counts_file_path(
                dataset_path))

        samples = pandas.Series([
            Gene_Expression_Dataset.get_sample_name(cell_name)
            for cell_name in list(data_frame.columns)])

        samples = set(samples)

        label_cells = {}

        if len(samples) > 1:
            for sample_name in samples:
                    label_cells[sample_name] = set()

            for cell_name in data_frame.columns:
                sample_name = Gene_Expression_Dataset.get_sample_name(cell_name)
                label_cells[sample_name].add(cell_name)

        Gene_Expression_Dataset.write_label_cells_to_file(
            label_cells,
            Gene_Expression_Dataset.get_cell_labels_file_path(dataset_path))

        cell_transcript_counts = data_frame.sum()

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

    def __init__(self, dataset_path, seed_matrix_file_path=None):

        if not os.path.exists(dataset_path):
            os.makedirs(dataset_path)
        elif os.path.isfile(dataset_path):
            raise Exception("Requested to initialize a dataset in '%s' \
                but this is a file, not a folder!" % dataset_path)

        if seed_matrix_file_path is not None:
            Gene_Expression_Dataset.initialize_dataset(
                dataset_path, seed_matrix_file_path)

        self._dataset_path = dataset_path
        self._data_frame = None
        self._pca = None
        self._transformed = {}
        self._label_cells = {}
        self._cell_transcript_counts = None

        self._load_dataset_from_path(dataset_path)

        self._data_mode = self.Data_Mode.READ_COUNT
        self._gene_count_threshold = 0

    def reload(self):

        self.__init__(self._dataset_path)

    def get_label_counts(self):

        label_counts = []

        for label, cells in self._label_cells:
            label_counts.append((label, len(cells)))

        return label_counts

    def plot_cell_transcript_count_eCDF(self):

        cell_transcript_counts = self._data_frame.sum(axis=0)
        value_counts = cell_transcript_counts.value_counts()
        eCDF = value_counts.sort_index().cumsum() * 1. / \
               self.num_cells

        plotting.plot_eCDF(eCDF)

    def filter_low_gene_counts(self, gene_count_threshold):

        self._gene_count_threshold = gene_count_threshold

        if gene_count_threshold <= 0:
            return

        genes_below_threshold =\
            (self._data_frame.max(axis=1) < gene_count_threshold)

        self._data_frame = self._data_frame[~genes_below_threshold]

    def filter_low_transcript_cells(self, transcript_count_threshold):

        cells_below_threshold = self._cell_transcript_counts <\
            transcript_count_threshold

        cells_below_threshold = cells_below_threshold[cells_below_threshold]

        cell_names_below_threshold = set(cells_below_threshold.index)

        self._data_frame = self._data_frame.drop(cell_names_below_threshold,
                                                 axis=1)

        for label in self._label_cells:
            self._label_cells[label] = self._label_cells[label].difference(
                cell_names_below_threshold)

        self.filter_low_gene_counts(self._gene_count_threshold)

        self._cell_transcript_counts = self._cell_transcript_counts.drop(
            cell_names_below_threshold)

    def convert_to_reads_per_million_transcripts(self):

        if self._data_mode == self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS:
            return

        self._data_frame = self._data_frame.div(self._cell_transcript_counts)

        self._data_frame *= 1e6

        self._data_mode = self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS

    def label_cells(self, label, cells):

        if label not in self._label_cells:
            self._label_cells[label] = set()

        self._label_cells[label] = self._label_cells[label].union(cells)

    def transform(self, method=Transformation_Method.PCA, num_dimensions=2):

        if method == self.Transformation_Method.PCA:

            self._pca = PCA(n_components=num_dimensions)

            transformed = self._pca.fit_transform(self._data_frame.transpose())

            self._transformed[method] = pandas.DataFrame(transformed)

            self._transformed[method].columns = \
                ["PC_%i" % i for i in range(1, num_dimensions + 1)]

        elif method == self.Transformation_Method.TSNE:

            if self.Transformation_Method.PCA in self._transformed:
                transformed = TSNE(verbose=True).fit_transform(
                    self._transformed[self.Transformation_Method.PCA])
            else:
                transformed = TSNE(verbose=True).fit_transform(
                    self._data_frame.transpose())

            self._transformed[self.Transformation_Method.TSNE] = \
                pandas.DataFrame(transformed)

            self._transformed[self.Transformation_Method.TSNE].columns = \
                ["tSNE_%i" % i for i in range(1, num_dimensions + 1)]

        self._transformed[method].index = self._data_frame.columns

    @property
    def num_cells(self):

        return self._data_frame.shape[1]

    def get_cell_gene_expression_by_label(self, transform=None):

        label_cells = {}

        for label in self._label_cells:
            if transform is None:
                label_cells[label] = self._data_frame[self._label_cells[label]]
            else:
                label_cells[label] = \
                    self._transformed[transform].loc[self._label_cells[label]]

        return label_cells

    def get_label_means(self, random_shuffle=False, is_median=False):

        if random_shuffle:
            label_counts = self.get_label_counts()
            labels = [label_count[0] for label_count in label_counts]
            label_weights = [label_count[1] for label_count in label_counts]

            label_cells = {}

            for label, _ in label_counts:
                label_cells[label] = set()

            for cell_name in self._data_frame.columns:

                sample = random.choices(labels, weights=label_weights)
                label_cells[sample[0]].add(cell_name)

        else:
            label_cells = self._label_cells

        label_means = pandas.DataFrame()

        for label in label_cells:

            label_data_frame = self._data_frame[list(label_cells[label])]

            if is_median:
                label_mean = label_data_frame.median(axis=1)
            else:
                label_mean = label_data_frame.mean(axis=1)

            label_mean.name = label
            label_means = label_means.append(label_mean)

        return label_means

    def compare_gene_expression(self, label_1, label_2=None):

        gene_value_scores = {}

        min_value = self._data_frame[self._data_frame > 0].min().min()

        if isinstance(label_1, str):
            label_1_cells = self._label_cells[label_1]
        else:
            label_1_cells = set(self._data_frame.columns)

            for label in label_1:
                label_1_cells = \
                    label_1_cells.intersection(self._label_cells[label])

        if label_2 is None:
            label_2_cells = set(self._data_frame.columns).difference(
                label_1_cells)
        elif isinstance(label_2, str):
            label_2_cells = self._label_cells[label_2]
        else:
            label_2_cells = set(self._data_frame.columns)

            for label in label_2:
                label_2_cells = \
                    label_2_cells.intersection(self._label_cells[label])

        for gene, gene_counts in self._data_frame.iterrows():

            sample_1_values = gene_counts[label_1_cells]
            sample_2_values = gene_counts[label_2_cells]

            sample_1_mean = sample_1_values.mean()
            sample_2_mean = sample_2_values.mean()

            if sample_1_mean == 0:
                sample_1_mean = min_value / 2
            if sample_2_mean == 0:
                sample_2_mean = min_value / 2

            log_2_fold_change = math.log2(sample_1_mean / sample_2_mean)

            _, p_value = stats.ttest_ind(sample_1_values, sample_2_values)

            gene_value_scores[gene] = (log_2_fold_change, p_value)

        return gene_value_scores

    def normalize(self, method=Normalization_Method.STD):

        self.convert_to_reads_per_million_transcripts()

        if method == self.Normalization_Method.STD:

            gene_stds = self._data_frame.std(axis=1)

            self._data_frame = self._data_frame.div(gene_stds, axis=0)

        elif method == self.Normalization_Method.ECDF:

            gene_index = 0
            for gene, gene_counts in self._data_frame.iterrows():

                value_counts = gene_counts.value_counts()
                eCDF = value_counts.sort_index().cumsum() * 1. / self.num_cells

                map = {}
                for i, j in eCDF.iteritems():
                    map[i] = j

                for cell, gene_count in gene_counts.iteritems():
                    gene_counts[cell] = map[gene_count]

                gene_index += 1

    def save(self):

        Gene_Expression_Dataset.write_pandas_csv(
            self._data_frame, self._get_gene_counts_file_path())

        for method, data_frame in self._transformed.items():

            method_name = \
                self.Transformation_Method(method).name

            file_path = os.path.join(self._dataset_path, "%s.csv" % method_name)

            Gene_Expression_Dataset.write_pandas_csv(data_frame, file_path)

        Gene_Expression_Dataset.write_label_cells_to_file(
            self._label_cells, self._get_cell_labels_file_path())

    def _load_dataset_from_path(self, dataset_path):

        gene_count_file_path = os.path.join(dataset_path, "gene_counts.csv")

        self._data_frame = Gene_Expression_Dataset.read_pandas_csv(
            gene_count_file_path)

        self._cell_transcript_counts = Gene_Expression_Dataset.read_pandas_csv(
            self.get_cell_transcript_counts_file_path(self._dataset_path))

        self._label_cells = Gene_Expression_Dataset.get_label_cells_from_file(
            self.get_cell_labels_file_path(self._dataset_path))

        for method_name, method in \
                self.Transformation_Method.__members__.items():

            file_name = "%s.csv" % method_name

            file_path = os.path.join(dataset_path, file_name)

            if os.path.isfile(file_path):
                self._transformed[method] = pandas.read_csv(
                    file_path, sep=",", header=0, index_col=0)

    def _get_cell_labels_file_path(self):
        return Gene_Expression_Dataset.get_cell_labels_file_path(
            self._dataset_path)

    def _get_cell_transcript_counts_file_path(self):
        return Gene_Expression_Dataset.get_cell_transcript_counts_file_path(
            self._dataset_path)

    def _get_gene_counts_file_path(self):
        return Gene_Expression_Dataset.get_gene_counts_file_path(
            self._dataset_path)

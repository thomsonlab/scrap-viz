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
from statsmodels.sandbox.stats.multicomp import multipletests


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
    def initialize_dataset(dataset_path, seed_matrix_file_path):

        data_frame = Gene_Expression_Dataset.\
            read_pandas_csv(seed_matrix_file_path)

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
        self._gene_counts = None
        self._zero_genes = None
        self._normalized_gene_counts = None
        self._pca = None
        self._transformed = {}
        self._label_cells = {}
        self._cell_transcript_counts = None

        self._load_dataset_from_path(dataset_path)

        self._data_mode = self.Data_Mode.READ_COUNT
        self._gene_count_threshold = 0

    def reload(self):

        self.__init__(self._dataset_path)

    def get_labels(self):

        return list(self._label_cells.keys())

    def get_label_counts(self):

        label_counts = []

        for label, cells in self._label_cells:
            label_counts.append((label, len(cells)))

        return label_counts

    def plot_cell_transcript_count_eCDF(self):

        cell_transcript_counts = self._gene_counts.sum(axis=0)
        value_counts = cell_transcript_counts.value_counts()
        eCDF = value_counts.sort_index().cumsum() * 1. / \
               self.num_cells

        plotting.plot_eCDF(eCDF)

    def filter_low_gene_counts(self, gene_count_threshold):

        self._gene_count_threshold = gene_count_threshold

        if gene_count_threshold <= 0:
            return

        genes_below_threshold =\
            (self._gene_counts.max(axis=1) < gene_count_threshold)

        self._gene_counts = self._gene_counts[~genes_below_threshold]

    def filter_low_transcript_cells(self, transcript_count_threshold):

        cells_below_threshold = \
            self._cell_transcript_counts["TOTAL_TRANSCRIPT_COUNT"] <\
            transcript_count_threshold

        cells_below_threshold = cells_below_threshold[cells_below_threshold]

        cell_names_below_threshold = set(cells_below_threshold.index)\
            .intersection(self._gene_counts.columns)

        self._gene_counts = self._gene_counts.drop(
            cell_names_below_threshold, axis=1)

        for label in self._label_cells:
            self._label_cells[label] = self._label_cells[label].difference(
                cell_names_below_threshold)

        self.filter_low_gene_counts(self._gene_count_threshold)

        self._cell_transcript_counts = self._cell_transcript_counts.drop(
            cell_names_below_threshold)

    def normalize_cells(
            self, data_mode=Data_Mode.READS_PER_MILLION_TRANSCRIPTS):

        if self._data_mode != self.Data_Mode.READ_COUNT:
            raise Exception("Cells already normalized!")

        if data_mode == self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS:

            self._gene_counts = self._gene_counts.div(
                self._cell_transcript_counts.loc[
                    self._gene_counts.columns]["TOTAL_TRANSCRIPT_COUNT"])

            self._gene_counts *= 1e6

        elif data_mode == self.Data_Mode.GENE_PROBABILITIES:

            self._zero_genes = self._gene_counts.apply(lambda x: x == 0)

            for cell, gene_counts in self._gene_counts.iteritems():
                value_counts = gene_counts.value_counts()
                zero_counts = value_counts[0]
                single_counts = value_counts[1]
                total_count = self._cell_transcript_counts[
                    "TOTAL_TRANSCRIPT_COUNT"][cell]
                zero_probability = single_counts/total_count

                gene_counts = gene_counts/total_count/(1+zero_probability)
                gene_counts[gene_counts == 0] = zero_probability/zero_counts

                self._gene_counts[cell] = gene_counts

        self._data_mode = data_mode

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
                transformed = TSNE(verbose=True, perplexity=30).fit_transform(
                    self._transformed[self.Transformation_Method.PCA])
            else:
                transformed = TSNE(verbose=True, perplexity=30).fit_transform(
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
                return self._label_cells[labels]
            else:
                cells = self._gene_counts.columns

                for label in labels:
                    cells = cells.intersection(self._label_cells[label])

                return cells

    def get_cell_gene_expression(self, transform=None):

        if transform is None:
            return self._gene_counts
        else:
            return self._transformed[transform]

    def get_cell_gene_expression_by_label(self, transform=None):

        label_cells = {}

        for label in self._label_cells:
            if transform is None:
                label_cells[label] = self._gene_counts[self._label_cells[label]]
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

            for cell_name in self._gene_counts.columns:

                sample = random.choices(labels, weights=label_weights)
                label_cells[sample[0]].add(cell_name)

        else:
            label_cells = self._label_cells

        label_means = pandas.DataFrame()

        for label in label_cells:

            label_data_frame = self._gene_counts[list(label_cells[label])]

            if is_median:
                label_mean = label_data_frame.median(axis=1)
            else:
                label_mean = label_data_frame.mean(axis=1)

            label_mean.name = label
            label_means = label_means.append(label_mean)

        return label_means

    def compare_gene_expression(self, label_1, label_2=None,
                                use_normalized=True):

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

            if self._data_mode == self.Data_Mode.GENE_PROBABILITIES:
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
                sample_1_mean = min_value / 2
            if sample_2_mean == 0:
                sample_2_mean = min_value / 2

            log_2_fold_change = math.log2(sample_1_mean / sample_2_mean)

            _, p_value = stats.ranksums(sample_1_values, sample_2_values)

            gene_value_scores[gene] = (log_2_fold_change, p_value,
                                       sample_1_mean, sample_1_SD,
                                       sample_2_mean, sample_2_SD)

        df = pandas.DataFrame.from_dict(gene_value_scores, orient="index")

        df.columns = ["Log2 Change", "p-value", "Group 1 Mean",
                      "Group 1 SD", "Group 2 Mean", "Group 2 SD"]

        p_values = df["p-value"]
        #
        # _, p_values, _, _ = multipletests(p_values)
        #
        # df["p-value"] = p_values

        return df

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

            num_cells = len(self._label_cells[label].intersection(cells))
            cell_ratio = num_cells / total_cells

            if num_cells > 0:
                label_counts[label] = (num_cells, cell_ratio)

        df = pandas.DataFrame.from_dict(label_counts, orient="index")
        df.columns = ["# Cells", "Ratio"]

        return df

    def normalize_genes(self, method=Normalization_Method.STD):

        if method == self.Normalization_Method.STD:

            gene_stds = self._gene_counts.std(axis=1)

            self._normalized_gene_counts = self._gene_counts.copy().div(
                gene_stds, axis=0)

        elif method == self.Normalization_Method.ECDF:

            self._normalized_gene_counts = self._gene_counts.copy()
            gene_index = 0
            for gene, gene_counts in self._normalized_gene_counts.iterrows():

                value_counts = gene_counts.value_counts()
                eCDF = value_counts.sort_index().cumsum() * 1. / self.num_cells

                map = {}
                for i, j in eCDF.iteritems():
                    map[i] = j

                for cell, gene_count in gene_counts.iteritems():
                    gene_counts[cell] = map[gene_count]

                gene_index += 1

    def save(self):

        # Gene_Expression_Dataset.write_pandas_csv(
        #     self._gene_counts, self._get_gene_counts_file_path())
        #
        # for method, data_frame in self._transformed.items():
        #
        #     method_name = \
        #         self.Transformation_Method(method).name
        #
        #     file_path = os.path.join(self._dataset_path, "%s.csv" % method_name)
        #
        #     Gene_Expression_Dataset.write_pandas_csv(data_frame, file_path)

        Gene_Expression_Dataset.write_label_cells_to_file(
            self._label_cells, self._get_cell_labels_file_path())

    def _load_dataset_from_path(self, dataset_path):

        gene_count_file_path = os.path.join(dataset_path, "gene_counts.csv")

        self._gene_counts = Gene_Expression_Dataset.read_pandas_csv(
            gene_count_file_path)

        # duplicate_aggregator = {}
        # for cell in self._gene_counts.columns:
        #     duplicate_aggregator[cell] = sum

        # self._gene_counts = self._gene_counts.groupby(level=0).sum()

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

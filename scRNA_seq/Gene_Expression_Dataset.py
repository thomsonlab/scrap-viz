import pandas
from enum import Enum
from . import plotting
import time
import random


class Gene_Expression_Dataset:

    class Data_Mode(Enum):
        READ_COUNT = 0
        READS_PER_MILLION_TRANSCRIPTS = 1

    class Normalization_Method(Enum):
        STD = 0
        ECDF = 1

    def get_sample_name(cell_name):
        return cell_name[:cell_name.find("_")]

    def get_random_sample(self, cell_name):

        if cell_name not in self._random_sample_map:
            self._random_sample_map[cell_name] = random.choices(
                self._sample_names, weights=self._sample_weights)

        return self._random_sample_map[cell_name]

    def __init__(self, file_path=None):

        if file_path is not None:
            self._load_dataset_from_file(file_path)
        else:
            self._data_frame = pandas.DataFrame()

        self._data_mode = self.Data_Mode.READ_COUNT

    def get_sample_counts(self):

        sample_names = pandas.Series([
            cell_name[:cell_name.find("_")]
            for cell_name in list(self._data_frame.columns)])

        return sample_names.value_counts()

    def plot_cell_transcript_count_eCDF(self):

        cell_transcript_counts = self._data_frame.sum(axis=0)
        value_counts = cell_transcript_counts.value_counts()
        eCDF = value_counts.sort_index().cumsum() * 1. / \
               self.num_cells

        plotting.plot_eCDF(eCDF)

    def filter_low_gene_counts(self, gene_count_threshold):

        if gene_count_threshold <= 0:
            return

        genes_below_threshold =\
            (self._data_frame.max(axis=1) < gene_count_threshold)

        self._data_frame = self._data_frame[~genes_below_threshold]

    def filter_low_transcript_cells(self, transcript_count_threshold):


    def convert_to_reads_per_million_transcripts(self):

        if self._data_mode == self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS:
            return

        self._data_frame /= self._cell_transcript_counts

        self._data_frame *= 1e6

        self._data_mode = self.Data_Mode.READS_PER_MILLION_TRANSCRIPTS

    @property
    def num_cells(self):

        return self._data_frame.shape[1]

    def get_sample_means(self, random_shuffle=False, is_median=False):

        if random_shuffle:
            sample_counts = self.get_sample_counts()
            self._sample_names = list(sample_counts.index)
            self._sample_weights = list(sample_counts.values)
            self._random_sample_map = {}

            group_by_function = self.get_random_sample
        else:
            group_by_function = Gene_Expression_Dataset.get_sample_name

        sample_data_frames = self._data_frame.groupby(group_by_function, axis=1)

        sample_means = pandas.DataFrame()

        for sample, data_frame in sample_data_frames:

            if is_median:
                sample_mean = data_frame.median(axis=1)
            else:
                sample_mean = data_frame.mean(axis=1)
            sample_mean.name = sample
            sample_means = sample_means.append(sample_mean)

        return sample_means

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

    def save(self, file_path):

        self._data_frame = self._data_frame.append(self._cell_transcript_counts)
        self._data_frame.to_csv(
            file_path, sep=',', encoding='utf-8', chunksize=1000)
        self._data_frame = self._data_frame.drop(["_TOTAL_TRANSCRIPT_COUNT"])

    def _load_dataset_from_file(self, file_path):

        self._data_frame = pandas.read_csv(
            file_path, sep=",", header=0, index_col=0)

        if "_TOTAL_TRANSCRIPT_COUNT" not in self._data_frame.index:
            self._initialize_dataset()
        else:
            self._cell_transcript_counts =\
                self._data_frame.loc["_TOTAL_TRANSCRIPT_COUNT"]
            self._data_frame = self._data_frame.drop(
                ["_TOTAL_TRANSCRIPT_COUNT"])

    def _initialize_dataset(self):

        self._cell_transcript_counts = pandas.DataFrame(
            {"_TOTAL_TRANSCRIPT_COUNT": self._data_frame.sum()}).transpose()

        self._samples = pandas.Series(
            [name.split('_')[0] for name in list(self._data_frame.index)],
            index=self._data_frame.index)

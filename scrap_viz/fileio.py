import os
import h5py
from scipy.sparse import csc_matrix
from scipy.io import mmread
import pandas
import numpy


def convert_mtx_to_csv(mtx_directory_path, csv_file_path):

    mtx_file_path = os.path.join(mtx_directory_path, "matrix.mtx")

    sparse_matrix = mmread(mtx_file_path)

    cell_barcodes_file = open(os.path.join(mtx_directory_path, "barcodes.tsv"))
    barcodes = cell_barcodes_file.read().splitlines()
    cell_barcodes_file.close()

    genes_file = open(os.path.join(mtx_directory_path, "genes.tsv"))
    genes = genes_file.read().splitlines()
    genes_file.close()

    short_genes_set = set([gene.split("\t")[1] for gene in genes])
    short_genes = []

    for gene in genes:

        _, gene_short = gene.split("\t")

        if gene_short in short_genes_set:

            index = 1

            while True:
                gene_short_plus_index = "%s_%i" % (gene_short, index)

                if gene_short_plus_index not in short_genes_set:
                    gene_short = gene_short_plus_index
                    short_genes_set.add(gene_short)
                    break
                else:
                    index += 1

        short_genes.append(gene_short)

    data_frame = pandas.SparseDataFrame(sparse_matrix)
    data_frame = data_frame.to_dense()
    data_frame = data_frame.fillna(0)
    data_frame = data_frame.astype(numpy.uint32)
    data_frame.index = short_genes
    data_frame.columns = barcodes

    write_pandas_csv(data_frame, csv_file_path)


def convert_h5_to_csv(h5_file_path, csv_file_path):

    hdf = h5py.File(h5_file_path)

    cellranger_version = 2

    if "matrix" in hdf:
        cellranger_version = 3

    matrix_name = None

    if cellranger_version == 2:
        for key, value in hdf.items():
            matrix_name = key
            break
    else:
        matrix_name = "matrix"

    sparse = csc_matrix(
        (hdf[matrix_name]["data"], hdf[matrix_name]["indices"], hdf[matrix_name]["indptr"]))

    df = pandas.SparseDataFrame(sparse)
    df = df.to_dense()
    df = df.fillna(0)
    df = df.astype(numpy.uint32)
    num_genes_present = df.shape[0]

    if cellranger_version == 2:
        num_genes = len(list(hdf[matrix_name]["genes"]))
        gene_names = [x.decode("UTF-8") for x in list(hdf[matrix_name]["gene_names"])]
    else:
        num_genes = len(list(hdf[matrix_name]["features"]))
        gene_names = [x.decode("UTF-8") for x in list(hdf[matrix_name]["features"]["name"])]
    barcodes = [x.decode("UTF-8") for x in list(hdf[matrix_name]["barcodes"])]

    if len(gene_names) != num_genes_present:
        print("Missing some genes: data frame is %i gene(s) but index has %i. \
               Filling rest with zeros." % (len(gene_names), len(df.index)))

    df.index = gene_names[0:num_genes_present]
    df.columns = barcodes
    if num_genes != num_genes_present:
        df = df.append(pandas.DataFrame(0, index=gene_names[num_genes_present:], columns=barcodes))
    write_pandas_csv(df, csv_file_path)


def write_pandas_csv(data_frame, file_path):
    pandas.DataFrame(data_frame)\
        .to_csv(file_path, sep=',', encoding='utf-8', chunksize=1000)

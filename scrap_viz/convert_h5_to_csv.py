import os
import h5py
from scipy.sparse import csc_matrix
import pandas
import numpy
from scrap_viz import Gene_Expression_Dataset


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

    df = pandas.SparseDataFrame(sparse, dtype=numpy.uint32)
    df = df.to_dense()
    df = df.fillna(0)
    df = df.astype(int)
    num_genes_present = df.shape[0]

    if cellranger_version == 2:
        num_genes = len(list(hdf[matrix_name]["genes"]))
        gene_names = [x.decode("UTF-8") for x in list(hdf[matrix_name]["gene_names"])]
    else:
        num_genes = len(list(hdf[matrix_name]["features"]))
        gene_names = [x.decode("UTF-8") for x in list(hdf[matrix_name]["features"]["name"])]
    barcodes = [x.decode("UTF-8") for x in list(hdf[matrix_name]["barcodes"])]
    df.index = gene_names[0:num_genes_present]
    df.columns = barcodes
    if num_genes != num_genes_present:
        df = df.append(pandas.DataFrame(0, index=gene_names[num_genes_present:], columns=barcodes))
    Gene_Expression_Dataset.write_pandas_csv(
        df, csv_file_path)

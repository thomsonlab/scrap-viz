import os
import h5py
from scipy.sparse import csc_matrix
import pandas
import numpy
from scRNA_seq import Gene_Expression_Dataset

workspace_path = os.getcwd()
transcript_file_path = os.path.join(workspace_path, "transcripts", "filtered_feature_bc_matrix.h5")

hdf = h5py.File(transcript_file_path)

sparse = csc_matrix(
    (hdf["matrix"]["data"], hdf["matrix"]["indices"], hdf["matrix"]["indptr"]))

df = pandas.SparseDataFrame(sparse, dtype=numpy.uint32)
df = df.to_dense()
df = df.fillna(0)
df = df.astype(int)
num_genes_present = df.shape[0]
num_genes = len(list(hdf["matrix"]["features"]))
gene_names = [x.decode("UTF-8") for x in list(hdf["matrix"]["features"]["name"])]
barcodes = [x.decode("UTF-8") for x in list(hdf["matrix"]["barcodes"])]
df.index = gene_names[0:num_genes_present]
df.columns = barcodes
if num_genes != num_genes_present:
    df = df.append(pandas.DataFrame(0, index=gene_names[num_genes_present:], columns=barcodes))
Gene_Expression_Dataset.write_pandas_csv(
    df, os.path.join(workspace_path, "transcripts", "gene_counts.csv"))

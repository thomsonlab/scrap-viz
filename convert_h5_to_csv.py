import os
import h5py
from scipy.sparse import csc_matrix
import pandas
import numpy
from scRNA_seq import Gene_Expression_Dataset

workspace_path = '/media/dibidave/mobinoodle/Virus_Transcriptomics/DRG'
transcript_file_path = os.path.join(workspace_path, "transcripts", "filtered_gene_bc_matrices_h5.h5")

hdf = h5py.File(transcript_file_path)

sparse = csc_matrix(
    (hdf["mm10"]["data"], hdf["mm10"]["indices"], hdf["mm10"]["indptr"]))

df = pandas.SparseDataFrame(sparse, dtype=numpy.uint32)
df = df.to_dense()
df = df.fillna(0)
df = df.astype(int)
df.index = [x.decode("UTF-8") for x in list(hdf["mm10"]["gene_names"])]
df.columns = [x.decode("UTF-8") for x in list(hdf["mm10"]["barcodes"])]
Gene_Expression_Dataset.write_pandas_csv(
    df, os.path.join(workspace_path, "transcripts", "gene_counts.csv"))

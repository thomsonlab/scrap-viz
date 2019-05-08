import os
import argparse
import shutil

from scrap_viz import fileio
from scrap_viz import Gene_Expression_Dataset


def get_arguments():

    argparser = argparse.ArgumentParser(
        description="Initialize a SCRAP workspace")
    argparser.add_argument("--workspace_path", "-w",
                           help="Path to the workspace", default=None)
    argparser.add_argument("--h5_file_path", "-f",
                           help="Path to the H5 file", default=None)
    argparser.add_argument("--csv_file_path", "-c",
                           help="Path to a CSV file", default=None)
    argparser.add_argument("--mtx_directory_path", "-m",
                           help="Path to an mtx directory", default=None)

    args = argparser.parse_args()

    return args


def initialize_dataset():

    args = get_arguments()

    if args.workspace_path is None:
        workspace_path = os.path.join(os.getcwd(), "workspace")
    else:
        workspace_path = args.workspace_path

    if os.path.isfile(workspace_path):
        raise ValueError("Cannot use file as a workspace")

    if not os.path.isdir(workspace_path):
        os.makedirs(workspace_path, exist_ok=True)

    if args.h5_file_path is None and \
            args.csv_file_path is None and \
            args.mtx_directory_path is None:
        raise ValueError("Must specify path to file with -f, -c, or -m")

    raw_gene_counts_csv_path = os.path.join(
        workspace_path, "raw_gene_counts.csv")

    if args.h5_file_path is not None:
        fileio.convert_h5_to_csv(args.h5_file_path, raw_gene_counts_csv_path)
    elif args.csv_file_path is not None:
        shutil.copyfile(args.csv_file_path, raw_gene_counts_csv_path)
    else:
        fileio.convert_mtx_to_csv(args.mtx_directory_path,
                                  raw_gene_counts_csv_path)

    seed_matrices_file_paths = [raw_gene_counts_csv_path]

    Gene_Expression_Dataset(
        workspace_path,
        seed_matrix_file_path=seed_matrices_file_paths
    )

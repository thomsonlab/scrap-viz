import os
import argparse

from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot


def get_arguments():

    argparser = argparse.ArgumentParser(
        description="Launch an scRNA-seq plotting server")
    argparser.add_argument("--pipeline_name", "-p",
                           help="Name of the pipeline to use", default=None)
    argparser.add_argument("--workspace_path", "-w",
                           help="Path to the workspace", default=None)

    args = argparser.parse_args()

    return args


def launch_server():

    args = get_arguments()

    if args.workspace_path is None:
        args.workspace_path = os.getcwd()

    pipeline_name = args.pipeline_name

    print("Loading dataset...")
    gene_expression_dataset = Gene_Expression_Dataset(
        args.workspace_path,
        name=pipeline_name
    )

    plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

    plot.start()

if __name__ == "__main__":
    launch_server()

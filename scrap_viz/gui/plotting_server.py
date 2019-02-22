#!/usr/bin/env python3

import os
import argparse

from scrap_viz import Gene_Expression_Dataset
from scrap_viz import Gene_Expression_Dataset_Plot
from scrap_viz import Gene_Metadata


def get_arguments():

    argparser = argparse.ArgumentParser(
        description="Launch an scRNA-seq plotting server")
    argparser.add_argument("--pipeline_name", "-p",
                           help="Name of the pipeline to use", default=None)
    argparser.add_argument("--workspace_path", "-w",
                           help="Path to the workspace", default=None)
    argparser.add_argument("--host_port", "-o",
                           help="What port to host on", default=8050, type=int)

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

    gene_metadata = Gene_Metadata(dir_path=args.workspace_path)

    print("Dataset loaded, creating plot...")

    plot = Gene_Expression_Dataset_Plot(gene_expression_dataset,
                                        gene_metadata,
                                        port=args.host_port)

    plot.start()


if __name__ == "__main__":
    launch_server()

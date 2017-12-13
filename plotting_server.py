import os
from scRNA_seq import Gene_Expression_Dataset
from scRNA_seq import Gene_Expression_Dataset_Plot


dataset_path = os.path.expanduser(
    os.path.join(
        "~", "Virus_Transcriptomics", "workspace"
    )
)

pipeline_name = "5_1000_RPM_SD"

print("Loading dataset...")
gene_expression_dataset = Gene_Expression_Dataset(
    dataset_path,
    name=pipeline_name
)

plot = Gene_Expression_Dataset_Plot(gene_expression_dataset)

plot.start()


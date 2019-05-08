# SCRAP-viz

## Installing

To install (or upgrade), in your conda or python environment:

```
pip install --upgrade git+https://github.com/ThomsonLab/SCRAP-viz.git
```

## Running

### Initializing a dataset
For each new dataset to analyze, you need to initialize a SCRAP-viz workspace for it. This workspace folder will automatically save labels you make and is where differential expression analyses will be exported to.

To initialize a dataset, from your command prompt, terminal, or conda prompt:
```
scrap-init-ws -w [path to the workspace you want to create] -f [path to your H5 file]
```

For example:
```
scrap-init-ws -w test_workspace -f /home/dibidave/Downloads/filtered_feature_bc_matrix.h5
```
Will create a folder "test_workspace" in the current folder I'm in on my command prompt, and will initialize it to be a SCRAP-viz workspace based on the gene count data in the given H5 file.

### Preprocessing a dataset
```
scrap-preprocess -w [path to the workspace] -p [what to call this preprocessing pipeline]
```
For example:
```
scrap-preprocess -w test_workspace -p my_preprocessing
```
Will process the gene counts in the "test_workspace" folder, and will name the preprocessing pipeline "my_preprocessing"

### Running the GUI
```
scrap-viz -w [path to the workspace] -p [name of the preprocessing pipeline to use]
```
For example:
```
scrap-viz -w test_workspace -p my_preprocessing
```
Will launch a GUI server from the data in the "my_preprocessing" preprocessing pipeline

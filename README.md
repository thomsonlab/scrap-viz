# SCRAP-viz

## Installing

To install, in your conda or python environment:
pip install git+https://github.com/ThomsonLab/SCRAP-viz.git

## Running

### Intializing a dataset
For each new dataset to analyze, you need to initialize a SCRAP-viz workspace for it. This workspace folder will automatically save labels you make and is where differential expression analyses will be exported to.

To initialize a dataset, from your command prompt, terminal, or conda prompt:
```
scrap-init-ws -w [path to the workspace you want to create] -f [path to your H5 file]
```

### Preprocessing a dataset
```
scrap-preprocess -w [path to the workspace] -p [what to call this preprocessing pipeline]
```

### Running the GUI
```
scrap-viz -w [path to the workspace] -p [name of the preprocessing pipeline to use]
```

from setuptools import setup

setup(
    name='scrap-viz',
    version='0.2.0',
    description='Single-cell RNA sequencing analysis tools and visualization',
    url='https://github.com/ThomsonLab/scrap-viz',
    author='David Brown',
    author_email='dibidave@gmail.com',
    license='MIT',
    packages=[
        "scrap_viz",
        "scrap_viz.gui"
    ],
    python_requires='~=3.6',
    entry_points={
        "console_scripts": [
            "scrap-viz=scrap_viz.gui.plotting_server:launch_server",
            "scrap-preprocess=scrap_viz.preprocessing:preprocess",
            "scrap-init-ws=scrap_viz.initialize_dataset:initialize_dataset"
        ]
    },
    install_requires=[
        "pandas",
        "numpy",
        "sklearn",
        "scipy",
        "statsmodels",
        "scvi",
        "torch",
        "dash>=0.37.0",
        "dash-core-components",
        "dash-html-components",
        "dash-renderer",
        "plotly",
        "flask"
    ]

)

from setuptools import setup

setup(
    name='scrap-viz',
    version='0.1.0',
    description='Single-cell RNA sequencing analysis tools and visualization',
    url='https://github.com/ThomsonLab/scrap-viz',
    author='David Brown',
    author_email='dibidave@gmail.com',
    license='MIT',
    packages=['scrap_viz'],
    python_requires='~=3.6',
    entry_points={
        'console_scripts': [
            'plotting_server=plotting_server:main',
        ],
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

from setuptools import setup

setup(
    name='scRNA-seq',
    version='0.1.0',
    description='Single-cell RNA sequencing analysis tools and visualization',
    url='https://github.com/GradinaruLab/scRNA-seq',
    author='David Brown',
    author_email='dibidave@gmail.com',
    license='MIT',
    packages=['scRNA_seq'],
    python_requires='~=3.6',
    entry_points={
        'console_scripts': [
            'plotting_server=plotting_server:main',
        ],
    },

)

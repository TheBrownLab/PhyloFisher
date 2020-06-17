from setuptools import setup, find_packages

setup(
    name='phylofisher',
    version='0.1.14',
    packages=find_packages(),
    scripts=['phylofisher/fisher.py',
             'phylofisher/config.py',
             'phylofisher/fishing_net.py',
             'phylofisher/forest.py',
             'phylofisher/forge.py',
             'phylofisher/informant.py',
             'phylofisher/lumberjack.py',
             'phylofisher/purge.py',
             'phylofisher/build_dataset.py',
             'phylofisher/install_deps.py',
             'phylofisher/single_gene_tree_constructor.py',
             'phylofisher/missing_data.py',
             'phylofisher/utilities/fast_site_removal.py',
             'phylofisher/utilities/mammal_modeler.py',
             'phylofisher/utilities/bipartition_examiner.py',
             'phylofisher/utilities/fast_tax_removal.py',
             'phylofisher/utilities/aa_comp_calculator.py',
             'phylofisher/utilities/SR4_class_recoder.py',
             'phylofisher/utilities/taxon_collapser.py',
             'phylofisher/utilities/genetic_code.py',
             'phylofisher/utilities/heterotachy.py',
             'phylofisher/utilities/random_sample_iteration.py'
             ],
    python_requires='>=3.7',
    install_requires=['biopython==1.76',
                      'pyqt5>=5',
                      'ete3==3.1.1',
                      'pandas==1.0.3',
                      'matplotlib==3.1.3',
                      'scipy'
                      ],
    url='https://github.com/TheBrownLab/PhyloFisher',
    license='MIT',
    author='David Zihala',
    author_email='zihaladavid@gmail.com',
    description='PhyloFisher is a software package for the creation, analysis, and visualization of phylogenomic '
                'datasets that consist of protein sequences from eukaryotic organisms.'
)

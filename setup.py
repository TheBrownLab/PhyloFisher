from setuptools import setup

setup(
    name='phylofisher',
    version='0.1',
    packages=['phylofisher'],
    scripts=['phylofisher/fisher.py',
             'phylofisher/config.py',
             'phylofisher/fishing_net.py',
             'phylofisher/forest.py',
             'phylofisher/forge.py',
             'phylofisher/informant.py',
             'phylofisher/lumberjack.py',
             'phylofisher/purge.py',
             'phylofisher/build_dataset.py',
             'phylofisher/single_gene_tree_constructor.py',
             'phylofisher/utilities/missing_data.py',
             'phylofisher/utilities/fast_site_removal.py',
             'phylofisher/utilities/bipartition_examiner.py',
             'phylofisher/utilities/fast_tax_removal.py',
             'phylofisher/utilities/heteroevolving_sites.py',
             'phylofisher/utilities/aa_comp_calculator.py',
             'phylofisher/utilities/SR4_class_recoder.py',
             'phylofisher/utilities/taxon_collapser.py',
             'phylofisher/utilities/genetic_code.py'
             ],

    url='',
    license='MIT',
    author='david',
    author_email='zihaladavid@gmail.com',
    description=''
)

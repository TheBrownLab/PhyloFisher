Installation guide
------------------

- install conda via Miniconda or Anaconda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/


- prepare conda virtual environment:

.. code-block:: bash

    cd PhyloFisher
    conda env create -f fisher_env.yml

- activate phylo_fisher environment:

.. code-block:: bash

    conda activate phylo_fisher

- intall phylofisher:

.. code-block:: bash

    pip install .

- after you finish using PhyloFisher use `conda deactivate` to deactivate phylo_fisher env

.. code-block:: bash

    conda deactivate
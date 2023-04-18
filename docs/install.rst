Install
=======

Conda
-----

.. code-block:: console

   $ conda install -c bioconda hmnfusion

pip
---

.. code-block:: console

   $ wget https://github.com/guillaume-gricourt/HmnFusion/releases/download/{{ version }}/pip.zip
   $ unzip pip.zip
   $ pip install hmnfusion-{{ version }}-py3-none-any.whl
   $ rm pip.zip hmnfusion-{{ version }}-py3-none-any.whl hmnfusion-{{ version }}.tar.gz

Docker
------
Pull hmnfusion

.. code-block:: console
    $ docker pull ghcr.io/guillaume-gricourt/hmnfusion:{{ version }}

Pull hmmnfusion-align to create BAM files

.. code-block:: console
    $ docker pull ggricourt/hmnfusion-align:{{ version }}

The reference files use to build BAM files could be cite with this DOI:
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6619597.svg
   :target: https://doi.org/10.5281/zenodo.6619597

.. note::
    The size of the image is near to 15Gb


Dependencies
------------

Genefuse and Lumpy are available in the ``docker`` image.
But they are not installed with ``conda`` and ``pip``.
You can install them in the current environment with this:

.. code-block:: console
    $ hmnfusion install-software

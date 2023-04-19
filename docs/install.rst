Install
=======

Conda
-----

.. code-block:: console

    $ conda install -c bioconda hmnfusion

Pip
---

.. code-block:: console
   :substitutions:

    $ wget https://github.com/guillaume-gricourt/HmnFusion/releases/download/|version|/pip.zip
    $ unzip pip.zip
    $ pip install hmnfusion-|version|-py3-none-any.whl
    $ rm pip.zip hmnfusion-|version|-py3-none-any.whl hmnfusion-|version|.tar.gz

Docker
------

Pull ``hmnfusion``. It contains Genefuse, Lumpy and HmnFusion.

.. code-block:: console
   :substitutions:

    $ docker pull ghcr.io/guillaume-gricourt/hmnfusion:|version|

Pull ``hmmnfusion-align`` to create BAM files

.. code-block:: console
   :substitutions:

    $ docker pull ggricourt/hmnfusion-align:|version|

The reference files used to build BAM files refer to the DOI 10.5281/zenodo.6619597

.. warning::
    The size of the image is almost 15Go


Dependencies
------------

Genefuse and Lumpy are available in the ``docker`` image.
But they are not installed with ``conda`` and ``pip``.
You can install them in the current environment with this:

.. code-block:: console

    $ hmnfusion install-software

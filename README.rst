======
CLSify
======

Classification of *Candidatus Liberibacter solanacearum* following IPPC (International Plant Protection Convention) standard `DP 21: Candidatus Liberibacter solanacearum <https://www.ippc.int/en/publications/84157>`_.

-----------
Quick Facts
-----------

- License: MIT
- Programming Language Python

------------------------------
Input / Output - What it Does!
------------------------------

This program takes as the input Sanger sequences from the 16S, 16S-23S, and 50S primers from the IPPC standard DP21.
It then aligns them to the GenBank reference sequences ``EU812559`` and ``EU834131`` (as specified in DP21).
Based on the alignments and the document DP21, sequence identity is computed and haplotyping is performed, yielding:

- sequence identity to ``EU822559`` for identifying the species *C. Liberibacter solanacearum*, and
- haplotyping of the read based on variation from the reference sequence.

Sample names can be inferred from the read names or from a separate mapping TSV file.

-----------
Quick Start
-----------

This is gonna be really quick!

Installation
============

We recommend using `Bioconda <https://bioconda.github.io>`_.

First `install Bioconda <https://bioconda.github.io/user/install.html#getting-started>`_.
Then:

.. code-block:: bash

    # conda install -y clsify

And -- tadaa -- you're ready to go!

Running
=======

You can have one FASTA (or FASTQ) file with all of your reads or one file for each.
If you have a single sequence per FASTA (or FASTQ) file then you can use the file name instead of the sequence name.

.. code-block:: bash

    # clsify -o result.tsv INPUT.fasta
    ## OR
    # clsify [--use-file-name] -o result.tsv INPUT1.fasta INPUT2.fasta [...]
    ## e.g.,
    # clsify [--use-file-name] -o result.tsv INPUT*.fasta

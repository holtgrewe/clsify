.. overview_cli:

======================
Command Line Interface
======================

You can alos run Haplotype-Lso locally from the command line.

.. note:: This section needs some extension.

.. code-block:: shell

    $ hlso cli \
        [--sample-name-from-file] \
        [--sample-regex REGEX] \
        [--output OUTPUT] \
        seq_file [seq_file ...]

This will read all sequence files ``seq_file`` (can be FASTA, FASTQ, AB1, SCF), perform conversion to FASTA (if needed) and then perform a haplotyping.
When provided, the result will be written to the XLSX file ``OUTPUT``.

You can override the regular expression to extract the sample name and region from the query name with ``--sample-regex``.

By default, the query sequence names are taken from their identifier.
As this is hard for the binary files AB1 and SCF, you can also configure Haplotype-Lso to use the file names (without extension) as the sample names.
This is the behaviour from the web frontend.

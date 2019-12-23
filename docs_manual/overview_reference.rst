.. _overview-reference:

======================
Reference Construction
======================

This section describes how the reference set (reference sequences and haplotyping position) can be generated.

.. note:: This section needs extension and more explanation.

-----
Input
-----

The overall input is a TSV file ``seeds_accessions`` that lists for each haplotype and region a GenBank accession with the prototype sequence.

.. code-block:: shell

    species              haplotype  region   accession   source

    Ca. L. solanacearum  A          16S      FJ498802.1  .
    Ca. L. solanacearum  B          16S      FJ939136.1  .
    Ca. L. solanacearum  C          16S      GU373048.1  .
    Ca. L. solanacearum  D          16S      HQ454302.1  .
    Ca. L. solanacearum  E          16S      KF737348.1  .

    Ca. L. solanacearum  A          16S-23S  FJ830690.1  .
    Ca. L. solanacearum  B          16S-23S  FJ830700.1  .
    Ca. L. solanacearum  C          16S-23S  JX280523.1  .
    Ca. L. solanacearum  D          16S-23S  JX308304.1  .
    Ca. L. solanacearum  E          16S-23S  KF737347.1  .

    Ca. L. solanacearum  A          50S      EU834131.1  .
    Ca. L. solanacearum  B          50S      FJ498807.1  .
    Ca. L. solanacearum  C          50S      GU373051.1  .
    Ca. L. solanacearum  D          50S      HQ454317.1  .
    Ca. L. solanacearum  E          50S      KY777461.1  .

-----------------------
Download Seed Sequences
-----------------------

.. code-block:: shell

    $ hlso cli ref_download \
        path/to/seeds_accession.tsv \
        path/to/seeds_paths.tsv

This will download sequences by accession, download them next to the ``seeds_accession.tsv`` file.
It will write the file ``seeds_paths.tsv`` with the names of the downloaded files:

.. code-block:: shell

    species              haplotype  region   accession   path
    Ca. L. solanacearum  A          16S      FJ498802.1  FJ498802.1.fasta
    Ca. L. solanacearum  B          16S      FJ939136.1  GU373048.1.fasta
    Ca. L. solanacearum  C          16S      GU373048.1  GU373048.1.fasta
    Ca. L. solanacearum  D          16S      HQ454302.1  HQ454302.1.fasta
    Ca. L. solanacearum  E          16S      KF737348.1  KF737348.1.fasta
    Ca. L. solanacearum  A          16S-23S  FJ830690.1  FJ830690.1.fasta
    Ca. L. solanacearum  B          16S-23S  FJ830700.1  FJ830700.1.fasta
    Ca. L. solanacearum  C          16S-23S  JX280523.1  JX280523.1.fasta
    Ca. L. solanacearum  D          16S-23S  JX308304.1  JX308304.1.fasta
    Ca. L. solanacearum  E          16S-23S  KF737347.1  KF737347.1.fasta
    Ca. L. solanacearum  A          50S      EU834131.1  EU834131.1.fasta
    Ca. L. solanacearum  B          50S      FJ498807.1  FJ498807.1.fasta
    Ca. L. solanacearum  C          50S      GU373051.1  GU373051.1.fasta
    Ca. L. solanacearum  D          50S      HQ454317.1  HQ454317.1.fasta
    Ca. L. solanacearum  E          50S      KY777461.1  KY777461.1.fasta

-----------------------------
Performing Seed BLAST Queries
-----------------------------

The next step is to perform a BLAST search via NCBI WWWBLAST to obtain sequences similar to the seeds.

.. code-block:: shell

    $ hlso ref_blast path/to/seeds_paths.tsv

For each seeed query ``accession.fasta``, a file ``accession.blast.xml`` will be generated with the BLAST results.

----------------------------
Consensus and Table Creation
----------------------------

Finally, compute consensus sequences and the haplotyping table.

.. code-block:: shell

    $ hls ref_consensus path/to/seeds_path.tsv \
        --output-table haplotype_table.txt

This will perform a consensus computation of the seeds, generate a haplotype-specific sequence for each region and each haplotype, and create a haplotyping table.

The file ``haplotype_table.txt`` can then be used for the haplotyping of sequences themselves.

.. code-block::

    # --------  -------  ------  ------  ------  ---------------  ------  ------  ------  ------  ------
    reference   region   pos     ref     alt     description      A       B       C       D       E
    # --------  -------  ------  ------  ------  ---------------  ------  ------  ------  ------  ------
    EU812559.1  16S      108     AT      A       n.109delT        AT      AT      AT      AT      A
    EU812559.1  16S      115     A       G       n.115A>G         A       A       A       A       G
    EU812559.1  16S      116     C       T       n.116C>T         C       C       C       T       C
    EU812559.1  16S      151     A       G       n.151A>G         A       A       A       A       G
    EU812559.1  16S      212     T       G       n.212T>G         T       G       T       T       T
    EU812559.1  16S      581     T       C       n.581T>C         T       C       T       T       T
    EU812559.1  16S      959     C       T       n.959C>T         C       C       C       C       T
    EU812559.1  16S      1039    A       G       n.1039A>G        A       A       G       G       A
    EU812559.1  16S      1039    AC      A       n.1040delC       AC      AC      AC      AC      A
    EU812559.1  16S      1073    G       A       n.1073G>A        G       G       G       A       G
    # --------  -------  ------  ------  ------  ---------------  ------  ------  ------  ------  ------

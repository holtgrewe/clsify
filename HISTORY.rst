=======
History
=======

------
v0.4.3
------

- Fixing issue with AB1/SCF support.

------
v0.4.2
------

- Fix for multi-region sequences.
- Fixing collision, variant must be identified by sequence, position, and reference.

------
v0.4.1
------

- Fixing bug with display of BLAST match.
- Adding link-out to NCBI WWWBLAST submission.
- Removing some cruft from repository.
- Fixing bug in phylogenetic coputation in case of BLAST all-to-all matches not a square number.

------
v0.4.0
------

- Compute region-wise phylogenetic tree.

------
v0.3.4
------

- Support for uploading to PyPi.

------
v0.3.3
------

- Fixing regular expression.

------
v0.3.2
------

- Changing file name pattern do dot-separated scheme.
- Greatly extending documentation.

------
v0.3.1
------

- Starting out with tutorial and manual.

------
v0.3.0
------

- Adding tutorial.
- Rebranding as "Haplotype-LSO" (``hlso``).
- Properly normalizing indels according to Tan et al. (2015).
- Adding support for haplotyping with indels.

------
v0.2.0
------

- Removing dependency on bcftools.
  Haplotyping is done from BLAST match now.
- Adding tests for ``blast`` module.
- Rewrite of the whole BLAST and haplotyping interface and architecture.

------
v0.1.1
------

- Zapping gremlins in haplotype table.
- Change formatting of README.

------
v0.1.0
------

Initial release.

- Everything is new!

# Test Data

The data accession numbers were obtained from the following sources:

- Nelson, W.R., Sengoda, V.G., Alfaro-Fernandez, A.O. et al.
  *A new haplotype of "Candidatus Liberibacter solanacearum" identified in the Mediterranean region.*
  Eur J Plant Pathol (2013) 135: 633.
  https://doi.org/10.1007/s10658-012-0121-3

The table fit the information regarding the data can be found in `samples.tsv`.

```bash
REF=../../hlso/data/ref_seqs.fasta 
for f in Nelson2012/*.fa; do
  blastn -db $REF -query $f -outfmt 15 \
  | jq \
  > ${f%.fa}.json
done
```

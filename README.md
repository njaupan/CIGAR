# CIGAR_SV
Detect TE insertion polymorphisms from long reads using CIGAR

## Map long reads to the reference using minimap2

```bash
minimap2 -cx map-ont -t $thread $reference $reads -o $paf
```
## Check CIGAR info using the script

```
python CIGAR_SV -i $paf
```

# CIGAR_SV
Detect TE insertion polymorphisms from long reads using CIGAR

## 1) Map long reads to the reference using minimap2

```bash
minimap2 -cx map-ont -t $thread $reference $reads -o $paf
```
## 2) Check CIGAR info using the script

```
python CIGAR_SV -i $paf

usage: python CIGAR_SV.py [-h] [-i PAF] [-m MIN__ALIGNED_LENGTH]
                            [-mq MIN__QUERY_LENGTH] [-p]

Statstics of minimap2 alignment results(.paf files)

optional arguments:
  -h, --help            show this help message and exit
  -i PAF, --paf PAF     paf file name
  -m MIN__ALIGNED_LENGTH, --min__aligned_length MIN__ALIGNED_LENGTH
  -mq MIN__QUERY_LENGTH, --min__query_length MIN__QUERY_LENGTH
  -p, --primary_aligntment
                        Exclude secondary and inversion alignments.


```
# Detected SV in csv format.


|Col |Type  |Description                               |
|---:|:----:|:-----------------------------------------|
|1   |string|query sequence name                       |
|2   |int   |query length                              |
|2   |int   |query start                               |
|3   |int   |query end                                 |
|4   |int   |direction                                 |
|5   |int   |reference sequence name                   |
|6   |int   |reference  length                         |
|7   |int   |reference  start                          |
|8   |int   |reference   end                           |
|9   |int   |match base                                |
|10  |int   |insertion length                          |
|11  |int   |insertion position                        |
|12  |int   |deletion length                           |
|13  |int   |deletion position                         |
|14  |int   |identity                                  |
|15  |int   |primary/secondary alignment               |
|16  |int   |query covrage                             |

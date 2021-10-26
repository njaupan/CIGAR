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
|5   |string|reference sequence name                   |
|6   |int   |reference  length                         |
|7   |int   |sv start  position on reference           |
|8   |int   |sv end  position on reference             |
|9   |int   |sv type (insertion or deletion) (ins/del) |
|10  |int   |mapping quality                           |
|11  |string|primary alignment                         |
|12  |int   |blast identity                            |
|13  |string|corresponding sv fasta sequence           |


# Contact
CIGAR_SV

Copyright © 2021 Panpan Zhang (njaupanpan@gmail.com)

Any question, concern, or bug report about the program should be posted as an Issue on GitHub. Before posting, please check previous issues (both Open and Closed) to see if your issue has been addressed already. Also, please follow these good GitHub practices.

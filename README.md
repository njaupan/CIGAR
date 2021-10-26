# CIGAR_SV
Detect TE insertion polymorphisms from long reads using CIGAR

## 1) Map long reads to the reference using minimap2

```bash
minimap2 -t 8 --cs -cx map-ont reference.fa reads.fasta/fastq > output.paf

```

## 2) Check CIGAR info using the script

```
python cigar_sv_fasta.py -i output.paf -r reference.fa -l reads.fasta/fastq ```


usage: python cigar_sv.py -i <input.paf> -r <reference.fa> -l <reads.fastq/fasta>  (option)

A tool to detect large structural variants (SVs) using long read sequencing

optional arguments:
  -h, --help            show this help message and exit
  -i <input.paf>        read alignment paf file with --cs -cx tag
  -r <reference.fa>     reference genome fasta (uncompressed or bgzipped)
  -l <reads.fastq/fasta>
                        query long read fastq/fasta file (uncompressed or bgzipped)
  -o PATH               output directory [./CIGAR_output]

filter options:
  -m INT                minimum length of a structural variant for detection [50]
  -ma INT               minimum length of the alignment block
  -mq INT               minimum length of the query read

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

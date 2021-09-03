#!/usr/bin/env python3

"""
   Metrics from paf file 
   Copyright 2020 Panpan Zhang (njaupanpan@gmail.com)

   This script is to get the statistics from paf file, such as Insertion, Deletion, Subsititution and Gap-compressed Identity by the same definition from minimap2 developer Heng Li's blog: http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
		       
	Tag  Type                       Description                      

	tp   A    Type of aln: P/primary, S/secondary and I,i/inversion 
	cm   i    Number of minimizers on the chain                     
    s1   i    Chaining score
    s2   i    Chaining score of the best secondary chain
    NM   i    Total number of mismatches and gaps in the alignment
    MD   Z    To generate the ref sequence in the alignment
    AS   i    DP alignment score
    ms   i    DP score of the max scoring segment in the alignment
    nn   i    Number of ambiguous bases in the alignment
    ts   A    Transcript strand (splice mode only)
    cg   Z    CIGAR string (only in PAF) M:MATCH; I:iNSERTION; D:DELETION
    cs   Z    Difference string

    Example commands:
    paf_identity.py -i input.paf (option -m / -mq/-p )

"""

import re
import sys
import argparse 

parser = argparse.ArgumentParser(description="Statstics of minimap2 alignment results(.paf files)")

parser.add_argument("-i","--paf",help="paf file name")
parser.add_argument('-m', '--min__aligned_length', type=int, default=None)
parser.add_argument('-mq', '--min__query_length', type=int, default=None)
parser.add_argument('-p', '--primary_aligntment',action='store_false',
                    help='Exclude secondary and inversion alignments.')
args=parser.parse_args()
# Define the search patterns:
pafFileNamePattern = re.compile('(.*)\.paf')
def statsFromPaf(pafFile):
    paf=open(pafFile)
    OutfilePrefix1 =pafFileNamePattern.match(pafFile) # define outfile 
    OutfilePrefix=OutfilePrefix1.group(1)      
    outfile=open(OutfilePrefix + '.paf.stats','w') 
    headers = ['queryID', 'qlen','qstart', 'qend' ,'direction', 'refID','rlen','rstart', 'rend','allmatch','ins','ins_start','delt','delt_start','blast_iden','PS','qcov']
    
    outfile.write('\t'.join(headers))
    outfile.write('\n')
    total=0
    n_primary=0
    #parts = pd.DataFrame([i.split('\t') for i in pafFile]), if install pandas module
    for line in paf:         
        parts = line.strip().split("\t")
        total=total+1
        #get tag "cg" for cigar
        cg=line.strip().split("cg:Z:")[1].split("\t")[0]
        PS=line.strip().split("tp:A:")[1].split("\t")[0]
        #DE=line.strip().split("de:f:")[1].split("\t")[0]
        #value=re.findall(r'(\d+)(\w)', cg)
        M= re.compile(r'(\d+)[M]').findall(cg)
        I= re.compile(r'(\d+)[I]').findall(cg)
        D= re.compile(r'(\d+)[D]').findall(cg)
        #L= re.compile(r'(\d+)[:]').findall(cs)

        match = sum(list(map(int, M)))
        ins = max(list(map(int, I)) or [0]) 

        l=list(map(int, I))
        from itertools import takewhile
        ins_a=sum(1 for _ in takewhile(lambda x: x< ins,l))
        il = l[0: ins_a]
        il_sum = sum(list(map(int, il)))
        il_pos = il_sum +int(parts[7])
        
        ld=list(map(int, D))
        delt = max(list(map(int, D)) or [0])
        delt_a=sum(1 for _ in takewhile(lambda y: y< delt,ld))
        dl = ld[0: delt_a]       
        dl_sum = sum(list(map(int, dl)))
        dl_pos = dl_sum +int(parts[7])
        
        # NM is edit distance: NM = INS + DEL + SUB,
        # but for minimap2 NM = #mismatches + #I + #D + #ambiguous_bases, NM does not count reference skip `N`.
        NM= int(parts[10])- int(parts[9])
        sub=NM-ins-delt
        ins_compress = cg.count("I")
        delt_compress = cg.count("D")
        blast_iden =  100.0 *int((parts[9]))/ int((parts[10]))
        gap_compress = ins_compress +delt_compress+sub
        gap_compress_iden = 100.0 *int((parts[9]))/(int(parts[10])-gap_compress)
        qcov=100.0 *int((parts[9]))/ int((parts[1]))
        resultss = {
		 "queryID": parts[0],
		 "qlen":  int(parts[1]),
		 "qstart": int(parts[2]),
		 "qend": int(parts[3]),
		 "direction": parts[4],
		 "refID": parts[5],
		 "rlen": int(parts[6]),
		 "rstart": int(parts[7]),
		 "rend": int(parts[8]),
		 "allmatch": int(parts[9]),
		 "qcov" : qcov,
                 "PS":PS,
                 #"divergence":DE,
		 "ins": ins,
		 "delt": delt,
		 "ins_start" : il_pos,
		 "delt_start" : dl_pos,
		 "sub": sub,
		 "gap_compress":gap_compress,
		 "blast_iden": blast_iden,
		 "gap_compress_iden": gap_compress_iden,
         }
        if args.min__aligned_length is not None and resultss['allmatch'] < args.min__aligned_length:
            n_primary=n_primary+1
            continue
        if args.min__query_length is not None and resultss['qlen'] < args.min__query_length:
            continue
        if not args.primary_aligntment:
            if resultss['primary/secondary'] is 'S':
               continue
        out_row = (str(resultss[x]) for x in headers)
        outfile.write('\t'.join(out_row))
        outfile.write('\n')
statsFromPaf(str(args.paf))

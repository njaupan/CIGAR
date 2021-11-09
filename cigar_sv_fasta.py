#!/usr/bin/env python3

"""
   Metrics from paf file 
   Copyright 2021 Panpan Zhang (njaupanpan@gmail.com)

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
    python cigar_sv.py -i input.paf (option -m / -mq/-p )

"""

__version__ = '1.0'

import os
import re
import sys
import argparse 
import subprocess
import multiprocessing

def statsFromPaf(inPaf,OutfilePrefix,output_path,min_align,min_query,min_sv):
    paf=open(inPaf) 
    outfile=open(output_path+OutfilePrefix + '.paf.CIGAR_SV.csv','w') 
    #headers = ['queryID', 'qlen','qstart', 'qend' ,'direction', 
    # 'refID','rlen','rstart', 'rend','allmatch','ins','ins_start','delt','delt_start',
    # 'blast_iden','PS','qcov']
    headers = ['queryID','qlen','qstart', 'qend','direction','refID','rlen','sv_start','sv_end','sv_type','sv_size','MQ','PS','blast_iden','sv_seq']
    outfile.write('\t'.join(headers))
    outfile.write('\n')

    total=0
    n_primary=0
    #parts = pd.DataFrame([i.split('\t') for i in pafFile]), if install pandas module
    for line in paf:         
        parts = line.strip().split("\t")
        total=total+1
        #get tag "cg" for cigar
        print( "Analysing read: " + str(parts[0]))
        PS=line.strip().split("tp:A:")[1].split("\t")[0]
        if PS == 'P':
            cs=line.strip().split("cs:Z:")[1].split("\t")[0].upper()
            Regex= re.compile(r'[A,C,G,T,N]+').findall(cs)
            #print(Regex)
            #print max sequence
            sv_seq=max(Regex, key=len, default='0')
            sv_size=len(sv_seq)

            cg=line.strip().split("cg:Z:")[1].split("\t")[0]
            #DE=line.strip().split("de:f:")[1].split("\t")[0]
            #find M,I,D
            M= re.compile(r'(\d+)[M]').findall(cg)
            I= re.compile(r'(\d+)[I]').findall(cg)
            D= re.compile(r'(\d+)[D]').findall(cg)

            #split by string and value 
            L=cg.replace("M", "M-").replace("I", "I-").replace("D", "D-").replace("H", "H-").replace("S", "S-")
            L1=L.split("-")
            
            ins = max(list(map(int, I)) or [0]) 
            delt = max(list(map(int, D)) or [0]) 

            if sv_size == ins:
                sv_type=str("ins")
                svID=str(sv_size)+"I"
                #find the position of the biggest sv size in the list
                sv_post=L1.index(svID)
                #sub-list until the index of the biggest sv
                svl = L1[0: sv_post]
                #extract the value of matches and deletions
                matchers=['M','D']
                match_d=[s for s in svl if any(xs in s for xs in matchers)]
                match_d=[x.replace('M','').replace('D','').replace('I','') for x in match_d]
                #sum up value of all matches and deletions
                svl_sum = sum(list(map(int,match_d)))
                #find sv position in chrom
                sv_start = int(svl_sum) +int(parts[7]) 
                ##print( "Large SV: " + str(parts[5]) + " " +str(dl_pos))
                sv_end =  sv_start

            elif sv_size == delt:
                sv_type=str("del")
                svID=str(sv_size)+"D"
                #find the position of the biggest sv size in the list
                sv_post=L1.index(svID)
                #sub-list until the index of the biggest sv
                svl = L1[0: sv_post]
                #extract the value of matches and deletions
                matchers=['M','D']
                match_d=[s for s in svl if any(xs in s for xs in matchers)]
                match_d=[x.replace('M','').replace('D','').replace('I','') for x in match_d]
                #sum up value of all matches and deletions
                svl_sum = sum(list(map(int,match_d)))
                #find sv position in chrom
                sv_start = int(svl_sum) +int(parts[7]) 
                ##print( "Large SV: " + str(parts[5]) + " " +str(dl_pos))
                sv_end =  sv_start + sv_size

            blast_iden =  100.0 *int((parts[9]))/ int((parts[10]))
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
            "match_block": int(parts[10]),
		    "qcov" : qcov,
            "PS":PS,
            "sv_type" : sv_type,
		    "sv_size": int(sv_size),
		    "sv_start" : sv_start,
            "sv_end" : sv_end,
            "sv_seq" : sv_seq,
            "MQ" : int(parts[11]),
		    "blast_iden": blast_iden,
            }


            if min_align is not None and resultss['match_block'] < min_align:
                n_primary=n_primary+1
                continue
            if min_query is not None and resultss['qlen'] < min_query:
                continue
            if sv_size is not None and resultss['sv_size'] < min_sv:
                continue

            out_row = (str(resultss[x]) for x in headers)
            outfile.write('\t'.join(out_row))
            outfile.write('\n')

def main():    

    description = "A tool to detect large structural variants (SVs) using long read sequencing"
    parser = argparse.ArgumentParser(description=description, usage="python cigar_sv.py -i <input.paf> -r <reference.fa> -l <reads.fastq/fasta>  (option)")

    parser.add_argument("-i",metavar="<input.paf>",help="read alignment paf file with --cs -cx tag")
    parser.add_argument("-r", metavar="<reference.fa>", type=str, default="", help="reference genome fasta (uncompressed or bgzipped)")
    parser.add_argument("-l",metavar="<reads.fastq/fasta>",help="query long read fastq/fasta file (uncompressed or bgzipped)")
    parser.add_argument("-o", metavar="PATH", type=str, default="CIGAR_output", help="output directory [./CIGAR_output]")

    filter_options = parser.add_argument_group("filter options") 

    filter_options.add_argument('-m', metavar="INT", type=int, default=50,  help="minimum length of a structural variant for detection [50]")
    filter_options.add_argument('-ma',  metavar="INT", type=int, default=None,help="minimum length of the alignment block")
    filter_options.add_argument('-mq', metavar="INT", type=int, default=None,help="minimum length of the query read")


    args=parser.parse_args()

    if not args.i or not args.l or not args.r:
        sys.exit("\n** The read alignment paf file, query long reads and reference genome files are required **")

    inPaf=args.i
    # Define the search patterns:
    pafFileNamePattern = re.compile('(.*)\.paf')
    OutfilePrefix1 =pafFileNamePattern.match(str(inPaf)) # define outfile 
    OutfilePrefix=OutfilePrefix1.group(1).split("/")[-1]   
    
    print(OutfilePrefix)

    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/" 

    reads=args.l
    reference=args.r
    min_sv=args.m
    min_align=args.ma
    min_query=args.mq


    print( "Detecting insertion and deletion positions from read alignment using CIGAR.")

    statsFromPaf(str(inPaf),OutfilePrefix,output_path,min_align,min_query,min_sv)

    print( "Insertion and deletion detection is done.")
    

    cmd="awk -v OFS='\\t' 'NR>1 {print \">\"$6\"_\"$8\"_\"$10\"_\"$11\"_\"$1\"_\"$2\"\\n\"$15}' " + output_path+OutfilePrefix + ".paf.CIGAR_SV.csv  > " +output_path+OutfilePrefix +".paf.sv.fa"
    #print(cmd)
    subprocess.call(cmd, shell=True)
    cmd="awk -v OFS='\\t' 'NR>1 {print $6,$8,$9,$1,$10,$11}' " + output_path+OutfilePrefix + ".paf.CIGAR_SV.csv|bedtools sort |bedtools merge -d 20 -c 4,5,6,4 -o count,distinct,min,distinct > " +output_path+OutfilePrefix +".paf.sv.bed"
    #print(cmd)
    subprocess.call(cmd, shell=True)

    print( "Extraction of insertion and deletion sequences is done.")
    print( "Goodbye ! Have a nice good ! ")
    
def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)

if __name__ == '__main__':
    main()
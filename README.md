# CIGAR_SV
Detect TE insertion polymorphisms from long reads using CIGAR

# 1st step, map long reads to the reference using minimap2

'''
minimap2 -cx map-ont -t $thread $reference $reads -o $paf

'''
# Check CIGAR info using the script

'''
python CIGAR_SV -i $paf
'''

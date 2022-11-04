#!/usr/bin/env python

try:
    import warnings
    import sys
    import os 
    from collections import defaultdict    
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get OTUs from mmseq2 result
LIBRARY = "38-7-R4A" #the name of the library that you are analyzing

# Step 1 Parse 'clusterRes_cluster.tsv' 
rep2all_reads = defaultdict(list) # rep => [all reads]
with open('%s/clusterRes_cluster.tsv' %LIBRARY, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        rep, read = tmp[0], tmp[1]
        rep2all_reads[rep].append(read)
lines.close()

# Step 2 Print OTU result
f = open('%s/OTU_result.txt' %LIBRARY , 'w')
for rep in rep2all_reads:
    line = rep + '\t' + ','.join(rep2all_reads[rep]) + '\n'
    f.write(line)
f.close()    
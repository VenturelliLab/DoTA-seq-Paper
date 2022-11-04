#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
 
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: To parse the blastn result for 16S rRNA gene classification    
    
def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                if (" " or "\t") in line: # Break at the first " " or "\t"
                    spliter = ""
                    for i in range(len(line)):
                        if line[i] == " " or line[i] == "\t":
                            spliter = line[i]
                            break 
                           
                    head = line.split(f'{spliter}', 1)[0]
                    seq_dict[head] = ""
                else:
                    head = line
                    seq_dict[head] = ""
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict

# Inputs:
LIBRARY = "38-7-R4A" #Name of library being analyzed
reads_fasta = '%s/clusterRes_rep_seq.fasta' %LIBRARY
blastn_result = '%s/%s.nochimera_rep.blastn_result.txt'%(LIBRARY, LIBRARY)
tax_ref = '/storage1/databases/QIIME2/qiime/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt'

# Step 1 Get query2length dict
query2length = {} # query => length of seq
seq = store_seq(reads_fasta)
for query_w_array in seq:
    query = query_w_array.replace('>','',1)
    length = len(seq[query_w_array])
    query2length[query] = length
    
# Step 2 Get query2hit
query2hit = {} # query => length of seq
with open(blastn_result, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        query, hit, iden, aligned_length = tmp[0], tmp[1], float(tmp[2]), float(tmp[3])
        if iden >= 95 and float(aligned_length / query2length[query]) > 0.9:
            query2hit[query] = hit
        
# Step 3 Get hit_map
hit_map = {} # hit => tax        
with open(tax_ref, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        hit, tax = tmp[0], tmp[1]
        hit_map[hit] = tax
        
# Step 4 Print result
f = open('%s/%s.nochimera_rep.blastn_result_parsed.txt'%(LIBRARY, LIBRARY), 'w')
for query in query2length:
    if query in query2hit:
        f.write(f"{query}\t{hit_map[query2hit[query]]}\n")
    else:
        f.write(f"{query}\tno significant hit\n")
f.close()        
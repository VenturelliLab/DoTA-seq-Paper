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
    
# Aim: To parse blastn result to get hit summary for each barcode (exclude calculating "no significant hit")
# The final result will be: barcode \t read_no (for hit1) \t read_no (for hit2) \t read_no (for hit3) ...

blastn_result_parsed = '38-7-R4A/38-7-R4A.nochimera_rep.blastn_result_parsed.txt'
otu_result_file = '38-7-R4A/OTU_result.txt'
output_file = '38-7-R4A/38-7-R4A.nochimera.hit_summary.wo_no_sig_hit_or_rare_OTU.txt'

# Step 1 Get rep_read2reads dict
rep_read2reads = {} # rep_read => [reads]
read2rep_read = {} # read => rep_read
with open(otu_result_file, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        rep_read, reads = tmp[0], tmp[1]
        reads_list = reads.split(',')
        rep_read2reads[rep_read] = reads_list
        
        for read in reads_list:
            read2rep_read[read] = rep_read
lines.close()  

# Step 2 Parse the input to get dicts
rep_read2hit = {} # rep_read => hit
with open(blastn_result_parsed, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        rep_read, hit = tmp[0], tmp[1]
        # Only use the tax string down to the family level    
        hit_at_family_level = hit.split(';D_5__')[0]
        rep_read2hit[rep_read] = hit_at_family_level
lines.close()

# Step 3 Get otu info for each barcode
barcode2reads_raw = defaultdict(list) # barcode => [reads]
barcode2otus = defaultdict(list) # barcode => [[reads, reads], [reads], [reads, reads, reads]]
for read in read2rep_read:
    barcode = read.rsplit('-')[0]
    barcode2reads_raw[barcode].append(read)
    
## Only keep barcodes that contain >= 30 reads
barcode2reads = {} 
for barcode in barcode2reads_raw:
    reads_no = len(barcode2reads_raw[barcode])
    if reads_no >= 30:
        barcode2reads[barcode] = barcode2reads_raw[barcode]
        
for barcode in barcode2reads:
    reads = barcode2reads[barcode]
    reads_set = set(reads)
    rep_reads_set = set()
    for read in reads_set:
        rep_reads_set.add(read2rep_read[read])
    for rep_read in rep_reads_set:
        all_otu_reads_set = set(rep_read2reads[rep_read])
        reads_set_for_this_barcode_and_otu = reads_set & all_otu_reads_set
        barcode2otus[barcode].append(list(reads_set_for_this_barcode_and_otu))

# Step 3 Get the result dict
barcode2info = defaultdict(dict) # barcode => hit => read_no
hits = set()
for barcode in barcode2reads:
    otus = barcode2otus[barcode]
    
    # Exclude OTUs with no significant hits
    otus_exclude_no_sig_hits = []
    reads_exclude_no_sig_hits = []
    hit2otu_list_1st_exclusion = defaultdict(list) # hit => [otus]
    for otu in otus:
        read_1st = otu[0]
        rep_read = read2rep_read[read_1st]
        hit = rep_read2hit[rep_read]
        if hit != 'no significant hit':
            hit2otu_list_1st_exclusion[hit].append(otu)
            otus_exclude_no_sig_hits.append(otu)
            for read in otu:
                reads_exclude_no_sig_hits.append(read)
            
    otus_exclude_no_sig_hits_no = len(otus_exclude_no_sig_hits)
    reads_exclude_no_sig_hits_no = len(reads_exclude_no_sig_hits)
    
    # Exclude rare OTUs (OTUs that contain <= 5% of total reads)
    rare_cutoff = 0.05
    otus_exclude_no_sig_hits_and_rare_otus = []
    reads_exclude_no_sig_hits_and_rare_otus = []
    hit2otu_list_2nd_exclusion = defaultdict(list) # hit => [otus]
    for hit in hit2otu_list_1st_exclusion:
        otu_list_1st_exclusion = hit2otu_list_1st_exclusion[hit]
        for otu in otu_list_1st_exclusion:
            read_perc = float(len(otu) / reads_exclude_no_sig_hits_no)
            if read_perc >= rare_cutoff:
                hit2otu_list_2nd_exclusion[hit].append(otu)
                otus_exclude_no_sig_hits_and_rare_otus.append(otu)
                for read in otu:
                    reads_exclude_no_sig_hits_and_rare_otus.append(read)
    
    for hit in hit2otu_list_2nd_exclusion:
        otu_list_2nd_exclusion = hit2otu_list_2nd_exclusion[hit]
        if otu_list_2nd_exclusion:
            hits.add(hit)
        read_no = 0
        for otu in otu_list_2nd_exclusion:
            for read in otu:
                read_no = read_no + 1
        barcode2info[barcode][hit] = read_no

hits_list = sorted(list(hits))
for barcode in barcode2info:
    for hit in hits_list:
        barcode2info[barcode][hit] = barcode2info[barcode].get(hit, 0)
    
# Step 4 Write down the result
f = open(output_file, 'w')
header = 'head' + '\t' + '\t'.join(hits_list) + '\n'
f.write(header)
for barcode in barcode2info:
    line_list = []
    line_list.append(barcode)
    for hit in hits_list:
        line_list.append(str(barcode2info[barcode][hit]))
    line =  '\t'.join(line_list) + '\n'
    f.write(line)   
f.close()    

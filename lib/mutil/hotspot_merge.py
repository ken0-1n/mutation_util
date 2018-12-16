import sys, os, subprocess, gzip
import pysam, glob

def merge_hotspot_list(in_hotspot_mutation, in_genomon_mutation, output_file, skip_h_header, skip_g_header):

    hotspot_hash = {}
    is_header = skip_h_header
    with open(in_hotspot_mutation, 'r') as hin:
        for line in hin:
            # skip header line
            if is_header:
                is_header = False
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chr = F[0]
            start = F[1]
            end = F[2]
            ref = F[3]
            alt = F[4]
            hotspot_hash[chr+"\t"+start+"\t"+end+"\t"+ref+"\t"+alt] = line

    hout = open(output_file, 'w')
    is_header = skip_g_header
    with open(in_genomon_mutation, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            # skip header line
            if is_header:
                is_header = False
                print >> hout, line +"\tscore(hotspot)"
                continue
            F = line.split('\t')
            chr = F[0]
            start = F[1]
            end = F[2]
            ref = F[3]
            alt = F[4]
            key = chr+"\t"+start+"\t"+end+"\t"+ref+"\t"+alt

            if key in hotspot_hash: 
                del hotspot_hash[key]
            print >> hout, line +"\t---"

    for v in hotspot_hash.values():
        print >> hout, v

    hout.close()


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

def merge_hotspot_vcf(in_hotspot_mutation, in_genomon_mutation, output_file):

    with open(output_file, 'w') as hout:
        hotspot_hash = {}
        with open(in_hotspot_mutation, 'r') as hin:
            for line in hin:
                line = line.rstrip('\n')
                # print metadata & header line
                if line.startswith("#"):
                    print >> hout, line
                    continue
                F = line.split('\t')
                chrom, pos, ids, ref, alt = F[0:5]
                hotspot_hash[chrom+"\t"+pos+"\t"+ref+"\t"+alt] = line
    
        with open(in_genomon_mutation, 'r') as hin:
            for line in hin:
                line = line.rstrip('\n')
                # skip header line
                if line.startswith("#"): continue
                F = line.split('\t')
                chrom, pos, ids, ref, alt, qual, filters, infos = F[0:8]
                key = chrom+"\t"+pos+"\t"+ref+"\t"+alt
                if key in hotspot_hash: 
                    hotspot_line =  hotspot_hash[key]
                    hotspot_infos = hotspot_line.split("\t")[7]
                    hotspot_infoF = hotspot_infos.split(";")
                    hotspot_lod_socre = ""
                    for hotspot_info_val in hotspot_infoF:
                        if hotspot_info_val.startswith("LS="):
                            hotspot_lod_socre=hotspot_info_val
                    del hotspot_hash[key]
                    print >> hout, "\t".join(F[0:7]) +"\t"+ infos+";"+hotspot_lod_socre +"\t"+ "\t".join(F[8:])
                else:
                    print >> hout, line

        for v in hotspot_hash.values():
            print >> hout, v


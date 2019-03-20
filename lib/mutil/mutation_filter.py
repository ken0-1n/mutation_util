import sys, os, subprocess, gzip
import vcf
from genomon_header_info import Genomon_header_info


###############################################
def filter_mutation_list(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, v_count, hotspot_database):

    # genomon header idx infomation object
    ghi = Genomon_header_info()
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            # print meta data
            if not line.startswith("#chr"):
                if line.startswith("#"):
                    print >> hout, line.rstrip('\n')
                    continue
            # get header line
            header = line
            ghi.set_header_information(header)
            print >> hout, header.rstrip('\n')
            break

    hotspot_list = []
    if os.path.exists(hotspot_database):
        with open(hotspot_database, "r") as HI:
            for line in HI:
                line = line.rstrip('\n')
                F = line.split("\t")
                chr = F[0]
                start = F[1]
                end = F[2]
                ref= F[3]
                alt = F[4]
                key = str(chr) +"\t"+ str(start) +"\t"+ str(end) +"\t"+ ref +"\t"+ alt
                hotspot_list.append(key)

    is_header = True
    with open(input_file, 'r') as hin:
        for line in hin:
            # skip meta data
            if not line.startswith("#chr"):
                if line.startswith("#"):
                    continue
            # skip header line
            if is_header:
                is_header = False
                continue

            line = line.rstrip('\n')
            F = line.split('\t')
            key = str(F[ghi.chr]) +"\t"+ str(F[ghi.start]) +"\t"+ str(F[ghi.end]) +"\t"+ F[ghi.ref].upper() +"\t"+ F[ghi.alt].upper()

            if (( ghi.fisher == -1 or float(F[ghi.fisher]) >= float(fishpval)) and \
                ( ghi.ebcall == -1 or float(F[ghi.ebcall]) >= float(ebpval))   and \
                ( ghi.realign == -1 or (F[ghi.realign] != "---" and float(F[ghi.realign]) >= float(realignpval))) and \
                ( ghi.tcount == -1 or  (F[ghi.tcount] != "---" and int(F[ghi.tcount]) >= int(tcount))) and \
                ( ghi.ncount == -1 or  (F[ghi.ncount] != "---" and int(F[ghi.ncount]) <= int(ncount))) and \
                ( ghi.post10q == -1 or float(F[ghi.post10q]) >= float(post10q)) and \
                ( ghi.r_post10q == -1 or (F[ghi.r_post10q] != "---" and float(F[ghi.r_post10q]) >= float(r_post10q))) and \
                ( ghi.v_count == -1 or   (F[ghi.v_count] != "---" and int(F[ghi.v_count]) >= int(v_count)))):
                print >> hout, line

            elif ( ghi.score_hotspot != -1 and (F[ghi.score_hotspot]) != "---") :
                print >> hout, line

            elif key in hotspot_list:
                print >> hout, line

    hout.close()

###############################################
def filter_mutation_vcf(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, sample1, sample2):

    vcf_reader = vcf.Reader(filename = input_file)
    vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)

    for record in vcf_reader:
        # TODO EBCall
        # FP:  Fisher's Test
        # FPR: Fisher's Test processed with Realignment
        # B10: 10% posterior quantile 
        # B1R: 10% posterior quantile Processed with Realignment
        # NAR: Number of allelic reads (Tumor)|(Sample1)
        # NAR: Number of allelic reads (Normal)|(Sample2)
        if (( "FP"  not in record.INFO or record.INFO["FP"]  >= float(fishpval))   and \
            ( "FPR" not in record.INFO or record.INFO["FPR"] >= float(realignpval)) and \
            ( "B10" not in record.INFO or record.INFO["B10"] >= float(post10q))    and \
            ( "B1R" not in record.INFO or record.INFO["B1R"] >= float(r_post10q))  and \
            ( sample1 == None or record.genotype(sample1)["NAR"] >= int(tcount))   and \
            ( sample2 == None or record.genotype(sample2)["NAR"] <= int(ncount))):
            vcf_writer.write_record(record)

        # hotspot(LOD of score) is exist
        elif "LS" in record.INFO: 
            vcf_writer.write_record(record)

    vcf_writer.close()


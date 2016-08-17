import sys, os, subprocess, gzip
import pysam, glob
from genomon_header_info import Genomon_header_info

def make_tabix_db(in_genomon_mutation, output_prefix):

    ghi = Genomon_header_info()
    hout = open(output_prefix + ".tmp.bed", 'w')
    
    is_header = True
    with open(in_genomon_mutation, 'r') as hin:
        for line in hin:
            # skip meta data
            if line.startswith("#"): continue
            # skip header line
            if is_header:
                ghi.set_header_information(line)
                is_header = False
                continue
    
            line = line.rstrip('\n')
            F = line.split('\t')
    
            # skip indel 
            if F[ghi.ref] == '-' or F[ghi.alt] == '-': continue
    
            start = int(F[ghi.start]) - 1
            print >> hout, F[ghi.chr] +'\t'+ str(start) +'\t'+ F[ghi.end] +'\t'+ F[ghi.ref] +'\t'+ F[ghi.alt] +'\t'+ F[ghi.tdepth] +'\t'+ F[ghi.tvariant] +'\t'+ F[ghi.ndepth] +'\t'+F[ghi.nvariant]+'\t'+ F[ghi.fisher] +'\t'+ F[ghi.ebcall]
    
    hout.close()
    
    hout = open(output_prefix + ".bed", 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_prefix + ".tmp.bed"], stdout = hout)
    hout.close()

    # compress and index
    subprocess.call(["bgzip", "-f", output_prefix + ".bed"])
    subprocess.call(["tabix", "-p", "bed", output_prefix + ".bed.gz"])

    os.unlink(output_prefix + ".tmp.bed")


def compare_hotspot_list(in_hotspot_mutation, in_genomon_mutation, output_dir, is_header):

    if not os.path.isdir(output_dir): os.mkdir(output_dir)

    base, ext = os.path.splitext( os.path.basename(in_genomon_mutation) )
    tabix_prefix = output_dir +"/"+ base

    make_tabix_db(in_genomon_mutation, tabix_prefix)

    base, ext = os.path.splitext( os.path.basename(in_hotspot_mutation) )
    output_prefix = output_dir +"/"+ base

    print_header_flag = True
    tb = pysam.TabixFile(tabix_prefix + ".bed.gz")
    hout = open(output_prefix + ".genomon_mutation.txt", 'w')
    with open(in_hotspot_mutation, 'r') as hin:
        for line in hin:
            if is_header and print_header_flag:
                line = line.rstrip('\n')
                print >> hout,line + "\tGenomonDepthTumor\tGenomonVariantTumor\tGenomonDepthNormal\tGenomonVariantNormal\tGenomonPval(Fisher)\tGenomonPval(EBCall)"
                print_header_flag = False
                continue

            line = line.rstrip('\n')
            F = line.split('\t')
            chr = F[0]
            start = F[1]
            end = F[2]
            ref = F[3]
            alt = F[4]

            result = ""
            try:
                records = tb.fetch(chr, (int(start) - 1), int(end))
                for record_line in records:
                    record = record_line.split('\t')
                    ref_tb = record[3]
                    alt_tb = record[4]

                    if ref == ref_tb and alt == alt_tb:
                        result = "\t".join(record[5:])
             
            except ValueError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]

            if result != "":
                print >> hout, line +"\t"+ result
            else:
                print >> hout, line +"\t---\t---\t---\t---\t---\t---"

    hout.close()


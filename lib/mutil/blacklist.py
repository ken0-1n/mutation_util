import sys
import os
import pysam
from genomon_header_info import Genomon_header_info
import glob

def filter(in_genomon_mutation_glob, blacklist_tabix, output_dir, min_candidate):

    tb = pysam.TabixFile(blacklist_tabix)
    ghi = Genomon_header_info()

    files = glob.glob(in_genomon_mutation_glob) 
    for in_genomon_mutation in files:

        base, ext = os.path.splitext( os.path.basename(in_genomon_mutation) )
        barcode = base.split('_')[0]
        dirname = os.path.dirname(in_genomon_mutation)
        disease = dirname.split('/')[-3]
        if not os.path.exists(output_dir) :os.mkdir(output_dir)
        if not os.path.exists(output_dir + '/' + disease) :os.mkdir(output_dir + '/' + disease)
        result_file = output_dir +'/'+disease+'/'+barcode + '.genomon_mutation.result.filt.blacklist_filtered.txt'
        black_file = output_dir +'/'+disease+'/'+barcode + '.genomon_mutation.result.filt.error_list.txt'
        
        hout = open(result_file, 'w')
        hout_black = open(black_file, 'w')
        is_header = True
        with open(in_genomon_mutation, 'r') as hin:
            for line in hin:

                print_flag = True
                # skip meta data
                if line.startswith("#"):
                    print >> hout, line.rstrip('\n')
                    continue
                if is_header:
                    header = line
                    ghi.set_header_information(header)
                    print >> hout, header.rstrip('\n')
                    is_header = False
                    continue
                try:
                    F = line.rstrip('\n').split('\t')
                    records = tb.fetch(F[ghi.chr], (int(F[ghi.start]) - 1), int(F[ghi.end]))

                    for record_line in records:
                        record = record_line.split('\t')
                        ref_tb = record[3]
                        alt_tb = record[4]
                        sample_num = record[5]

                        if F[ghi.ref] == ref_tb and F[ghi.alt] == alt_tb and int(sample_num) >= min_candidate:
                            print >> hout_black, line.rstrip("\n")
                            print_flag = False
    
                except ValueError:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    # print >> sys.stderr, fname

                if print_flag:
                    print >> hout, line.rstrip("\n")


import sys, os, subprocess, gzip
import pysam, glob
from subprocess import check_call
from genomon_header_info import Genomon_header_info
from maf_header_info import Maf_header_info


###############################################
def check_gene_region (input_file, output_file, target, hi):

    selected_tabix_file = ""
    if ( target == "exon"):
        selected_tabix_file = "resource/refExon.bed.gz"
    elif ( target == "coding"):
        selected_tabix_file = "resource/refCoding.bed.gz"
    elif ( target == "pancan"):
        selected_tabix_file = "resource/refPancan.bed.gz"

    if (selected_tabix_file != ""):
        tb = pysam.TabixFile(selected_tabix_file)
        hout = open(output_file, 'w')
        is_header = True
        with open(input_file, 'r') as hin:
            for line in hin:
                # skip meta data
                if line.startswith("#"):
                    continue
                if is_header:
                    header = line
                    hi.set_header_information(header)
                    print >> hout, header.rstrip('\n')
                    is_header = False
                    continue
                refFlag = False
                try:
                    F = line.rstrip('\n').split('\t')
                    chr = F[hi.chr]
                    if not chr.startswith("chr"):
                        chr = "chr" + chr
                    records = tb.fetch(chr, (int(F[hi.start]) - 1), int(F[hi.end]))
                    for record_line in records:
                        refFlag = True
                except ValueError:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                if refFlag:
                    print >> hout, line.rstrip('\n')
        hout.close()

        return_filter_input = output_file
    else:
        return_filter_input = input_file
    
    return return_filter_input

###############################################
def check_hotspot (input_file, output_file, target, hi):
    
    selected_tabix_file = ""
    if ( target == "hotspot"):
        selected_tabix_file = "resource/hg19_cosmic70.bed.gz"

    print hi.__module__

    if (selected_tabix_file != ""):
        tb = pysam.TabixFile(selected_tabix_file)
        hout = open(output_file, 'w')
        is_header = True
        with open(input_file, 'r') as hin:
            for line in hin:
                # skip meta data
                if line.startswith("#"):
                    continue
                if is_header:
                    header = line
                    hi.set_header_information(header)
                    print >> hout, header.rstrip('\n')
                    is_header = False
                    continue
                refFlag = False
                try:
                    F = line.rstrip('\n').split('\t')
   
                    alt1 = ""
                    alt2 = ""
                    if hi.__module__ == "mutil.maf_header_info":
                         alt1 = F[hi.t_allele1]
                         alt2 = F[hi.t_allele2]
                    else:
                         alt1 = F[hi.alt]
                         alt2 = F[hi.alt]

                    records = ""
                    if F[hi.ref] == '-' or alt2 == '-':
                        records = tb.fetch(F[hi.chr], (int(F[hi.start]) - 11), (int(F[hi.end]) + 10))
                    else:
                        records = tb.fetch(F[hi.chr], (int(F[hi.start]) - 1), int(F[hi.end]))

                    for record_line in records:
                        record = record_line.split('\t')
                        ref_tb = record[3]
                        alt_tb = record[4]

                        # ins
                        if F[hi.ref] == '-' and ref_tb == "-":
                            score1 = exact_alignment(alt2, alt_tb)
                            score2 = exact_alignment(alt_tb, alt2)
                            if (float(score1) / float(len(alt2))) >= 0.8 and (float(score2) / float(len(alt_tb))) >= 0.8: 
                                refFlag = True

                        # del
                        elif alt2 == '-' and alt_tb == "-":
                            score1 = exact_alignment(F[hi.ref], ref_tb)
                            score2 = exact_alignment(ref_tb, F[hi.ref])
                            if (float(score1) / float(len(F[hi.ref]))) >= 0.8 and (float(score2) / float(len(ref_tb))) >= 0.8: 
                                refFlag = True

                        # SNV
                        elif F[hi.ref] == ref_tb and (alt1 == alt_tb or alt2 == alt_tb):
                            refFlag = True
    
                except ValueError:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]

                if refFlag:
                    print >> hout, line.rstrip('\n')
        hout.close()

        return_filter_input = output_file
    else:
        return_filter_input = input_file
    
    return return_filter_input


###############################################
def lift_over (input_file, output_prefix, map_chain):

    # check header line and NCBI_build
    liftOverFlag = False
    is_header = True
    mhi = Maf_header_info()
    with open(input_file, 'r') as hin:
        for line in hin:
            if is_header:
                header = line
                mhi.set_header_information(header)
                is_header = False
                continue

            F = line.rstrip('\n').split('\t')
            NCBI_Build = F[mhi.ncbi_build]
            if (NCBI_Build == "37" or NCBI_Build == "hg19" or NCBI_Build == "GRCh37" or NCBI_Build == "GRCh37-lite"):
                liftOverFlag = False
            elif (NCBI_Build == "36" or NCBI_Build == "hg18"):
                liftOverFlag = True
            else:
               raise ValueError("An unexptected NCBI_Build code: " + NCBI_Build + " file: " + input_file)
            break

    return_liftover_file =""

    if (liftOverFlag):
        is_header = True
        hout = open(output_prefix + ".liftover_input.bed", 'w')
        with open(input_file, 'r') as hin:
            for line in hin:
                # skip header line
                if is_header:
                    is_header = False
                    continue

                F = line.rstrip('\n').split('\t')
                start = F[mhi.start]
                if F[mhi.variant_type] in ('SNP','DEL',):
                    start = (int(start) - 1)

                print >> hout, "chr" + F[mhi.chr] +'\t'+ str(start) +'\t'+ F[mhi.end] +'\t'+ F[mhi.chr] +','+ F[mhi.start] +','+ F[mhi.end]
        hout.close()

        subprocess.call(["liftOver", output_prefix +".liftover_input.bed", map_chain, output_prefix +".liftover_output.bed", output_prefix + ".liftover_unmap.bed"])

        lift_over_result_dict = {}
        with open(output_prefix +".liftover_output.bed", 'r') as hin:
            for line in hin:
                record = line.rstrip('\n').split('\t')
                lift_over_result_dict[record[3]] = record[0] +"\t"+ record[1] +"\t"+ record[2]

        is_header = True
        hout = open(output_prefix + ".liftover_input.txt", 'w')
        with open(input_file, 'r') as hin:
            for line in hin:
                # skip header line
                if is_header:
                    header = line
                    print >> hout, header.rstrip('\n')
                    is_header = False
                    continue

                F = line.rstrip('\n').split('\t')
                key =  F[mhi.chr] +","+ F[mhi.start] +","+ F[mhi.end]

                liftover_position = (lift_over_result_dict[key]).split('\t')
                chr = liftover_position[0].replace('chr','')
                start = liftover_position[1]
                end = liftover_position[2] 

                if F[mhi.variant_type] in ('SNP','DEL',):
                    start = int(start) + 1

                print >> hout, "\t".join(F[0:4]) +"\t"+ chr +"\t"+ str(start) +"\t"+ end +"\t"+ "\t".join(F[7:])
        hout.close()

        return_liftover_file = output_prefix +".liftover_input.txt"
    else:
        return_liftover_file = input_file

    return return_liftover_file


###############################################
def make_tabix_db(input_file, output_prefix):

    mhi = Maf_header_info()
    with open(input_file, 'r') as hin:
        for line in hin:
            header = line
            mhi.set_header_information(header)
            break
   
    # convert the firehose Maf format to bed format
    is_header = True
    hout = open(output_prefix + ".tmp.bed", 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            # skip header line
            if is_header:
                is_header = False
                continue

            line = line.rstrip('\n')
            F = line.split('\t')

            t_ref_count = F[mhi.t_ref_count] if mhi.t_ref_count != -1 else ""
            t_alt_count = F[mhi.t_alt_count] if mhi.t_alt_count != -1 else ""
          
            if F[mhi.variant_type] == 'DNP':
                start = int(F[mhi.start]) - 1
                end = int(F[mhi.end])
                allist1 = list(F[mhi.t_allele1])
                allist2 = list(F[mhi.t_allele2])
                for i, value in enumerate(list(F[mhi.ref])):
                    if value != allist1[i]:
                        print >> hout, F[mhi.chr] +'\t'+ str(start + i) +'\t'+ str(end + i) +'\t'+ value +'\t'+ allist1[i] +'\t'+ F[mhi.t_barcode] +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
                    if value != allist2[i]:
                        print >> hout, F[mhi.chr] +'\t'+ str(start + i) +'\t'+ str(end + i) +'\t'+ value +'\t'+ allist2[i] +'\t'+ F[mhi.t_barcode] +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
            else:
                start = F[mhi.start]
                if F[mhi.variant_type] in ('SNP','DEL'):
                    start = int(F[mhi.start]) - 1

                if F[mhi.ref] != F[mhi.t_allele1]:
                    print >> hout, F[mhi.chr] +'\t'+ str(start) +'\t'+ F[mhi.end] +'\t'+ F[mhi.ref] +'\t'+ F[mhi.t_allele1] +'\t'+ F[mhi.t_barcode] +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
                if F[mhi.t_allele1] != F[mhi.t_allele2] and F[mhi.ref] != F[mhi.t_allele2]:
                    print >> hout, F[mhi.chr] +'\t'+ str(start) +'\t'+ F[mhi.end] +'\t'+ F[mhi.ref] +'\t'+ F[mhi.t_allele2] +'\t'+ F[mhi.t_barcode] +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
          
    hout.close()
   
    hout = open(output_prefix + ".bed", 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_prefix + ".tmp.bed"], stdout = hout)
    hout.close()
    
    # compress and index
    subprocess.call(["bgzip", "-f", output_prefix + ".bed"])
    subprocess.call(["tabix", "-p", "bed", output_prefix + ".bed.gz"])
    
    # remove intermediate file
    # subprocess.call(["rm", "-rf", output_prefix + ".tmp.bed"])


###############################################
def exact_alignment(seq1, seq2):

    A = [[0]* len(seq2) for i in range(len(seq1))]

    best_score = 0
    opt_coord = (0, 0)

    for i in range(len(seq1)):
        for j in range(len(seq2)):

            if i == 0:
                A[i][j] = 1 if seq1[i] == seq2[j] else 0
            elif j == 0:
                A[i][j] = 1 if seq1[i] == seq2[j] else 0
            else:
                A[i][j] = A[i - 1][j - 1] + 1 if A[i - 1][j - 1] > 0 and seq1[i] == seq2[j] else 0

            if A[i][j] > best_score:
                best_score = A[i][j]
                opt_coord = (i, j)

    return best_score


###############################################
def compare_list(in_genomon_mutation, output_dir, data_file, map_chain, ebpval, fishpval, realignpval, tcount, ncount, gene_ref, post10q, r_post10q, v_count, print_graph, pancan, hotspot):

    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    
    base, ext = os.path.splitext( os.path.basename(data_file) )
    output_prefix = output_dir +"/"+ base
    
    data_file = lift_over(data_file, output_prefix, map_chain)

    mhi = Maf_header_info()
    data_file = check_gene_region(data_file, output_prefix + ".ref_input.txt", gene_ref, mhi)
    mhi = Maf_header_info()
    data_file = check_gene_region(data_file, output_prefix + ".pancan_input.txt", pancan, mhi)
    mhi = Maf_header_info()
    data_file = check_hotspot(data_file, output_prefix + ".hotspot_input.txt", hotspot, mhi)

    make_tabix_db(data_file, output_prefix)
    tb = pysam.TabixFile(output_prefix +".bed.gz")
    db_file = output_prefix +".tmp.bed"

    base, ext = os.path.splitext( os.path.basename(in_genomon_mutation) )
    genomon_output_prefix = output_dir +"/"+ base
    ghi = Genomon_header_info()
    in_genomon_mutation= check_gene_region(in_genomon_mutation, genomon_output_prefix + ".ref_input.txt", gene_ref, ghi)
    ghi = Genomon_header_info()
    in_genomon_mutation= check_gene_region(in_genomon_mutation, genomon_output_prefix + ".pancan_input.txt", pancan, ghi)
    ghi = Genomon_header_info()
    in_genomon_mutation= check_hotspot(in_genomon_mutation, genomon_output_prefix + ".hotspot_input.txt", hotspot, ghi)
    result_genomon = output_dir +"/"+ base +"_firehose.txt"
    result_firehose = output_prefix +"_firehose_only.maf.txt"
 
    position_db_dict = {}
    hout = open(result_genomon, 'w')

    additional_header = "Tumor_Sample_Barcode"
    ghi = Genomon_header_info()
    # output meta and header line
    with open(in_genomon_mutation, 'r') as hin:
        for line in hin:
            # print meta data
            if line.startswith("#"):
                print >> hout, line.rstrip('\n') + "\t" + additional_header
                continue
            # get header line
            header = line
            ghi.set_header_information(header)
            print >> hout, header.rstrip('\n') + "\t" + additional_header
            break

    is_header = True
    with open(in_genomon_mutation, 'r') as hin:
        for line in hin:
            # skip meta data
            if line.startswith("#"):
                continue
            # skip header line
            if is_header:
                is_header = False
                continue
    
            line = line.rstrip('\n')
            F = line.split('\t')
            
            result = ""
            record_key = ""
            try:
                records = ""
                if F[ghi.ref] == '-' or F[ghi.alt] == '-':
                    records = tb.fetch(F[ghi.chr], (int(F[ghi.start]) - 11), (int(F[ghi.end]) + 10))
                else:
                    records = tb.fetch(F[ghi.chr], (int(F[ghi.start]) - 1), int(F[ghi.end]))

                for record_line in records:
                    record = record_line.split('\t')
                    tmp_record_key = record[0] +"\t"+ record[1] + "\t" + record[2]
                    tmp_result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
                    ref_tb = record[3]
                    alt_tb = record[4]
                    type_tb = record[17]
    
                    # ins
                    if F[ghi.ref] == '-' and type_tb == "INS":
                        score1 = exact_alignment(F[ghi.alt], alt_tb)
                        score2 = exact_alignment(alt_tb, F[ghi.alt])
                        if (float(score1) / float(len(F[ghi.alt]))) >= 0.8 and (float(score2) / float(len(alt_tb))) >= 0.8: 
                            record_key = tmp_record_key
                            result = tmp_result

                    # del
                    elif F[ghi.alt] == '-' and type_tb == "DEL":
                        score1 = exact_alignment(F[ghi.ref], ref_tb)
                        score2 = exact_alignment(ref_tb, F[ghi.ref])
                        if (float(score1) / float(len(F[ghi.ref]))) >= 0.8 and (float(score2) / float(len(ref_tb))) >= 0.8: 
                            record_key = tmp_record_key
                            result = tmp_result

                    # SNV
                    elif F[ghi.ref] == ref_tb and F[ghi.alt] == alt_tb:
                        record_key = tmp_record_key
                        result = tmp_result
    
            except ValueError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
  
#           if (( gene_ref == "" or F[ghi.gene] in gene_ref.split(',')) and \
#           if (( func_ref == "" or F[ghi.func] in func_ref.split(',')) and \
            if (( ghi.fisher == -1 or float(F[ghi.fisher]) >= float(fishpval)) and \
                ( ghi.ebcall == -1 or float(F[ghi.ebcall]) >= float(ebpval))   and \
                ( ghi.realign == -1 or F[ghi.realign] == "---" or float(F[ghi.realign]) >= float(realignpval)) and \
                ( ghi.tcount == -1 or F[ghi.tcount] == "---" or int(F[ghi.tcount]) >= int(tcount)) and \
                ( ghi.ncount == -1 or F[ghi.ncount] == "---" or int(F[ghi.ncount]) <= int(ncount)) and \
                ( ghi.post10q == -1 or float(F[ghi.post10q]) >= float(post10q)) and \
                ( ghi.r_post10q == -1 or F[ghi.r_post10q] == "---" or float(F[ghi.r_post10q]) >= float(r_post10q)) and \
                ( ghi.v_count == -1 or F[ghi.v_count] == "---" or int(F[ghi.v_count]) >= int(v_count))):

                 # if record_key != "":
                 position_db_dict[record_key] = "Exists\t"+ F[ghi.fisher] +"\t"+ F[ghi.ebcall] +"\t"+ F[ghi.realign] +"\t"+ F[ghi.tcount] +"\t"+ F[ghi.ncount]
                 if result == "": result = "\t\t\t"
                 print >> hout, line + result

            else:
                 # if record_key != "":
                 position_db_dict[record_key] = "Filtered\t"+ F[ghi.fisher] +"\t"+ F[ghi.ebcall] +"\t"+ F[ghi.realign] +"\t"+ F[ghi.tcount] +"\t"+ F[ghi.ncount]
    
    hout.close()

    hout = open(result_firehose, 'w')

    additional_header = "merge_status\tP-value(fisher)\tP-value(EBCall)\tP-value(fhsher_realignment)\tvariantPairNum_tumor\tvariantPairNum_normal"
    with open(data_file, 'r') as hin:
        header = hin.readline()
        header = header.rstrip()
        print >> hout, header + "\t" + additional_header

    pre_key = ""
    with open(db_file, 'r') as hin:

        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            key = F[0] +"\t"+ F[1] + "\t" + F[2]

            if key == pre_key: continue

#           gene_name = F[8]
#           if (gene_ref == "" or gene_name in gene_ref.split(',')):

            if key in position_db_dict:
                print >> hout, "\t".join(F[8:]) +"\t"+ position_db_dict[key]
            else:
                print >> hout, "\t".join(F[8:])

            pre_key = key

    hout.close()

    # make Venn diagram
    ####################################################
    genomon_total = 0
    firehose_total = 0
    count = 1
    hout_snv = open( output_dir + '/print_R_SNV_tmp.txt', 'w')
    hout_indel = open( output_dir + '/print_R_INDEL_tmp.txt', 'w')

    mhi = Maf_header_info()
    is_header = True
    with open(result_firehose, 'r') as hin:
        for line in hin:
            # skip header line
            if is_header:
                # get header line
                header = line
                mhi.set_header_information(header)
                is_header = False
                continue

            F = line.rstrip('\n').split('\t')
            genomon_count = str(count) if len(F) > mhi.merge_status and F[mhi.merge_status] == "Exists" else "0"
    
            if F[mhi.variant_type] == "INS" or F[mhi.variant_type] == "DEL":
                print >> hout_indel, str(count) +"\t"+ genomon_count
            elif F[mhi.variant_type] == "SNP" or F[mhi.variant_type] == "DNP":
                print >> hout_snv, str(count) +"\t"+ genomon_count
            count += 1
            firehose_total += 1
            if genomon_count != 0: genomon_total += 1

    ghi = Genomon_header_info()
    is_header = True
    with open(result_genomon, 'r') as hin:
        for line in hin:
            # skip header line
            if is_header:
                # get header line
                header = line
                ghi.set_header_information(header)
                is_header = False
                continue

            F = line.rstrip('\n').split('\t')
            genomon_count = str(count)
        
            if F[ghi.tumor_barcode] == "":
                if F[ghi.ref] == "-" or F[ghi.alt] == "-":
                    print >> hout_indel, "0"+"\t"+ genomon_count
                else:
                    print >> hout_snv, "0"+"\t"+ genomon_count
                count += 1
                genomon_total += 1

    hout_snv.close()
    hout_indel.close()

    if print_graph:
        if os.path.getsize(output_dir+'/print_R_SNV_tmp.txt'):
            os.system('R --vanilla --slave  --args ' +output_dir+'/print_R_SNV_tmp.txt ' +output_prefix+'.snv.tiff '+base+' '+str(fishpval)+' '+str(ebpval)+' '+str(realignpval)+' '+str(tcount)+' '+str(ncount)+' '+str(genomon_total)+' '+str(firehose_total)+' < script/venn_mutation.R')
        if os.path.getsize(output_dir+'/print_R_INDEL_tmp.txt'):
            os.system('R --vanilla --slave  --args ' +output_dir+'/print_R_INDEL_tmp.txt ' +output_prefix+'.indel.tiff '+base+' '+str(fishpval)+' '+str(ebpval)+' '+str(realignpval)+' '+str(tcount)+' '+str(ncount)+' '+str(genomon_total)+' '+str(firehose_total)+' < script/venn_mutation.R')

###############################################
def filt_mutation_result(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, v_count, hotspot_database):

    # genomon header idx infomation object
    ghi = Genomon_header_info()
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            # print meta data
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
def compare_all(in_genomon_mutation_glob, output_dir, data_file_dir, map_chain, ebpval, fishpval, realignpval, tcount, ncount, gene_ref, post10q, r_post10q, v_count, pancan, hotspot):

    if not os.path.exists(output_dir): os.mkdir(output_dir)
    
    files = glob.glob(in_genomon_mutation_glob) 
    for in_genomon_mutation in files:
        base, ext = os.path.splitext( os.path.basename(in_genomon_mutation) )
        barcode = base.split('_')[0]
        data_file = data_file_dir + '/' + barcode + '.maf.txt'
        if os.path.exists(data_file):
            print in_genomon_mutation
            print data_file
            compare_list(in_genomon_mutation, output_dir, data_file, map_chain, ebpval, fishpval, realignpval, tcount, ncount, gene_ref, post10q, r_post10q, v_count, False, pancan, hotspot) 

    # make Venn diagram
    ####################################################
    count = 1
    genomon_total = 0
    firehose_total = 0
    basename = os.path.basename(output_dir)
    hout_snv = open( output_dir + '/print_R_SNV_tmp.txt', 'w')
    hout_indel = open( output_dir + '/print_R_INDEL_tmp.txt', 'w')
    hout_genomon = open( output_dir +"/"+ basename +".mutation.result_firehose.txt", 'w')
    hout_firehose = open( output_dir + '/'+basename+".firehose_result.txt", 'w')
    
    files = glob.glob(output_dir +"/TCGA*_firehose_only.maf.txt") 
    merge_header = True
    for out_firehose in files:
        print out_firehose
        mhi = Maf_header_info()

        fname = os.path.basename(out_firehose)
        sample = fname.split(".")[0]

        is_header = True
        with open(out_firehose, 'r') as hin:
            for line in hin:
                # skip header line
                if is_header:
                    # get header line
                    header = line
                    mhi.set_header_information(header)
                    if merge_header:
                        print >> hout_firehose, "sample\t" + header.rstrip('\n')
                        merge_header = False
                    is_header = False
                    continue

                F = line.rstrip('\n').split('\t')
                genomon_count = str(count) if len(F) > mhi.merge_status and F[mhi.merge_status] == "Exists" else "0"
    
                if F[mhi.variant_type] == "INS" or F[mhi.variant_type] == "DEL":
                    print >> hout_indel, str(count) +"\t"+ genomon_count
                elif F[mhi.variant_type] == "SNP" or F[mhi.variant_type] == "DNP":
                    print >> hout_snv, str(count) +"\t"+ genomon_count
                print >> hout_firehose, sample +"\t"+line.rstrip('\n')
                count += 1
                firehose_total += 1
                if genomon_count != 0: genomon_total += 1

    files = glob.glob(output_dir +"/TCGA*mutation.result_firehose.txt") 
    merge_header = True
    for out_genomon_mutation in files:
        print out_genomon_mutation

        fname = os.path.basename(out_genomon_mutation)
        sample = fname.split(".")[0]

        ghi = Genomon_header_info()
        is_header = True
        with open(out_genomon_mutation, 'r') as hin:
            for line in hin:
                # skip header line
                if is_header:
                    # get header line
                    header = line
                    ghi.set_header_information(header)
                    if merge_header:
                        print >> hout_genomon, "sample\t" + header.rstrip('\n')
                        merge_header = False
                    is_header = False
                    continue

                F = line.rstrip('\n').split('\t')
                genomon_count = str(count)
        
                if F[ghi.tumor_barcode] == "":
                    if F[ghi.ref] == "-" or F[ghi.alt] == "-":
                        print >> hout_indel, "0"+"\t"+ genomon_count
                    else:
                        print >> hout_snv, "0"+"\t"+ genomon_count
                    count += 1
                    genomon_total += 1
                print >> hout_genomon, sample +"\t"+ line.rstrip('\n')

    hout_snv.close()
    hout_indel.close()

    base = os.path.basename(output_dir)
    if os.path.getsize(output_dir+'/print_R_SNV_tmp.txt'):
        os.system('R --vanilla --slave  --args ' +output_dir+'/print_R_SNV_tmp.txt ' +output_dir+'/'+base+'.snv.tiff '+base+' '+str(fishpval)+' '+str(ebpval)+' '+str(realignpval)+' '+str(tcount)+' '+str(ncount)+' '+str(genomon_total)+' '+str(firehose_total)+' < script/venn_mutation.R')
    if os.path.getsize(output_dir+'/print_R_INDEL_tmp.txt'):
        os.system('R --vanilla --slave  --args ' +output_dir+'/print_R_INDEL_tmp.txt ' +output_dir+'/'+base+'.indel.tiff '+base+' '+str(fishpval)+' '+str(ebpval)+' '+str(realignpval)+' '+str(tcount)+' '+str(ncount)+' '+str(genomon_total)+' '+str(firehose_total)+' < script/venn_mutation.R')



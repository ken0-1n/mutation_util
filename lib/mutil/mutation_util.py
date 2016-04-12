import sys, os, subprocess, gzip
import pysam
from genomon_header_info import Genomon_header_info


###############################################
def lift_over (input_file, output_prefix, map_chain):

    count = 1
    liftOverFlag = False

    # check header line and NCBI_build
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if (count == 2):
                NCBI_Build = F[3]
                if (NCBI_Build == "37" or NCBI_Build == "hg19" or NCBI_Build == "GRCh37" or NCBI_Build == "GRCh37-lite"):
                    liftOverFlag = False
                    break
                elif (NCBI_Build == "36" or NCBI_Build == "hg18"):
                    liftOverFlag = True
                    break
                else:
                   raise ValueError("An unexptected NCBI_Build code: " + NCBI_Build + " file: " + input_file)
            count += 1

    return_input_file =""

    if (liftOverFlag):
        hout = open(output_prefix + ".liftover_input.bed", 'w')
        with open(input_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                Chromosome = F[4]
                if Chromosome == "Chromosome": continue 
                Start_position = F[5]
                End_position = F[6]
                Variant_Type = F[9]
                
                start = Start_position
                if Variant_Type in ('SNP','DEL',):
                    start = int(Start_position) - 1

                print >> hout, "chr" + Chromosome +'\t'+ str(start) +'\t'+ End_position +'\t'+ Chromosome +','+ Start_position +','+ End_position
        hout.close()

        subprocess.call(["liftOver", output_prefix +".liftover_input.bed", map_chain, output_prefix +".liftover_output.bed", output_prefix + ".liftover_unmap.bed"])
        lift_over_result_dict = {}
        with open(output_prefix +".liftover_output.bed", 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                lift_over_result_dict[F[3]] = F[0] +"\t"+ F[1] +"\t"+ F[2]

        hout = open(output_prefix + ".new_input.txt", 'w')
        with open(input_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                Chromosome = F[4]
                if Chromosome == "Chromosome": continue 
                Start_position = F[5]
                End_position = F[6]
                Variant_Type = F[9]

                key =  F[4] +","+ F[5] +","+ F[6]

                liftF = (lift_over_result_dict[key]).split('\t')
                newChr = liftF[0].replace('chr','')
                newStart = liftF[1]
                newEnd = liftF[2] 

                if Variant_Type in ('SNP','DEL',):
                    newStart = int(newStart) + 1

                print >> hout, "\t".join(F[0:4]) +"\t"+ newChr +"\t"+ str(newStart) +"\t"+ newEnd +"\t"+ "\t".join(F[7:])
        hout.close()

        return_input_file = output_prefix +".new_input.txt"
    else:
        return_input_file = input_file

    return return_input_file


###############################################
def make_tabix_db(input_file, output_prefix):

    t_ref_count_idx = -1
    t_alt_count_idx = -1

    # check header line and NCBI_build
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            for i, v in enumerate(F):
                if v == "t_ref_count":
                    t_ref_count_idx = i
                elif v == "t_alt_count":
                    t_alt_count_idx = i
            break

    # convert the firehose Maf format to bed format
    hout = open(output_prefix + ".tmp.bed", 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            Chromosome = F[4]
            if Chromosome == "Chromosome": continue 

            Start_position = F[5]
            End_position = F[6]
            Variant_Type = F[9]
            Reference_Allele = F[10]
            Tumor_Seq_Allele1 = F[11]
            Tumor_Seq_Allele2 = F[12]
            Tumor_Sample_Barcode = F[15]
            t_ref_count = ""
            if t_ref_count_idx != -1:
                t_ref_count = F[t_ref_count_idx]
            t_alt_count = ""
            if t_alt_count_idx != -1:
                t_alt_count = F[t_alt_count_idx]
          
            if Variant_Type == 'DNP':
                start = int(Start_position) - 1
                end = int(End_position)
                allist1 = list(Tumor_Seq_Allele1)
                allist2 = list(Tumor_Seq_Allele2)
                for i, value in enumerate(list(Reference_Allele)):
                    if value != allist1[i]:
                        print >> hout, Chromosome +'\t'+ str(start + i) +'\t'+ str(end + i) +'\t'+ value +'\t'+ allist1[i] +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
                    if value != allist2[i]:
                        print >> hout, Chromosome +'\t'+ str(start + i) +'\t'+ str(end + i) +'\t'+ value +'\t'+ allist2[i] +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
            else:
                start = Start_position
                if Variant_Type in ('SNP','DEL'):
                    start = int(Start_position) - 1

                if Reference_Allele != Tumor_Seq_Allele1:
                    print >> hout, Chromosome +'\t'+ str(start) +'\t'+ End_position +'\t'+ Reference_Allele +'\t'+ Tumor_Seq_Allele1 +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
                if Reference_Allele != Tumor_Seq_Allele2:
                    print >> hout, Chromosome +'\t'+ str(start) +'\t'+ End_position +'\t'+ Reference_Allele +'\t'+ Tumor_Seq_Allele2 +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +"\t"+ line
          
        hout.close()
   
    
    hout = open(output_prefix + ".bed", 'w')
    subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_prefix + ".tmp.bed"], stdout = hout)
    hout.close()
    
    # compress and index
    subprocess.call(["bgzip", "-f", output_prefix + ".bed"])
    subprocess.call(["tabix", "-p", "bed", output_prefix + ".bed.gz"])
    
    # remove intermediate file
    # subprocess.call(["rm", "-rf", output_prefix + ".tmp.bed"])
    # subprocess.call(["rm", "-rf", output_prefix + ".tmp2.bed"])
    # subprocess.call(["rm", "-rf", output_prefix + ".tmp_chr.bed"])
    # subprocess.call(["rm", "-rf", output_prefix + ".unmapped.bed"])


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
def compare_list(in_genomon_mutation, output_dir, data_file, map_chain, ebpval, fishpval, realignpval, tcount, ncount, func_ref, gene_ref, post10q, r_post10q, v_count):

    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    
    base, ext = os.path.splitext( os.path.basename(data_file) )
    output_prefix = output_dir +"/"+ base
    
    new_data_file = lift_over(data_file, output_prefix, map_chain)
   
    make_tabix_db(new_data_file, output_prefix)
    tb = pysam.TabixFile(output_prefix +".bed.gz")
    db_file = output_prefix +".tmp.bed"
    
    base, ext = os.path.splitext( os.path.basename(in_genomon_mutation) )
    result_genomon = output_dir +"/"+ base +"_firehose.txt"
    result_firehose = output_prefix +"firehose_only.maf.txt"
 
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
                if F[ghi.ref_idx] == '-' or F[ghi.alt_idx] == '-':
                    records = tb.fetch(F[ghi.chr_idx], (int(F[ghi.start_idx]) - 11), (int(F[ghi.end_idx]) + 10))
                else:
                    records = tb.fetch(F[ghi.chr_idx], (int(F[ghi.start_idx]) - 1), int(F[ghi.end_idx]))

                for record_line in records:
                    record = record_line.split('\t')
                    tmp_record_key = record[0] +"\t"+ record[1] + "\t" + record[2]
                    tmp_result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
                    ref_tb = record[3]
                    alt_tb = record[4]
                    type_tb = record[17]
    
                    # ins
                    if F[ghi.ref_idx] == '-' and type_tb == "INS":
                        score1 = exact_alignment(F[ghi.alt_idx], alt_tb)
                        score2 = exact_alignment(alt_tb, F[ghi.alt_idx])
                        if (float(score1) / float(len(F[ghi.alt_idx]))) >= 0.8 and (float(score2) / float(len(alt_tb))) >= 0.8: 
                            record_key = tmp_record_key
                            result = tmp_result

                    # del
                    elif F[ghi.alt_idx] == '-' and type_tb == "DEL":
                        score1 = exact_alignment(F[ghi.ref_idx], ref_tb)
                        score2 = exact_alignment(ref_tb, F[ghi.ref_idx])
                        if (float(score1) / float(len(F[ghi.ref_idx]))) >= 0.8 and (float(score2) / float(len(ref_tb))) >= 0.8: 
                            record_key = tmp_record_key
                            result = tmp_result

                    # SNV
                    elif F[ghi.ref_idx] == ref_tb and F[ghi.alt_idx] == alt_tb:
                        record_key = tmp_record_key
                        result = tmp_result
    
            except ValueError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
  
            if (( gene_ref == "" or F[ghi.gene_idx] in gene_ref.split(',')) and \
                ( func_ref == "" or F[ghi.func_idx] in func_ref.split(',')) and \
                ( ghi.fisher_idx == -1 or float(F[ghi.fisher_idx]) >= float(fishpval)) and \
                ( ghi.ebcall_idx == -1 or float(F[ghi.ebcall_idx]) >= float(ebpval))   and \
                ( ghi.realign_idx == -1 or F[ghi.realign_idx] == "---" or float(F[ghi.realign_idx]) >= float(realignpval)) and \
                ( ghi.tcount_idx == -1 or F[ghi.tcount_idx] == "---" or int(F[ghi.tcount_idx]) >= int(tcount)) and \
                ( ghi.ncount_idx == -1 or F[ghi.ncount_idx] == "---" or int(F[ghi.ncount_idx]) <= int(ncount)) and \
                ( ghi.post10q_idx == -1 or float(F[ghi.post10q_idx]) >= float(post10q)) and \
                ( ghi.r_post10q_idx == -1 or F[ghi.r_post10q_idx] == "---" or float(F[ghi.r_post10q_idx]) >= float(r_post10q)) and \
                ( ghi.v_count_idx == -1 or F[ghi.v_count_idx] == "---" or int(F[ghi.v_count_idx]) >= int(v_count))):

                 # if record_key != "":
                 position_db_dict[record_key] = "Exists\t"+ F[ghi.fisher_idx] +"\t"+ F[ghi.ebcall_idx] +"\t"+ F[ghi.realign_idx] +"\t"+ F[ghi.tcount_idx] +"\t"+ F[ghi.ncount_idx]
                 if result == "": result = "\t\t\t"
                 print >> hout, line + result

            else:
                 # if record_key != "":
                 position_db_dict[record_key] = "Filtered\t"+ F[ghi.fisher_idx] +"\t"+ F[ghi.ebcall_idx] +"\t"+ F[ghi.realign_idx] +"\t"+ F[ghi.tcount_idx] +"\t"+ F[ghi.ncount_idx]
    
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

            gene_name = F[8]
            if (gene_ref == "" or gene_name in gene_ref.split(',')):

                if key in position_db_dict:
                    print >> hout, "\t".join(F[8:]) +"\t"+ position_db_dict[key]
                else:
                    print >> hout, "\t".join(F[8:])

                pre_key = key

    hout.close()


###############################################
def filt_mutation_result(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, v_count):

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

            if (( ghi.fisher_idx == -1 or float(F[ghi.fisher_idx]) >= float(fishpval)) and \
                ( ghi.ebcall_idx == -1 or float(F[ghi.ebcall_idx]) >= float(ebpval))   and \
                ( ghi.realign_idx == -1 or F[ghi.realign_idx] == "---" or float(F[ghi.realign_idx]) >= float(realignpval)) and \
                ( ghi.tcount_idx == -1 or F[ghi.tcount_idx] == "---" or int(F[ghi.tcount_idx]) >= int(tcount)) and \
                ( ghi.ncount_idx == -1 or F[ghi.ncount_idx] == "---" or int(F[ghi.ncount_idx]) <= int(ncount)) and \
                ( ghi.post10q_idx == -1 or float(F[ghi.post10q_idx]) >= float(post10q)) and \
                ( ghi.r_post10q_idx == -1 or F[ghi.r_post10q_idx] == "---" or float(F[ghi.r_post10q_idx]) >= float(r_post10q)) and \
                ( ghi.v_count_idx == -1 or F[ghi.v_count_idx] == "---" or int(F[ghi.v_count_idx]) >= int(v_count))):
                    
                print >> hout, line
    hout.close()


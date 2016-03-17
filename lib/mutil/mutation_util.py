import sys, os, subprocess, gzip
import pysam, logging, shutil


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
def compare_list(input_file, output_dir, data_file, map_chain, ebpval, fishpval, realignpval, tcount, ncount, func_ref, gene_ref):

    # if not os.path.exists(input_file):
    #     raise ValueError("file not exists: " + input_file)

    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    
    base, ext = os.path.splitext( os.path.basename(data_file) )
    output_prefix = output_dir +"/"+ base
    
    new_data_file = lift_over(data_file, output_prefix, map_chain)
   
    make_tabix_db(new_data_file, output_prefix)
    tb = pysam.TabixFile(output_prefix +".bed.gz")
    db_file = output_prefix +".tmp.bed"
    
    base, ext = os.path.splitext( os.path.basename(input_file) )
    result_file = output_dir +"/"+ base +"_firehose.txt"
    result_file3 = output_dir +"/"+ base +"_firehose_filtered.txt"
    result_file2 = output_prefix +"firehose_only.maf.txt"
    result_file4 = output_prefix +"firehose_only.maf_filtered.txt"
 
    fisher_idx = -1
    ebcall_idx = -1
    realign_idx = -1
    tcount_idx = -1
    ncount_idx = -1
    funcgene_idx = -1
    gene_idx = -1
    # check header line and NCBI_build
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            for i, v in enumerate(F):
                if v == "P-value(fisher)":
                    fisher_idx = i
                elif v == "P-value(EBCall)":
                    ebcall_idx = i
                elif v == "P-value(fhsher_realignment)":
                    realign_idx = i
                elif v == "variantPairNum_tumor":
                    tcount_idx = i
                elif v == "variantPairNum_normal":
                    ncount_idx = i
                elif v == "Func.refGene":
                    funcgene_idx = i
                elif v == "Gene.refGene":
                    gene_idx = i
            break

    position_db_dict = {}
    position_db_dict_filtered = {}
    hout = open(result_file, 'w')
    hout_filt = open(result_file3, 'w')
    with open(input_file, 'r') as hin:
        hin.readline()
        for line in hin:
    
            line = line.rstrip()
            itemlist = line.split('\t')
            
            chr = itemlist[0]
            start = (int(itemlist[1]) - 1)
            end = int(itemlist[2])
            ref = itemlist[3]
            alt = itemlist[4]
            
            result = "\t\t\t"
            # insertion and deletion
            if ref == '-' or alt == '-':
                start = start - 10 
                end = end + 10

            record_key = ""
            try:
                records = tb.fetch(chr, start, end)
                for record_line in records:
                    record = record_line.split('\t')
                    ref_db = record[3]
                    alt_db = record[4]
                    type_db = record[17]
    
                    # ins
                    if ref == '-' and type_db == "INS":
                        score1 = exact_alignment(alt, alt_db)
                        score2 = exact_alignment(alt_db, alt)
                        if (float(score1) / float(len(alt))) >= 0.8 and (float(score2) / float(len(alt_db))) >= 0.8: 
                            record_key = record[0] +"\t"+ record[1] + "\t" + record[2]
                            result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
                    # del
                    elif alt == '-' and type_db == "DEL":
                        score1 = exact_alignment(ref, ref_db)
                        score2 = exact_alignment(ref_db, ref)
                        if (float(score1) / float(len(ref))) >= 0.8 and (float(score2) / float(len(ref_db))) >= 0.8: 
                            record_key = record[0] +"\t"+ record[1] + "\t" + record[2]
                            result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]

                    # SNV
                    elif ref == ref_db and alt == alt_db:
                        record_key = record[0] +"\t"+ record[1] + "\t" + record[2]
                        result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
    
            except ValueError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                # logging.error( ("{0}:{1}:{2} {3}:{4}-{5}".format( exc_type, fname, exc_tb.tb_lineno, chr, start ,end) ) )
  
            if (func_ref != "" and itemlist[funcgene_idx] not in func_ref.split(',')) or \
               (gene_ref != "" and itemlist[gene_idx] not in gene_ref.split(',')) or \
               (fishpval > 0.0 and float(itemlist[fisher_idx]) < fishpval) or \
               (ebpval > 0.0 and float(itemlist[ebcall_idx]) < ebpval) or \
               (itemlist[realign_idx] != "---" and realignpval > 0.0 and float(itemlist[realign_idx]) < realignpval) or \
               (itemlist[tcount_idx] != "---" and tcount > -1 and int(itemlist[tcount_idx]) < tcount) or \
               (itemlist[ncount_idx] != "---" and ncount > -1 and int(itemlist[ncount_idx]) > ncount): \
               
                 # if record_key != "":
                 position_db_dict_filtered[record_key] = itemlist[fisher_idx] +"\t"+ itemlist[ebcall_idx] +"\t"+ itemlist[realign_idx] +"\t"+ itemlist[tcount_idx] +"\t"+ itemlist[ncount_idx]
                 print >> hout_filt, line + result
            else:
                 # if record_key != "":
                 position_db_dict[record_key] = itemlist[fisher_idx] +"\t"+ itemlist[ebcall_idx] +"\t"+ itemlist[realign_idx] +"\t"+ itemlist[tcount_idx] +"\t"+ itemlist[ncount_idx]
                 print >> hout, line + result
    
    hout.close()
    hout_filt.close()

    pre_key = ""
    hout = open(result_file2, 'w')
    with open(db_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            key = F[0] +"\t"+ F[1] + "\t" + F[2]

            if key == pre_key: continue

            newline = "\t".join(F[8:])
            if len(F) == 63:
                newline = newline +"\t\t\t"

            if position_db_dict.has_key(key):
                print >> hout, newline +"\t"+ position_db_dict[key]
            else:
                print >> hout, newline

            pre_key = key

    pre_key = ""
    hout = open(result_file4, 'w')
    with open(db_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            key = F[0] +"\t"+ F[1] + "\t" + F[2]

            if key == pre_key: continue

            newline = "\t".join(F[8:])
            if len(F) == 63:
                newline = newline +"\t\t\t"

            if position_db_dict_filtered.has_key(key):
                print >> hout, newline +"\t"+ position_db_dict_filtered[key]

            pre_key = key






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
                if (NCBI_Build == "37" or NCBI_Build == "hg19"):
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
            F = line.rstrip('\n').split('\t')
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
          
            start = Start_position
            if Variant_Type in ('SNP','DEL',):
                start = int(Start_position) - 1
    
            if Reference_Allele != Tumor_Seq_Allele1:
                print >> hout, Chromosome +'\t'+ str(start) +'\t'+ End_position +'\t'+ Reference_Allele +'\t'+ Tumor_Seq_Allele1 +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +'\t'
    
            if Reference_Allele != Tumor_Seq_Allele2:
                print >> hout, Chromosome +'\t'+ str(start) +'\t'+ End_position +'\t'+ Reference_Allele +'\t'+ Tumor_Seq_Allele2 +'\t'+ Tumor_Sample_Barcode +"\t"+ t_ref_count +'\t'+ t_alt_count +'\t'
           
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
def compare_list(input_file, output_dir, data_file, map_chain, ebpval, fishpval, realignpval, tcount, ncount):

    # if not os.path.exists(input_file):
    #     raise ValueError("file not exists: " + input_file)
    
    base, ext = os.path.splitext( os.path.basename(data_file) )
    output_prefix = output_dir +"/"+ base
    
    new_data_file = lift_over(data_file, output_prefix, map_chain)
   
    make_tabix_db(new_data_file, output_prefix)
    tb = pysam.TabixFile(output_prefix +".bed.gz")
    
    base, ext = os.path.splitext( os.path.basename(input_file) )
    result_file = output_dir +"/"+ base +"_firehose.txt"
    hout = open(result_file, 'w')
   
    fisher_idx = -1
    ebcall_idx = -1
    realign_idx = -1
    tcount_idx = -1
    ncount_idx = -1
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
            break

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
            
            if (fishpval > 0.0 and float(itemlist[fisher_idx]) < fishpval): continue
            if (ebpval > 0.0 and float(itemlist[ebcall_idx]) < ebpval): continue
            if (realignpval > 0.0 and float(itemlist[realign_idx]) < realignpval): continue 
            if (tcount > -1 and float(itemlist[tcount_idx]) < tcount): continue
            if (ncount > -1 and float(itemlist[ncount_idx]) > ncount): continue
    
            result = "\t\t\t"
            # insertion and deletion
            if ref != '-' or alt != '-':
                start = start - 10 
                end = end + 10

            # SNV
            try:
                records = tb.fetch(chr, start, end)
                for record_line in records:
                    record = record_line.split('\t')
                    ref_db = record[3]
                    alt_db = record[4]
    
                    if ref == '-':
                        score1 = exact_alignment(alt, alt_db)
                        score2 = exact_alignment(alt_db, alt)
                        # print "score:" +str(score1) + "/ "+ str(len(alt))
                        # print "score:" +str(score2) + "/ "+ str(len(alt_db))
                        # print float(score1) / float(len(alt))
                        # print float(score2) / float(len(alt_db))
                        if (float(score1) / float(len(alt))) >= 0.8 and (float(score2) / float(len(alt_db))) >= 0.8: 
                            result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
                    elif alt == '-':
                        score1 = exact_alignment(ref, ref_db)
                        score2 = exact_alignment(ref_db, ref)
                        if (float(score1) / float(len(ref))) >= 0.8 and (float(score2) / float(len(ref_db))) >= 0.8: 
                            result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
                    elif ref == ref_db and alt == alt_db:
                        result = "\t"+ record[5] +"\t"+ record[6] +"\t"+ record[7]
    
            except ValueError:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                logging.error( ("{0}:{1}:{2} {3}:{4}-{5}".format( exc_type, fname, exc_tb.tb_lineno, chr, start ,end) ) )
                print >> hout, line + result
                continue
    
            print >> hout, line + result
    hout.close()



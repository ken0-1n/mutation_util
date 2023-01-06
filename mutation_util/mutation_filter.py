import sys, os, subprocess, gzip
import vcf


###############################################
def filter_mutation_list(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, v_count, hotspot_database, ghi, flag_mis_base_0, fishpval_base_0, realignpval_base_0, ncount_other, ncount_depth):

    # genomon header idx infomation object
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
            print(header.rstrip('\n'),file=hout)
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

            if ( ghi.score_hotspot != None and (F[ghi.score_hotspot]) != "---") :
                print(line, file=hout)

            elif key in hotspot_list:
                print(line, file=hout)

            else:
                # single mode 
                if ghi.fisher == None:
                    if (( ghi.ebcall == None or float(F[ghi.ebcall]) >= ebpval)   and \
                        ( float(F[ghi.post10q]) >= post10q) and \
                        ( ghi.r_post10q != None and F[ghi.r_post10q] != "---" and float(F[ghi.r_post10q]) >= r_post10q) and \
                        ( ghi.v_count != None and F[ghi.v_count] != "---" and int(F[ghi.v_count]) >= v_count)):
                        print(line, file=hout)

                # pair mode
                else:
                    if (( ghi.ebcall == None or float(F[ghi.ebcall]) >= ebpval)   and \
                        ( float(F[ghi.fisher]) >= fishpval or \
                            ( int(F[ghi.nvariant]) == 0 and \
                              flag_mis_base_0 == True and \
                              float(F[ghi.fisher]) >= fishpval_base_0) \
                        ) and \
                        ( ghi.realign != None and F[ghi.realign] != "---" and \
                            ( float(F[ghi.realign]) >= realignpval or 
                                ( int(F[ghi.ncount]) == 0 and \
                                  flag_mis_base_0 == True and \
                                  float(F[ghi.realign]) >= realignpval_base_0)) \
                        ) and \
                        ( ghi.tcount != None and (F[ghi.tcount] != "---" and int(F[ghi.tcount]) >= tcount)) and \
                        ( ghi.ncount != None and (F[ghi.ncount] != "---" and int(F[ghi.ncount]) <= ncount)) and \
                        ( ghi.ncount_other != None and (F[ghi.ncount_other] != "---" and int(F[ghi.ncount_other]) <= ncount_other)) and \
                        ( ghi.ncount_ref != None and (F[ghi.ncount_ref] != "---" and int(F[ghi.ncount_ref]) + int(F[ghi.ncount]) >= ncount_depth))):
                        
                        print(line, file=hout)


    hout.close()

###############################################
def filter_mutation_vcf(input_file, output_file, ebpval, fishpval, realignpval, tcount, ncount, post10q, r_post10q, v_count, sample1, sample2, ghi, flag_mis_base_0, fishpval_base_0, realignpval_base_0, ncount_other, ncount_depth):

    vcf_reader = vcf.Reader(filename = input_file)
    vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)

    for record in vcf_reader:

        # hotspot(LOD of score) is exist
        if "LS" in record.INFO: 
            vcf_writer.write_record(record)

        # TODO EBCall
        # FP:  Fisher's Test
        # FPR: Fisher's Test processed with Realignment
        # B10: 10% posterior quantile 
        # B1R: 10% posterior quantile Processed with Realignment
        # NAR: Number of allelic reads (Tumor)|(Sample1)
        # NAR: Number of allelic reads (Normal)|(Sample2)
        
        
        else:
            # single mode 
            if sample2 == None:
            
                if \
                ( "EB" not in record.INFO or record.INFO["EB"] >= ebpval) and \
                ( record.INFO["B10"] >= float(post10q))  and \
                ( record.INFO["B1R"] != None and (record.INFO["B1R"] >= r_post10q))  and \
                ( record.genotype(sample1)["NAR"] != None and record.genotype(sample1)["NAR"] >= v_count):
                    vcf_writer.write_record(record)
            
            # pair mode
            else:
        
                if \
                ( "EB" not in record.INFO or record.INFO["EB"] >= ebpval) and \
                ( record.INFO["FP"] >= float(fishpval) or \
                    (record.genotype(sample2)["AD"] == 0 and \
                     flag_mis_base_0 == True and \
                     record.INFO["FP"] >= fishpval_base_0) \
                ) and \
                ( record.INFO["FPR"] != None and \
                    (record.INFO["FPR"] >= realignpval or \
                        (record.genotype(sample2)["NAR"] == 0 and \
                        flag_mis_base_0 == True and \
                        record.INFO["FPR"] >= realignpval_base_0) \
                    )
                ) and \
                ( record.genotype(sample1)["NAR"] != None and record.genotype(sample1)["NAR"] >= tcount)   and \
                ( record.genotype(sample2)["NAR"] != None and record.genotype(sample2)["NAR"] <= ncount)   and \
                ( record.genotype(sample2)["NOR"] != None and record.genotype(sample2)["NOR"] <= ncount_other)  and \
                ( record.genotype(sample2)["NNR"] != None and (record.genotype(sample2)["NNR"] + record.genotype(sample2)["NAR"]) >= int(ncount_depth)) \
                :
                    vcf_writer.write_record(record)

    vcf_writer.close()
    

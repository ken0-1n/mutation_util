
#! /usr/bin/env python

class Maf_header_info:

    def __init__(self):

        self.ncbi_build = -1
        self.chr = -1
        self.start = -1
        self.end = -1
        self.variant_class = -1
        self.variant_type = -1
        self.ref = -1
        self.t_allele1 = -1
        self.t_allele2 = -1
        self.t_barcode = -1
        self.t_ref_count = -1
        self.t_alt_count = -1
        self.merge_status = -1

    def set_header_information(self,header):
        F = header.rstrip('\n').split('\t')
        for i, v in enumerate(F):
            if v == "NCBI_Build":
                self.ncbi_build = i
            elif v == "Chromosome":
                self.chr = i
            elif v == "Start_position":
                self.start = i
            elif v == "End_position":
                self.end = i
            elif v == "Variant_Classification":
                self.variant_class = i
            elif v == "Variant_Type":
                self.variant_type = i
            elif v == "Reference_Allele":
                self.ref = i
            elif v == "Tumor_Seq_Allele1":
                self.t_allele1 = i
            elif v == "Tumor_Seq_Allele2":
                self.t_allele2 = i 
            elif v == "Tumor_Sample_Barcode":
                self.t_barcode = i
            elif v == "t_ref_count":
                self.t_ref_count = i 
            elif v == "t_alt_count":
                self.t_alt_count = i
            elif v == "merge_status":
                self.merge_status = i
    
    

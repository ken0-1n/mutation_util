
#! /usr/bin/env python

class Genomon_header_info:

    def __init__(self):
    
        self.fisher_idx = -1
        self.ebcall_idx = -1
        self.realign_idx = -1
        self.tcount_idx = -1
        self.ncount_idx = -1
        self.post10q_idx = -1
        self.r_post10q_idx = -1
        self.v_count_idx = -1
        self.chr_idx = -1
        self.start_idx = -1
        self.end_idx = -1
        self.ref_idx = -1
        self.alt_idx = -1
        self.func_idx = -1
        self.gene_idx = -1
   
    def set_header_information(self,header):
        F = header.rstrip('\n').split('\t')
        for i, v in enumerate(F):
            if v == "P-value(fisher)":
                self.fisher_idx = i
            elif v == "P-value(EBCall)":
                self.ebcall_idx = i
            elif v == "P-value(fisher_realignment)":
                self.realign_idx = i
            elif v == "variantPairNum_tumor" :
                self.tcount_idx = i
            elif v == "variantPairNum_normal":
                self.ncount_idx = i
            elif v == "10%_posterior_quantile":
                self.post10q_idx = i
            elif v == "10%_posterior_quantile(realignment)":
                self.r_post10q_idx = i
            elif v == "variantPairNum":
                self.v_count_idx = i
            elif v == "Chr":
                self.chr_idx = i
            elif v == "Start":
                self.start_idx = i
            elif v == "End":
                self.end_idx = i
            elif v == "Ref":
                self.ref_idx = i
            elif v == "Alt":
                self.alt_idx = i
            elif v == "Func.refGene":
                self.func_idx = i
            elif v == "Gene.refGene":
                self.gene_idx = i
    

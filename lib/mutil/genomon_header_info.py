
#! /usr/bin/env python

class Genomon_header_info:

    def __init__(self):
    
        self.fisher = -1
        self.ebcall = -1
        self.realign = -1
        self.tcount = -1
        self.ncount = -1
        self.post10q = -1
        self.r_post10q = -1
        self.v_count = -1
        self.chr = -1
        self.start = -1
        self.end = -1
        self.ref = -1
        self.alt = -1
        self.func = -1
        self.gene = -1
        self.tumor_barcode = -1
        self.tdepth = -1
        self.tvariant = -1
        self.ndepth = -1
        self.nvariant = -1
        self.score_hotspot = -1
   
    def set_header_information(self,header):
        F = header.rstrip('\n').split('\t')
        for i, v in enumerate(F):
            if v == "P-value(fisher)":
                self.fisher = i
            elif v == "P-value(EBCall)":
                self.ebcall = i
            elif v == "P-value(fisher_realignment)":
                self.realign = i
            elif v == "variantPairNum_tumor" :
                self.tcount = i
            elif v == "variantPairNum_normal":
                self.ncount = i
            elif v == "10%_posterior_quantile":
                self.post10q = i
            elif v == "10%_posterior_quantile(realignment)":
                self.r_post10q = i
            elif v == "variantPairNum":
                self.v_count = i
            elif v == "Chr":
                self.chr = i
            elif v == "Start":
                self.start = i
            elif v == "End":
                self.end = i
            elif v == "Ref":
                self.ref = i
            elif v == "Alt":
                self.alt = i
            elif v == "Func.refGene":
                self.func = i
            elif v == "Gene.refGene":
                self.gene = i
            elif v == "Tumor_Sample_Barcode":
                self.tumor_barcode = i
            elif v == "depth_tumor":
                self.tdepth = i
            elif v == "variantNum_tumor":
                self.tvariant = i
            elif v == "depth_normal":
                self.ndepth = i
            elif v == "variantNum_normal":
                self.nvariant = i
            elif v == "score(hotspot)":
                self.score_hotspot = i
    

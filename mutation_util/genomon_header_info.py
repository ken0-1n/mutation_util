class Genomon_header_info:

    def __init__(self):
    
        self.fisher = None
        self.ebcall = None
        self.realign = None
        self.tcount = None
        self.ncount = None
        self.ncount_ref = None
        self.ncount_other = None
        self.post10q = None
        self.r_post10q = None
        self.v_count = None
        self.chr = None
        self.start = None
        self.end = None
        self.ref = None
        self.alt = None
        self.func = None
        self.gene = None
        self.tumor_barcode = None
        self.tdepth = None
        self.tvariant = None
        self.ndepth = None
        self.nvariant = None
        self.score_hotspot = None
   
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
            elif v == "readPairNum_normal":
                self.ncount_ref= i
            elif v == "otherPairNum_normal":
                self.ncount_other = i
            elif v == "10%_posterior_quantile":
                self.post10q = i
            elif v == "10%_posterior_quantile(realignment)":
                self.r_post10q = i
            elif v == "variantPairNum":
                self.v_count = i
            elif v == "Chr" or v == "#chr":
                self.chr = i
            elif v == "Start" or v == "start":
                self.start = i
            elif v == "End" or v == "end":
                self.end = i
            elif v == "Ref" or v == "ref":
                self.ref = i
            elif v == "Alt" or v == "alt":
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
    

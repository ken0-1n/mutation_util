import sys
import os
import argparse
from . import mutation_filter
from . import hotspot_merge
from . import genomon_header_info

def run_filter(arg):
    ghi = genomon_header_info.Genomon_header_info()
    is_anno = True if arg.print_format == 'anno' else False
    if is_anno == True:
        mutation_filter.filter_mutation_list(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.count, arg.hotspot_db, ghi, arg.flag_mis_base_0)
    else:
        mutation_filter.filter_mutation_vcf(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.sample1, arg.sample2, ghi, arg.flag_mis_base_0)

def run_merge_hotspot(arg):
    is_vcf = True if arg.print_format == 'vcf' else False
    if is_vcf == True:
        hotspot_merge.merge_hotspot_vcf(arg.input_vcf, arg.hotspot_vcf, arg.merged_vcf)
    else:
        print("The anno format is not supported yet.",file=sys.stderr)


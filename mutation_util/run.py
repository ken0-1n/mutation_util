import sys
import os
import argparse
from . import mutation_filter
from . import genomon_header_info

def run_filter(arg):
    ghi = genomon_header_info.Genomon_header_info()
    is_anno = True if arg.print_format == 'anno' else False
    if is_anno == True:
        mutation_filter.filter_mutation_list(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.count, arg.hotspot_db, ghi)
    else:
        mutation_filter.filter_mutation_vcf(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.sample1, arg.sample2, ghi)

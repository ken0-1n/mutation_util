#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
import mutation_util
import count_bases


#
# Main
#
def run_compare(arg):
    mutation_util.compare_list(arg.input, arg.output, arg.database, arg.mapchain, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.gene_ref, arg.post10q, arg.r_post10q, arg.count, arg.print_graph, arg.pancan, arg.hotspot)

def run_filter(arg):
    mutation_util.filt_mutation_result(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.count)

def run_all(arg):
    mutation_util.compare_all(arg.input, arg.output, arg.database, arg.mapchain, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.gene_ref, arg.post10q, arg.r_post10q, arg.count, arg.pancan, arg.hotspot)

def run_count(arg):
    count_bases.pileup_and_count(arg.inbam, arg.ref_fa, arg.samtools_path, arg.region)


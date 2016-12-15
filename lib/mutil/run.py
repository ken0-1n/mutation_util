#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
import mutation_util
import count_bases
import hotspot_check
import add_annotation
import blacklist
import hotspot_merge


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

def run_hotspot(arg):
    hotspot_check.compare_hotspot_list(arg.in_hotspot_mutation, arg.in_genomon_mutation, arg.output, arg.noheader)

def run_annotate(arg):
    add_annotation.annotate(arg.in_mutation, arg.annovar_path, arg.noheader) 

def run_blacklist(arg):
    blacklist.filter(arg.in_mutation, arg.blacklist, arg.output, arg.min_candidate) 
    
def run_hotspot_merge(arg):
    hotspot_merge.merge_hotspot_list(arg.in_hotspot_mutation, arg.in_fisher_mutation, arg.output_file, arg.hotspot_header, arg.fisher_header)



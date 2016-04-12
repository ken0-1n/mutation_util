#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
import mutation_util


#
# Main
#
def run_compare(arg):

    mutation_util.compare_list(arg.input, arg.output, arg.database, arg.mapchain, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.func_ref, arg.gene_ref, arg.post10q, arg.r_post10q, arg.count)

def run_filter(arg):

    mutation_util.filt_mutation_result(arg.input, arg.output, arg.eb_pval, arg.fish_pval, arg.realign_pval, arg.tcount, arg.ncount, arg.post10q, arg.r_post10q, arg.count)


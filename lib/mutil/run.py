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

    mutation_util.compare_list(arg.input, arg.output, arg.database, arg.mapchain, arg.ebpval, arg.fishpval, arg.realignpval, arg.tcount, arg.ncount, arg.func_ref, arg.gene_ref)


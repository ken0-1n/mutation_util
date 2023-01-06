#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .run import run_filter
from .run import run_merge_hotspot

def create_parser():
    prog = "mutil"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_filter_parser(subparsers):
        
        filter_parser = subparsers.add_parser("filter")
        filter_parser.add_argument( '-i', '--input', help = 'mutation list', type = str, default = None, required = True )
        filter_parser.add_argument( '-o', '--output', help = 'output directory', type = str, default = None, required = True)
        filter_parser.add_argument( '-e', '--eb_pval', help = 'EBCall P-value', type = float, default = 0.0)
        filter_parser.add_argument( '-f', '--fish_pval', help = 'Fisher test P-value', type = float, default = 0.0)
        filter_parser.add_argument( '-r', '--realign_pval', help = 'Realignment Fisher test P-value', type = float, default = 0.0)
        filter_parser.add_argument( '-t', '--tcount', help = 'read count of tumor', type = int, default = -1)
        filter_parser.add_argument( '-n', '--ncount', help = 'read count of normal', type = int, default = 1000000)
        filter_parser.add_argument( '-p', '--post10q', help = 'Fisher 10 percent posterior quantile', type = float, default = 0.0)
        filter_parser.add_argument( '-q', '--r_post10q', help = 'Realignment 10 percent posterior quantile', type = float, default = 0.0)
        filter_parser.add_argument( '-c', '--count', help = 'read count', type = int, default = -1)
        filter_parser.add_argument( '-d', '--hotspot_db', help = 'hotspot_database_file', type = str, default = "")
        filter_parser.add_argument( '-1', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        filter_parser.add_argument( '-2', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        filter_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        filter_parser.add_argument( '-D', '--flag_mis_base_0', help = 'Ignore the fisher threshold when the number of mismatches is 0 in normal bam.', action = 'store_true', default = False )
        filter_parser.add_argument( '-F', '--fish_pval_base_0', help = 'Fisher test P-value when the number of mismatches is 0 in normal bam.', type = float, default = 0.0)
        filter_parser.add_argument( '-R', '--realign_pval_base_0', help = 'Realignment Fisher test P-value when the number of mismatches is 0 in normal bam.', type = float, default = 0.0)
        filter_parser.add_argument( '-N', '--ncount_depth', help = "minimum depth of normal", type = int, default = 0)
        filter_parser.add_argument( '-E', '--ncount_other', help = 'maximum other read count of normal', type = int, default = 1000000)
        

        return filter_parser
        
    def _create_merge_hotspot_parser(subparsers):
        
        merge_hotspot_parser = subparsers.add_parser("merge_hotspot")
        merge_hotspot_parser.add_argument( '-i', '--input_vcf', help = 'path of input vcf is bgziped', type = str, default = None, required = True )
        merge_hotspot_parser.add_argument( '-d', '--hotspot_vcf', help = 'path of hotspot vcf is bgziped', type = str, default = None, required = True)
        merge_hotspot_parser.add_argument( '-o', '--merged_vcf', help = 'path of output vcf. Do not use the gz extension', type = str, default = None, required = True)
        merge_hotspot_parser.add_argument( '-O', '--print_format', choices = ['vcf'], help = 'Print VCF or anno(TSV) format',  default = 'vcf' )
        return merge_hotspot_parser
        
    filter_parser = _create_filter_parser(subparsers)
    filter_parser.set_defaults(func = run_filter)
    merge_hotspot_parser = _create_merge_hotspot_parser(subparsers)
    merge_hotspot_parser.set_defaults(func = run_merge_hotspot)
    return parser

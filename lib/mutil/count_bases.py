#!/usr/bin/python
"""

fisher.py


"""
import sys
import os
import re
import pysam
import scipy.special
from scipy.stats import fisher_exact as fisher
import argparse
import logging
import math
import subprocess

#
# Globals
#
arg = None
target = None
remove_chr = None

#
# Class definitions
#

############################################################
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)

        except KeyError:
            value = self[item] = type(self)()

        return value

############################################################
def pileup_out( mpileup):

    global target
    global remove_chr

    #
    # Prepare mpileup data
    #
    mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
    mp_list_len = len( mp_list )
    ref_base_U = mp_list[ 2 ].upper()
    coordinate = mp_list[ 0:3 ]

    #
    # skip if depth is 0
    #
    if mp_list[ 3 ] == '0':
        return None

    chr =  mp_list[0]
    pos =  mp_list[1]
    depth = mp_list[3]
    read_bases = mp_list[4]
    qual_list = mp_list[5]

    indel = AutoVivification()
    #
    deleted = 0
    iter = target.finditer( read_bases )
    for m in iter:
        site = m.start()
        type = m.group( 1 )
        num = m.group( 2 )
        bases = m.group( 3 )[ 0:int( num ) ]
        if bases.islower():
            strand = ( '-', '+' )
        else:
            strand = ( '+', '-' )
    
        key = '\t'.join( coordinate + [ bases.upper() ] )
        if type in indel and key in indel[ type ]:
            indel[ type ][ key ][ strand[ 0 ] ] += 1
        else:
            indel[ type ][ key ][ strand[ 0 ] ] = 1
            indel[ type ][ key ][ strand[ 1 ] ] = 0
    
        read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
        deleted += 1 + len( num ) + int( num )
    
    #
    # Remove '^.' and '$'
    #
    read_bases = remove_chr.sub( '', read_bases )
    read_bases = read_bases.translate( None, '$' ) 

    #
    # Error check
    #
    if len( read_bases ) != len( qual_list ):
        logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
        return None

    # Count mismatch
    #
    read_bases = read_bases.replace( '.', ref_base_U )
    read_bases = read_bases.replace( ',', ref_base_U.lower() )
    
    depth = 0
    base_num = {
                "A": 0,
                "C": 0,
                "G": 0,
                "T": 0,
                "N": 0,
                "a": 0,
                "c": 0,
                "g": 0,
                "t": 0,
                "n": 0
               }
    for nuc, qual in zip( read_bases, qual_list ):
        if nuc in 'ATGCatgc':
            base_num[ nuc ] += 1
            depth += 1

    mis_num = 0
    mis_base_U = ref_base_U
    for nuc in ( 'A', 'C', 'G', 'T' ):
        if nuc != ref_base_U:
            tmp_mis_num = int(base_num[nuc]) + int(base_num[nuc.lower()]) 
            if tmp_mis_num > mis_num:
                mis_base_U = nuc
                mis_num = tmp_mis_num

    if mis_base_U == ref_base_U: return

    mismatch_qual_list = []
    reference_qual_list = []
    for nuc, qual in zip( read_bases, qual_list ):
        if nuc == mis_base_U or nuc == mis_base_U.lower():
            mismatch_qual_list.append(int(ord(qual))-33)
        if nuc == ref_base_U or nuc == ref_base_U.lower():
            reference_qual_list.append(int(ord(qual))-33)

    # print chr +"\t"+ pos +"\t"+ ref_base_U +"\t"+ mis_base_U +"\t"+ ",".join(map(str,reference_qual_list)) +"\t"+",".join(map(str,mismatch_qual_list))
    print chr +"\t"+ pos +"\t"+ ref_base_U +"\t"+ mis_base_U +"\t"+str(depth)+"\t"+ ",".join(map(str,reference_qual_list)) +"\t"+",".join(map(str,mismatch_qual_list))



############################################################
def pileup_and_count(
        in_bam,
        ref_fa,
        samtools,
        region
        ):

    global target
    global remove_chr

    #
    # Setup regular expression
    # ([\+\-])[0-9]+[ACGTNacgtn]+
    #
    target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
    remove_chr = re.compile( '\^.' )

    #
    # Open output file and write header
    #
    FNULL = open(os.devnull, 'w')
    #
    # Print header only for testing.
    #

    cmd_list = [samtools,'mpileup','-r',region,'-q','0','-BQ','0','-d','10000000','-f',ref_fa, in_bam ]
    pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
    end_of_pipe = pileup.stdout
    for mpileup in end_of_pipe:
        pileup_out( mpileup)

    FNULL.close()


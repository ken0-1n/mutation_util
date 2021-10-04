#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
from mutation_util import hotspot_merge as hm

class TestHotspotMerge(unittest.TestCase):

    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        in_vcf = cur_dir + "/../data/5929_hotspot_merge_test1-1.vcf.gz"
        hotspot_vcf = cur_dir + "/../data/5929_hotspot_output_test1-1.vcf.gz"
        output_vcf = cur_dir + "/../data/output_test1-1.vcf"

        hm.merge_hotspot_vcf( \
            in_vcf, \
            hotspot_vcf, \
            output_vcf)

        answer_file = cur_dir + "/../data/output_answer_test1-1.vcf"
        self.assertTrue(filecmp.cmp(output_vcf, answer_file, shallow=False))

    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        in_vcf = cur_dir + "/../data/5929_hotspot_merge_test1-2.vcf.gz"
        hotspot_vcf = cur_dir + "/../data/5929_hotspot_output_test1-2.vcf.gz"
        output_vcf = cur_dir + "/../data/output_test1-2.vcf"

        hm.merge_hotspot_vcf( \
            in_vcf, \
            hotspot_vcf, \
            output_vcf)

        answer_file = cur_dir + "/../data/output_answer_test1-2.vcf"
        self.assertTrue(filecmp.cmp(output_vcf, answer_file, shallow=False))

    def test1_3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        in_vcf = cur_dir + "/../data/5929_hotspot_merge_test1-3.vcf.gz"
        hotspot_vcf = cur_dir + "/../data/5929_hotspot_output_test1-3.vcf.gz"
        output_vcf = cur_dir + "/../data/output_test1-3.vcf"

        hm.merge_hotspot_vcf( \
            in_vcf, \
            hotspot_vcf, \
            output_vcf)

        answer_file = cur_dir + "/../data/output_answer_test1-3.vcf"
        self.assertTrue(filecmp.cmp(output_vcf, answer_file, shallow=False))

    def test2_1(self):
        
        cmd = ['mutil', '--version']
        subprocess.check_call(cmd)
        
        
    def test2_2(self):
        
        cmd = ['mutil', 'merge_hotspot', '--help']
        subprocess.check_call(cmd)
        
        
if __name__ == "__main__":
    unittest.main()


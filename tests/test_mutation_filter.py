#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
import subprocess
from mutation_util import mutation_filter as mf
from mutation_util import genomon_header_info

class TestMutationFilter(unittest.TestCase):

    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--fish_pval 1.0 --realign_pval 1.0 --eb_pval 4.0 --tcount 4 --ncount 2"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test1.txt"
        output = cur_dir + "/../data/5929_small_result_test1_1.txt"
        ebpval = 4.0 
        fishpval = 1.0
        realignpval = 1.0
        tcount = 4
        ncount = 2
        post10q = 100
        r_post10q = 100
        v_count = 100
        hotspot_database = ""
        flag_mis_base_0 = True
            
        mf.filter_mutation_list( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            hotspot_database, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test1_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test1_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--fish_pval 1.0 --realign_pval 1.0 --eb_pval 4.0 --tcount 4 --ncount 2"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test1.2.txt"
        output = cur_dir + "/../data/5929_small_result_test1_2.txt"
        ebpval = 4.0 
        fishpval = 1.0
        realignpval = 1.0
        tcount = 4
        ncount = 2
        post10q = 100
        r_post10q = 100
        v_count = 100
        hotspot_database = ""
        flag_mis_base_0 = True
            
        mf.filter_mutation_list( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            hotspot_database, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test1_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--post10q 0.1 --r_post10q 0.1 --count 4"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test2.txt"
        output = cur_dir + "/../data/5929_small_result_test2_1.txt"
        ebpval = 100 
        fishpval = 100
        realignpval = 100
        tcount = 100
        ncount = 100
        post10q = 0.1
        r_post10q = 0.1
        v_count = 4
        hotspot_database = ""
        flag_mis_base_0 = True
            
        mf.filter_mutation_list( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            hotspot_database, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test2_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test2_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--post10q 0.1 --r_post10q 0.1 --count 4"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test2.2.txt"
        output = cur_dir + "/../data/5929_small_result_answer_test2_2.txt"
        ebpval = 100 
        fishpval = 100
        realignpval = 100
        tcount = 100
        ncount = 100
        post10q = 0.1
        r_post10q = 0.1
        v_count = 4
        hotspot_database = ""
        flag_mis_base_0 = True
            
        mf.filter_mutation_list( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            hotspot_database, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test2_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    ######################################
    # Tumor/Normal Pair, VCF format
    ######################################
    def test3_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--fish_pval 1.0 --realign_pval 1.0 --eb_pval 4.0 --tcount 4 --ncount 2"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test3.txt"
        output = cur_dir + "/../data/5929_small_result_test3_1.txt"
        ebpval = 4.0 
        fishpval = 1.0
        realignpval = 1.0
        tcount = 4
        ncount = 2
        v_count = 4
        post10q = 100
        r_post10q = 100
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        flag_mis_base_0 = True
        
        mf.filter_mutation_vcf( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            sample1, \
            sample2, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test3_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test3_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--fish_pval 1.0 --realign_pval 1.0 --eb_pval 4.0 --tcount 4 --ncount 2"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test3.2.txt"
        output = cur_dir + "/../data/5929_small_result_test3_2.txt"
        ebpval = 4.0 
        fishpval = 1.0
        realignpval = 1.0
        tcount = 4
        ncount = 2
        v_count = 4
        post10q = 100
        r_post10q = 100
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        flag_mis_base_0 = True
        
        mf.filter_mutation_vcf( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            sample1, \
            sample2, \
            ghi, \
            flag_mis_base_0)

        answer_file = cur_dir + "/../data/5929_small_result_answer_test3_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))

    def test4_1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--post10q 0.1 --r_post10q 0.1 --count 4"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test4.txt"
        output = cur_dir + "/../data/5929_small_result_test4_1.txt"
        ebpval = 100 
        fishpval = 100
        realignpval = 100
        tcount = 4
        ncount = 100
        post10q = 0.1
        r_post10q = 0.1
        v_count = 4
        sample1 = "5929_tumor"
        sample2 = None
        flag_mis_base_0 = True
        
        mf.filter_mutation_vcf( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            sample1, \
            sample2, \
            ghi, \
            flag_mis_base_0)


        answer_file = cur_dir + "/../data/5929_small_result_test4_1.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test4_2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        "--post10q 0.1 --r_post10q 0.1 --count 4"
        ghi = genomon_header_info.Genomon_header_info()
        input = cur_dir + "/../data/5929_small_realignment_result_test4.2.txt"
        output = cur_dir + "/../data/5929_small_result_test4_2.txt"
        ebpval = 100 
        fishpval = 100
        realignpval = 100
        tcount = 4
        ncount = 100
        post10q = 0.1
        r_post10q = 0.1
        v_count = 4
        sample1 = "5929_tumor"
        sample2 = None
        flag_mis_base_0 = True
        
        mf.filter_mutation_vcf( \
            input, \
            output, \
            ebpval, \
            fishpval, \
            realignpval, \
            tcount, \
            ncount, \
            post10q, \
            r_post10q, \
            v_count, \
            sample1, \
            sample2, \
            ghi, \
            flag_mis_base_0)


        answer_file = cur_dir + "/../data/5929_small_result_test4_2.txt"
        self.assertTrue(filecmp.cmp(output, answer_file, shallow=False))


    def test5_1(self):
        
        cmd = ['mutil', '--version']
        subprocess.check_call(cmd)
        
        
    def test5_2(self):
        
        cmd = ['mutil', 'filter', '--help']
        subprocess.check_call(cmd)
        
        
if __name__ == "__main__":
    unittest.main()


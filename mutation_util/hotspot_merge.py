import sys, os, subprocess

def merge_hotspot_vcf(in_vcf, hotspot_vcf, merged_vcf):

    out_prefix, ext = os.path.splitext(merged_vcf)

    subprocess.check_call(['bcftools', 'isec', '-p', out_prefix, in_vcf, hotspot_vcf]) 

    subprocess.check_call(['bgzip', '-f', out_prefix+'/0001.vcf']) 

    subprocess.check_call(['tabix', '-p', 'vcf', out_prefix+'/0001.vcf.gz']) 
    
    subprocess.check_call(['bcftools', 'concat', '-a', '-O', 'vcf', '-o', merged_vcf, in_vcf, out_prefix+'/0001.vcf.gz']) 

    subprocess.check_call(['bgzip', '-f', merged_vcf]) 

    subprocess.check_call(['tabix', '-p', 'vcf', merged_vcf+'.gz']) 


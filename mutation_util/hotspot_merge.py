import sys, os, subprocess, gzip

def merge_hotspot_vcf(in_vcf, hotspot_vcf, merged_vcf):

    out_prefix, ext = os.path.splitext(merged_vcf)
    subprocess.check_call(['bcftools', 'isec', '-p', out_prefix, in_vcf, hotspot_vcf]) 


    with open(out_prefix+'_sorted.vcf', 'w') as hout:
        with gzip.open(in_vcf, 'rt') as hin:
            for line in hin:
                line = line.rstrip('\n')
                if line.startswith("##"):
                    print(line, file=hout)
                elif line.startswith('#'):
                    print('##INFO=<ID=LS,Number=1,Type=Float,Description="LOD Score of Hotspot Call">',file=hout)
                    print(line, file=hout)
                else:
                    print(line, file=hout)

        with open(out_prefix+'/0001.vcf', 'r') as hin:
            for line in hin:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                print(line, file=hout)

    subprocess.check_call(['bcftools', 'sort', '-O', 'vcf', '-o', merged_vcf, out_prefix+'_sorted.vcf']) 
    os.remove(out_prefix+'/0000.vcf')
    os.remove(out_prefix+'/0001.vcf')
    os.remove(out_prefix+'/0002.vcf')
    os.remove(out_prefix+'/0003.vcf')
    os.remove(out_prefix+'/README.txt')
    os.rmdir(out_prefix)
    os.remove(out_prefix+'_sorted.vcf')


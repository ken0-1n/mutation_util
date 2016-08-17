import sys, os, subprocess

def annotate(in_mutation, annovar_path, is_header):

    base, ext = os.path.splitext(in_mutation)

    meta_list = []
    header_list = []
    header_print_flg = True
    hout = open(base, 'w')
    with open(in_mutation, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                line = line.rstrip('\n')
                meta_list.append(line)
                continue
            if is_header and header_print_flg:
                line = line.rstrip('\n')
                header_list = line.split("\t")
                header_print_flg = False
                continue

            print >> hout, line
    hout.close()

    subprocess.call([annovar_path + "/table_annovar.pl", "--otherinfo", "--outfile", base, '-buildver','hg19', '-remove', '-protocol', 'refGene', '-operation', 'g', base, annovar_path+'/humandb'])

    hout = open(base+".result.txt", 'w')
    header_print_flg = True
    with open(base + ".hg19_multianno.txt", 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if is_header and header_print_flg:
                line = line.replace('Otherinfo','')
                print >> hout, line + "\t".join(header_list[5:])
                header_print_flg = False
                continue
            print >> hout, line
    hout.close()

    os.unlink(base)



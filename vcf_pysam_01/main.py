import vcf_pysam_helper
from datetime import datetime
import sys
import psutil
import os


def usage():
    process = psutil.Process(os.getpid())
    return process.memory_info()[0] / float(2 ** 20)

# ======================================================================


print 'memory usage: ', usage()
start = datetime.now()

# modify your input file here
# filename_in = '../../data/CAM_Realigned_snps.vcf'
filename_in = '../../data/sample_vcf_file.vcf'
# filename_in = '../../data/Pfeiffer.vcf'
filename_out = filename_in[0:len(filename_in)-3] + 'csv'

vcf_pysam_helper.vcf_get_all(filename_in, filename_out)
print 'Total runtime: ', datetime.now() - start

print 'memory usage: ', usage()

# ======================================================================


def main(filename_i):

    sta = datetime.now()
    filename_o = filename_i[0:len(filename_i)-3] + 'csv'

    vcf_pysam_helper.vcf_get_all(filename_i, filename_o)
    print 'Total runtime: ', datetime.now() - sta

    print 'memory usage: ', usage()
    

sys.exit(0)

import vcf_pysam_helper
from datetime import datetime
import sys
import gc
import os, psutil 

def usage():
    process = psutil.Process(os.getpid())
    return process.memory_info()[0] / float(2 ** 20)

##########################################################################    
print 'memory usage: ', usage()
start=datetime.now()

# modify your input file here
#filename_in = '../../data/CAM_Realigned_snps.vcf'
filename_in = '/sample_vcf_file.vcf'
#filename_in = '../../data/Pfeiffer.vcf'
filename_out = filename_in[0:len(filename_in)-3] + 'csv'

vcf_pysam_helper.vcf_get_all(filename_in, filename_out)
print 'Total runtime: ', datetime.now() - start

print 'memory usage: ',usage()

##########################################################################
def main(filename_in):
    start=datetime.now()
    filename_out = filename_in[0:len(filename_in)-3] + 'csv'
    start=datetime.now()

    vcf_pysam_helper.vcf_get_all(filename_in, filename_out)
    print 'Total runtime: ', datetime.now() - start

    print 'memory usage: ',usage()
    
sys.exit(0)

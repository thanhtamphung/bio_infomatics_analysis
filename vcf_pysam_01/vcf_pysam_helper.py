import pandas as pd
import numpy as np
from pysam import VariantFile as psVa
from datetime import datetime

#############################################################################
def vcf_get_record(filename):
    # read vcf file using pysam.VariantFile package
    # input = path of vcf file
    # output to list record (each element is like a row in dataset)
    
    vcf_reader = psVa(filename)
    record = []
    for i in vcf_reader:
        record.append(i)
    return record


#############################################################################
def vcf_get_all(filename_in, filename_out):
    # convert file vcf to any table file, for example csv here
    # input = path of vcf file
    # output = path of csv file
    
    # read vcf file
    start = datetime.now() 
    print 'Importing vcf file ...'
    record = vcf_get_record(filename_in)
    print 'Runtime of reading:', datetime.now()-start
    
    n_sample = len(record[0].samples)
    n = len(record)
    print 'Number of rows:', n
    
    print 'Converting to csv ...'
    start = datetime.now()
    
    # parse record_list to a dataframe
    df = pd.DataFrame(index = np.arange(n))
    tmp1 = []
    tmp2 = []
    tmp3 = []
    tmp4 = []
    tmp5 = []
    tmp6 = []
    tmp7 = []
    tmp8 = []
    tmp9 = []
    msamples = []
    for j in range(n_sample):
        msamples.append([])
    for i in range(n):
        tmp1.append(record[i].chrom)
        tmp2.append(record[i].pos)
        tmp3.append(record[i].id)
        tmp4.append(record[i].ref)
        tmp5.append(record[i].alts)
        tmp6.append(record[i].qual)
        tmp7.append(record[i].filter)
        tmp8.append(record[i].info)
        tmp9.append(record[i].format)
        for j in range(n_sample):
            msamples[j].append(record[i].samples[j].values())
    df['CHROM'] = tmp1
    df['POS'] = tmp2
    df['ID'] = tmp3
    df['REF'] = tmp4
    df['ALT'] = tmp5
    df['QUAL'] = tmp6
    df['FILTER'] = tmp7
    df['INFO'] = tmp8
    df['FORMAT'] = tmp9
    for j in range(n_sample):
        df[record[0].samples[j].name] = msamples[j]  
        
    # export dataframe to table file    
    df.to_csv(filename_out)
    print 'Runtime of converting:', datetime.now()-start
    print 'Finish!'
    return df

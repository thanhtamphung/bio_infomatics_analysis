# coding: utf-8

# # Burden test:
#     Date: 13 mar 2019
#     Last Edit: 29 May 2019
#     Install all dependencies
#     Ref: Rvtest https://genome.sph.umich.edu/wiki/Rvtests 
#          Exomiser https://www.sanger.ac.uk/science/tools/exomiser
#     Run cells step by step
#     Runtime: 30 minutes with input vcf 120MB
#     

# ## Update:
# Burden test with kernel SKAT
# 29 May: check spelling

# ## Dependencies:
# 
#     1/ python packages in import cell (pandas, numpy, pysam)
#     2/ bedtools
#     3/ vcftools (vcf-subset, vcf-sort, vcf-merge)
#     4/ htslib (bgzip, tabix)
#     5/ exomiser 8.0.0 https://github.com/exomiser/Exomiser/releases/tag/8.0.0
#     6/ rvtests https://github.com/zhanxw/rvtests

# # Burden test Framework

# ## 1/ Review data
# ## 2/ Filter position with bed file
# ## 3/ Split sample from vcf
# ## 4/ Filter read depth (min_base_count = 10) and maybe ref, alt
# ## 5/ Annotation with exomiser
# ## 6/ Filter function_class (protein affecting) and maybe freq
# ## 7/ Merge results
# ## 8/ Preprocessing and Burden test
# ## 9/ Output report
# 
# 

# ## =========================================================================

# In[1]:


# import cell
import subprocess
import os
import fileinput
import sys

import pandas as pd
from pysam import VariantFile as psVa

# ## ============= Modify these parameters =================================

# In[2]:


# folder containing all input: vcf, bedfile, correspondance_table, exome file
input_path = "/data/BurdenTest/Lithiasis_mar/"

# vcf and input file name. All input files must be in input folder
vcf_name = "Lithiasis_set_20082018-123-samples.vcf"
bedfile = "BED_Letav-sorted.bed"
ctable = "Correspondance_table.tsv"
exomefile = "exome.yml"
genefile = input_path + "refFlat.txt.gz"

# List of protein affection for filtering
af_pro = ["missense_variant", "splice_region_variant", "stop_gained", "frameshift_variant", "stop_lost",
          "inframe_deletion", "inframe_insertion"]

# Tool path: 
exom_path = "/home/phung/Documents/work/exomiser/exomiser_cli_8/exomiser-cli-8.0.0.jar"
rvtest_pwd = "/data/rvtests/executable/rvtest"

# ## ============= Run cells and do not modify ===========================

# In[3]:

print "Generating folders for output ..."
output = input_path + "output"
out_data = output + "/data"
out_exome = output + "/exome"
out_rvtest = output + "/rvtest"
out_report = output + "/report"
os.system("mkdir " + output)
os.system("mkdir " + out_data)
os.system("mkdir " + out_exome)
os.system("mkdir " + out_rvtest)
os.system("mkdir " + out_report)
print "***** All outputs are in the folder: ", output

# In[4]:


input_vcf = input_path + vcf_name
input_bed = input_path + bedfile
input_table = input_path + ctable
input_exome = input_path + exomefile

# In[5]:


name = vcf_name[0:len(vcf_name)-4]

# # ====================== 1/ Review input data ===========================================

# In[6]:
print "================================================================================"
print "========================= 1/ Review input data ================================="
print "================================================================================"
a = subprocess.check_output("grep -v '^#' " + input_vcf + " | cut -f 1 | sort", shell=True)
print "Number of variants (base pairs): \n", a.count('\n')

# In[7]:


a = subprocess.check_output("grep -v '^#' " + input_vcf + " | cut -f 1 | sort | uniq -c", shell=True)
print "Number of chrome: \n", a.count('\n')
print "Number of variants in each chrome: \n", a

# In[8]:


vcf_reader = psVa(input_vcf)
l_samples = list(vcf_reader.header.samples)
print "Number of samples:", len(l_samples)
print "Name of samples:\n", l_samples

# # ====================  2/ Filter with bed file  ======================================

# In[9]:
print "================================================================================"
print "===================== 2/ Filter with bed file =================================="
print "================================================================================"
# compare chrome names in vcf and bed file. Make sure the same format (or replace("X", "GL000192.1", 1))
df_bed = pd.read_csv(input_bed, sep="\t")
df_bed.head()

# In[10]:


vcf_bedfil = out_data + '/' + name + '_filterposition.vcf'
cmn = "bedtools intersect -a " + input_vcf + " -b " + input_bed + " -header > " + vcf_bedfil
print subprocess.check_output(cmn, shell=True)
print "***** Filter position with bedtools. Output is ", vcf_bedfil

# In[11]:


# check number of variants after filter position 
a = subprocess.check_output("grep -v '^#' " + vcf_bedfil + " | cut -f 1 | sort", shell=True)
print "Number of variants after filtering positions: ", a.count('\n')

# In[12]:


a = subprocess.check_output("grep -v '^#' " + vcf_bedfil + " | cut -f 1 | sort | uniq -c", shell=True)
print "Number of chrome: \n", a.count('\n')
print "Number of variants(Snps) in each chrome: \n", a

# In[13]:


# delete ##contig lines created by position filter. If not, it's false when running exomiser

ffinput = fileinput.input(vcf_bedfil, inplace=1)

for i, line in enumerate(ffinput):
    if line[0:8] == "##contig":
        newline = ""
    else:
        newline = line

    sys.stdout.write(newline)

ffinput.close()

# # ===============   3/ Split sample, Extract sample from multi-sample vcf =====================
#     To run exomiser, we have to split vcf file to be single-sample vcf

# In[14]:
print "================================================================================"
print "=========== 3/ Split sample, Extract sample from multi-sample vcf =============="
print "================================================================================"
# use vcftools to extract sample so it must have bgzip and tabix files\n\n\
print "***** Sort, zip, tabix file: ", vcf_bedfil
out_sort = vcf_bedfil[0:len(vcf_bedfil)-4] + "_sort.vcf"
bzip = out_sort + ".gz"
comm1 = "vcf-sort " + vcf_bedfil + " > " + out_sort
comm2 = "bgzip -c " + out_sort + " > " + bzip
comm3 = "tabix -p vcf " + bzip
subprocess.check_output(comm1, shell=True)
subprocess.check_output(comm2, shell=True)
subprocess.check_output(comm3, shell=True)

# In[15]:


print "***** Split samples ..."
for k in range(len(l_samples)):
    # monitor the process:
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d/%d" % ('=' * (1 + k * 20 / len(l_samples)), k+1, len(l_samples)))
    sys.stdout.flush()

    comm = "vcf-subset -c " + l_samples[k] + " " + bzip + " > " \
           + vcf_bedfil[0:len(vcf_bedfil)-4] + "_s" + str(k+1) + ".vcf"
    subprocess.check_output(comm, shell=True)
print "\n***** Done."

# # =================== 4/ Filter read depth (min_depth = 10, max_depth = 500)======================

# In[16]:
print "================================================================================"
print "========== 4/ Filter read depth (min_depth = 10, max_depth = 500)==============="
print "================================================================================"


def filter_DP(vcf_file):
    vcf_out = vcf_file[0:len(vcf_file)-4] + "_DP.vcf"
    f_in = open(vcf_file)
    f_out = open(vcf_out, "wb")

    read = f_in.readline()

    while read:
        new_line = read
        if "#" not in read[0]:
            row = new_line.split("\t")
            l_col = len(row)

            if l_col > 7:

                tmp_samp = row[l_col-1].split(":")
                if len(tmp_samp[2]) == 1:  # less than 10
                    new_line = ""
                else:
                    if int(tmp_samp[2]) > 500:
                        new_line = ""

        f_out.write(new_line)
        read = f_in.readline()

    f_in.close()
    f_out.close()


# In[17]:


print "***** Running filter read depth ..."
for k in range(len(l_samples)):
    file_vcf = vcf_bedfil[0:len(vcf_bedfil)-4] + "_s" + str(k+1) + ".vcf"
    filter_DP(file_vcf)
print "***** Done."

# In[18]:


print "***** File name after filter read depth: ", vcf_bedfil[0:len(vcf_bedfil)-4] + "_s" + str(1) + "_DP.vcf"

# # ======================= 5/ Annotation with Exomiser =============================

# In[19]:
print "================================================================================"
print "======================= 5/ Annotation by Exomiser ============================"
print "================================================================================"
print "***** Output for annotation in: ", out_exome
out_exome_conf = out_exome + "/conf"
out_exome_result = out_exome + "/result"
os.system("mkdir " + out_exome_conf)
os.system("mkdir " + out_exome_result)

# In[20]:


print "***** Generate conf files for running EXOMISER in ", out_exome_conf
for k in range(len(l_samples)):
    # print "sample:", k+1
    file_yml = out_exome_conf + "/exome_s" + str(k+1) + ".yml"
    os.system("cp " + input_exome + " " + file_yml)

    finput = fileinput.input(file_yml, inplace=1)
    file_vcf = vcf_bedfil[0:len(vcf_bedfil)-4] + "_s" + str(k+1) + "_DP.vcf"
    file_result = out_exome_result + '/' + name + "_exome_s" + str(k+1)

    for i, line in enumerate(finput):
        newline = line.replace("dataname", file_vcf)
        newline = newline.replace("resultname", file_result)
        sys.stdout.write(newline)

    finput.close()
print "***** Output this part will be in ", out_exome_result

# In[21]:


print "***** Annotation by EXOMISER *****"
"""
for k in range(len(l_samples)):
    #print "******************* Annotation sample " + str(k+1) + "/" + str(len(l_samples))
    # monitor the process:
    sys.stdout.write('\r') 
    sys.stdout.write("[%-20s] %d/%d" % ('='*(1+k*20/len(l_samples)), k+1,len(l_samples)))
    sys.stdout.flush()
    
    file_yml = out_exome_conf + "/exome_s" + str(k+1) + ".yml"
    comm = "java -Xms2g -Xmx10g -jar " + exom_path + " --analysis " + file_yml
    subprocess.check_output(comm, shell = True)
print "\n"   
"""

# In[22]:


# Parallel execute to speed up process of annotation EXOMISER
# EXOMISER also has paralle processing inside but if your computer has many threads, 
# Please try following code to speed up your process.

from joblib import Parallel, delayed

# print "***** Annotation by EXOMISER *****"


def exome_run(k):
    file_yml = out_exome_conf + "/exome_s" + str(k+1) + ".yml"
    comm = "java -Xms2g -Xmx10g -jar " + exom_path + " --analysis " + file_yml

    subprocess.check_output(comm, shell=True)
    print "Finish annotation sample", k+1


element_information = Parallel(n_jobs=-1)(delayed(exome_run)(k) for k in range(len(l_samples)))
print "***** Done step 5. \n"

# # ================  6/ Filter function_class (protein affecting) and maybe freq ===============

# In[22]:
print "================================================================================"
print "========  6/ Filter function_class (protein affecting) and maybe freq =========="
print "================================================================================"
# Remind list of protein affecting 
# af_pro = ["missense_variant", "splice_region_variant", "stop_gained", "frameshift_variant", "stop_lost", "inframe_deletion", "inframe_insertion"]
print "***** Filter genes with protein affecting: ", af_pro


def filter_fc(file_va, vcf_file):
    df_va = pd.read_csv(file_va, sep="\t")
    df_tmp = df_va[df_va.FUNCTIONAL_CLASS.isin(af_pro)]
    pos = list(df_tmp.POS)

    os.system("cp " + vcf_file + " " + vcf_file[0:len(vcf_file)-4] + "_fc.vcf")
    tmp = vcf_file[0:len(vcf_file)-4] + "_fc.vcf"
    tmpinput = fileinput.input(tmp, inplace=1)
    for i, line in enumerate(tmpinput):
        a = line.split("\t")

        if len(a) > 3:
            if a[1] == "POS":
                new_line = line
            else:
                if int(a[1]) in pos:
                    new_line = line
                else:
                    new_line = ""

        else:
            new_line = line
        sys.stdout.write(new_line)

        # sys.stderr.write(line +'\n')
    tmpinput.close()


# In[23]:


for k in range(len(l_samples)):
    file_va = out_exome_result + '/' + name + "_exome_s" + str(k+1) + ".variants.tsv"
    file_vcf = out_exome_result + '/' + name + "_exome_s" + str(k+1) + ".vcf"
    filter_fc(file_va, file_vcf)

print "***** Output file name: ", out_exome_result + '/' + name + "_exome_s" + str(1) + "_fc.vcf"

# # ================  7/ Merge results =============================================

# In[24]:
print "================================================================================"
print "================  7/ Merge results ============================================="
print "================================================================================"
# Merge vcf results
vcf_merge = out_exome + '/' + name + "_merge_exomiser_all.vcf"
# sort and bzip files
for k in range(len(l_samples)):
    result = out_exome_result + '/' + name + "_exome_s" + str(k+1) + ".vcf"
    out_sort = out_exome_result + '/' + name + "_exome_sort_s" + str(k+1) + ".vcf"
    bzip = out_exome_result + '/' + name + "_exome_bz_s" + str(k+1) + ".vcf.gz"

    comm1 = "vcf-sort " + result + " > " + out_sort
    comm2 = "bgzip -c " + out_sort + " > " + bzip
    comm3 = "tabix -p vcf " + bzip
    subprocess.check_output(comm1, shell=True)
    subprocess.check_output(comm2, shell=True)
    subprocess.check_output(comm3, shell=True)
comm4 = "vcf-merge"

for k in range(len(l_samples)):
    bzip = " " + out_exome_result + '/' + name + "_exome_bz_s" + str(k+1) + ".vcf.gz"
    comm4 = comm4 + bzip
comm4 = comm4 + " > " + vcf_merge
print subprocess.check_output(comm4, shell=True)

print "***** Merge results of exomiser. Output is ", vcf_merge

# In[25]:


# check number of variants after merging

a = subprocess.check_output("grep -v '^#' " + vcf_merge + " | cut -f 1 | sort", shell=True)
print "Number of variants after merging: ", a.count('\n')
print "Done step 7."

# # ==============  8/ Preprocessing and Burden test =============================
# 
#     Input to run rvtest-burden test:
#             - vcf file with genetype of all patients
#             - pheno file: phenotype of all patients
#             - genefile: mapping position with gene-names

# ## ==== 8.1/ Prepare input

# In[26]:
print "================================================================================"
print "==============  8/ Preprocessing and Burden test ==============================="
print "================================================================================"


def findchar(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


# In[28]:


# edit values, keep GT 1/1, 0/1 or 0/0
mfile = vcf_merge
print "***** Edit values, keep Genetype 1/1, 0/1 or 0/0"

os.system("cp " + mfile + " " + mfile[0:len(mfile)-4] + "_1.vcf")
mfile = mfile[0:len(mfile)-4] + "_1.vcf"
minput = fileinput.input(mfile, inplace=1)

for i, line in enumerate(minput):
    # sys.stderr.write(line +'\n')
    newline = line
    at1 = findchar(newline, "\t")
    if len(at1) > 5:
        line1 = newline[at1[8]+1:len(newline)]
        line1 = line1.replace("./", "0/")
        line1 = line1.replace("/.", "/0")
        line1 = line1.replace(".", "0/0")
        ac1 = findchar(line1, ":")
        while len(ac1) > 1:
            at2 = findchar(line1[ac1[0]:len(line1)], "\t")
            if len(at2) > 0:
                line1 = line1.replace(line1[ac1[0]:ac1[0]+at2[0]], "")
            else:
                line1 = line1.replace(line1[ac1[0]:len(line1)], "\n")

            ac1 = findchar(line1, ":")

        newline = newline[0:at1[8]+1] + line1
    sys.stdout.write(newline)
minput.close()
print "***** Output is ", mfile

# In[29]:


# bgzip and tabix files for burden test using rvtest
result = mfile
out_sort = mfile[0:len(mfile)-6] + "_sort.vcf"
bzip = out_sort[0:len(out_sort)-4] + "_bz.vcf.gz"
comm1 = "vcf-sort " + mfile + " > " + out_sort
comm2 = "bgzip -c " + out_sort + " > " + bzip
comm3 = "tabix -p vcf " + bzip
subprocess.check_output(comm1, shell=True)
subprocess.check_output(comm2, shell=True)
subprocess.check_output(comm3, shell=True)
print "***** Sort, zip, tabix file vcf to make input for burden test: ", bzip
print "***** Done step 8.1."

# ##  === 8.2/ Prepare phenotype values for burden test
# 
# extract samples ID in df_tab, then compare with ID 1 in df_bio to get values

# ### ===== View table and Check the correspondence table =====

# In[30]:


df_tab = pd.read_csv(input_table, sep="\t")

# In[31]:


df_tab.head()

# In[32]:


print "Verifying correspondence table ..."
flag = 0
for i in l_samples:
    if i not in list(df_tab["#HASH-ID"]):
        print i
        flag = flag + 1
if flag == 0:
    print "***** All samples are in the table."
else:
    print "***** Above samples are not in the table"


# In[33]:


def pheno_binary(idx1, idx2, phenoname):
    plq = []
    for i in range(len(df_tab)):
        tmp = df_tab["SAMPLES-NAMES"][i].split("-")
        if tmp[idx1][idx2] != 'n':
            plq.append(tmp[idx1][idx2])
        else:
            plq.append('0')
    df_tab[phenoname] = plq


# In[34]:


pheno_binary(3, 3, "PLq")
pheno_binary(4, 3, "HOX")
pheno_binary(5, 7, "HCaLciU")
pheno_binary(6, 6, "HCt3oL")
pheno_binary(7, 8, "RCt3oLOH")
print "***** Extract phenotype from sample names. Please check below table."


# In[35]:


def pheno_extract(idx1, phenoname):
    plq = []
    for i in range(len(df_tab)):
        tmp = df_tab["SAMPLES-NAMES"][i].split("-")
        plq.append(tmp[idx1])
    df_tab[phenoname] = plq


# In[36]:


pheno_extract(1, 'Sex')
pheno_extract(2, 'Age')
df_tab['Sex'] = df_tab['Sex'].map({'F': '2', 'M': '1'})

# In[37]:


df_tab.head()

# ### ===== Make phenotype file using df_tab =====

# In[38]:


# shell1 = "if [ ! -d '" +out_pheno+ "' ]; then\n\t" + "mkdir " + out_pheno + "\nfi"
# print subprocess.check_output(shell1, shell = True)
out_pheno = out_rvtest + "/pheno"
out_burden = out_rvtest + "/burden"
out_burden_gene = out_rvtest + "/burden_gene"

os.system("mkdir " + out_pheno)
os.system("mkdir " + out_burden_gene)
os.system("mkdir " + out_burden)

# In[39]:


print "***** Writing pheno file for burden test: ", out_pheno
f = open(out_pheno + '/pheno_all.ped', 'wb')
f.write("fid iid fatid matid sex age PLq HOX HCaLciU HCt3oL RCt3oLOH\n")

for k in range(len(df_tab)):
    name = df_tab["#HASH-ID"][k]

    text = name + " " + name + " 0 0 " + df_tab['Sex'][k] + ' ' + df_tab['Age'][k] + ' ' + df_tab['PLq'][k] \
           + ' ' + df_tab['HOX'][k] + ' ' + df_tab['HCaLciU'][k] + ' ' + df_tab['HCt3oL'][k] + ' ' \
           + df_tab['RCt3oLOH'][k] + "\n"
    f.write(text)

f.close()
print ('\nDone step 8.2')

# ## ==== 8.3/ Burden normal with genefile refFlat.txt.gz
# https://genome.sph.umich.edu/wiki/Rvtests

# ### 8.3.1/ Burden normal with only a genename 

# In[41]:


print "========= Burden normal with genefile refFlat.txt.gz"


def rvtest_gene(phenoname, gene):
    invcf = bzip
    pheno = out_pheno + '/pheno_all.ped'

    outrv = out_burden_gene + "/output_" + phenoname + "_" + gene

    comm = rvtest_pwd + " --inVcf " + invcf + " --pheno " + pheno + " --pheno-name " + phenoname \
           + " --covar " + pheno + " --covar-name age,sex --qtl --freqUpper 0.05 " + " --geneFile " \
           + genefile + " --gene " + gene + " --out " + outrv + " --kernel skat"
    print subprocess.check_output(comm, shell=True)


# ### ================== Modify this part for running burdentest function ===============
#     - rvtest_gene(phenoname, gene)
#     - phenoname: PLq, HOX, HCaLciU, HCt3oL, RCt3oLOH
#     - gene: name of gene

# In[42]:

gene_name = "ABCC6"
print "***** Burden test with some genes"
print "***** Output : " + out_burden_gene
print "----- For example: " + gene_name
rvtest_gene("PLq", gene_name)


# ### =====================================================================

# ### 8.3.2/ Burden test with all available genes

# In[46]:


def rvtest(phenoname):
    invcf = bzip
    pheno = out_pheno + '/pheno_all.ped'

    outrv = out_burden + "/output_" + phenoname

    comm = rvtest_pwd + " --inVcf " + invcf + " --pheno " + pheno + " --pheno-name " + phenoname + " --covar " \
           + pheno + " --covar-name age,sex --qtl --freqUpper 0.05" + " --geneFile " + genefile + " --out " \
           + outrv + " --kernel skat"  # adjust covar
    # print comm
    print subprocess.check_output(comm, shell=True)
    print "Please check your output in " + out_burden


# In[49]:


print "***** Output : ", out_burden
print "***** Running burden test ...\n"
rvtest("PLq")
rvtest("HOX")
rvtest("HCaLciU")
rvtest("HCt3oL")
rvtest("RCt3oLOH")
print "***** Done. Burden test completed!"

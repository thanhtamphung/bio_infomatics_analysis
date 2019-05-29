# coding: utf-8

# # MToolBox
#     Date: 8 Mar 2018
#     Last edit: 21 May 2019
#     Ref: https://github.com/mitoNGS/MToolBox/tree/master/MToolBox 
#     Please have a look the ref link above to understand parameters and more options
#     Runtime: 1H15' with input bam 8.2GB 

# ### Update
#     Edit bash file mtoolbox


import subprocess

#===================== PARAMETERS ===========================================
# path of the tool MToolBox folder which contains file MToolbox.sh
mtb_path = "/data/MToolBox/MToolBox-master/MToolBox/"
# the folder where you save your input file (fastq, bam ...)
input_path = "/data/MToolBox/tenon/S5/"
# names of inputs
f1 = "2052690574-005061_8008990406_S5.bam"
f2 = ""
# input type: fastq, bam, sam
input_type = "bam"
# please read doc MToolBox at ref-link to change options when running MToolBox
# for example: option = '-a "-z 0.8"' or option = ''
option = ''
#=============================================================================



a = f1.split('.')
name = a[0]
output_path = input_path + 'OUT_' + name + '/'

print ("***** Configure file and log file will be saved with names: %s.conf, %s.log in %s \n" %(name, name, input_path))
print ("***** The output vcf file will be saved with names: %s.vcf in %s \n" %(name, output_path))

conf_file = input_path + name + ".conf"
log_file = input_path + name + ".log"
list_file = input_path + name + ".lst"
mtb_bash = input_path + "mtoolbox"

#=============================================================================
## List file
f = open(list_file, 'wb')
f.write(f1 + '\n')
f.write(f2)
f.close()
print ("***** Generate list_file with all names of inputs for running MToolBox: " + list_file)

#=============================================================================
## Configure file
f = open(conf_file, 'wb')
f.write('#!/bin/bash\n')
content_conf = """
mtdb_fasta=chrM.fa
hg19_fasta=hg19RCRS.fa
mtdb=chrM
humandb=hg19RCRS
input_path={input_path}
output_name={output_path}
list={name}.lst
input_type={input_type}
ref=RCRS
vcf_name={name}
""".format(
input_path=input_path,
name=name,
output_path=output_path,
input_type=input_type
)
f.write(content_conf)
f.close()
print ("***** Generate conf_file for running MToolBox: " + conf_file)

#=============================================================================
## Bash/shell file to run MToolBox
f = open(mtb_bash, 'wb')
f.write('#!/bin/bash\n')
content_mtb = """
MTB={mtb_path}
INPUT={input_path}
NAME={name}
OPT={option}
export PATH="${{MTB}}:$PATH"
bash ${{MTB}}MToolBox.sh -i ${{INPUT}}${{NAME}}.conf ${{OPT}} &> ${{INPUT}}${{NAME}}.log
if [ $? -eq 0 ]; then 
    echo "Done! All output files are in ${{INPUT}}OUT_${{NAME}} ."
else
    echo "Not complete! Please check log file to figure out the problem."
fi
""".format(
mtb_path=mtb_path,
input_path=input_path,
name=name,
option=option
)
f.write(content_mtb)
f.close()
print ("***** Writing bash file to run MToolBox: " + mtb_bash)

#=============================================================================
## Run MToolBox
print ("***** Please view log file " + log_file + "\n ... still running ...") 
cmd = "bash " + mtb_bash
print subprocess.check_output(cmd, shell=True)




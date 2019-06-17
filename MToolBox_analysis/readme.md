# MToolBox
A pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing.

1/ Install MToolBox: https://github.com/mitoNGS/MToolBox/wiki/Installation 
  
2/ Download data here (in case you want a sample test): https://mega.nz/#!JTYAiIyJ 
          
          - 2052690574-005061_8008990406_S5.bam
          - 2052690574-005061_8008990406_S5.bam.bai
  
3/ Uncompress data and put it in the input folder
  
4/ Specify parameters in main.ipynb OR main_mtb.py

  - `mtb_path =`  : the path of the tool MToolBox folder which contains file MToolbox.sh

  - `input_path =`: the folder where you save your input file (fastq, bam ...)

  - `f1 =` , `f2 = `       : file input names, 
  
   for example `f1 = "bamfile.bam"`, `f2 = ""`, 
   if input files are fastq files, then `f1 = "myfile.R1.fastq.gz"`, `f2 = "myfile.R2.fastq.gz"`. 
   
   "FASTQ files MUST be renamed as `<sample_name>.R1.fastq`, `<sample_name>.R2.fastq` for PAIRED-END data and `<sample_name>.fastq` for SINGLE END data. FASTQ compressed input files could be accepted with *.fastq.gz extension."

  - `input_type =`: type of input file (fastq, bam, sam), for example `input_type = "bam"`
   
 5/ Run main.ipynb on Jupyter notebook
    OR run main_mtb.py on terminal: `time python main_mtb.py`
    
### Output
All output files are in the folder OUT_sample_name

- `../OUT_sample_name/sample_name.vcf`: prioritization analysis result

- `../OUT_sample_name/OUT_sample_name/sample_name.haplogroup.annotation.csv`: annotation result


### Processing

- Set up parameters

- Generate list file with all names of inputs

- Generate configure file for running MToolBox

- Generate bash file to run MToolBox

- Call bash file - run MToolBox 


For questions, please send to Phung Thanh Tam <thanhtam.phung@math.cnrs.fr> 

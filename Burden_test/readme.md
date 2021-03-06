# Annotation and Burden test

## Dependencies:

  1/ python packages in the import cell (pandas, numpy, pysam) in file .ipynb
  
  2/ bedtools
  
  3/ vcftools (vcf-subset, vcf-sort, vcf-merge)
  
  4/ htslib (bgzip, tabix)
  
  5/ exomiser 8.0.0 https://github.com/exomiser/Exomiser/releases/tag/8.0.0
  
  6/ rvtests https://github.com/zhanxw/rvtests

## Parameters
- input_path =: folder containing all input: vcf, bedfile, correspondance_table, exome file
- vcf_name =: vcf input file
- bedfile =: bed file
- ctable =: table link sample name id and phenotype
- exomefile =: model exomiser configure file
- genefile =: gene names with coodinates
- af_pro =: List of protein affections for filtering
- exom_path =: full path of "exomiser-cli-8.0.0.jar"
- rvtest_pwd =: full path of rvtests

## Burden test Framework

1/ Review data

2/ Filter position with bed file

3/ Split sample from vcf

4/ Filter read depth (min_base_count = 10) and maybe ref, alt

5/ Annotation with Exomiser

6/ Filter function_class (protein affecting) and maybe freq

7/ Merge results

8/ Preprocessing and Burden test

9/ Output report

<p align="center">
<img src="https://github.com/thanhtamphung/bio_infomatics_analysis/blob/master/Burden_test/process1.jpg" width="750">
</p>

## Run Burden test

1/ Install all dependencies

2/ Download data here: https://mega.nz/#!ZXJkXSYR

      - Lithiasis_set_20082018-123-samples.vcf

3/ Uncompress data and put it into the input folder (Lithiasis)

4/ Specify parameters and run burdentest_24apr.py OR burdentest_24apr.ipynb

## Output

- N_INFORMATIVE: Number of samples that are analyzed for association.
- NumVar: number of variants in the gene(or site)
- NumPolyVar: number of Polymorphic Genotypes
- Pvalue: permutation P-value

For questions, please send to Phung Thanh Tam <thanhtam.phung@math.cnrs.fr> 

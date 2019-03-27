# Annotation and Burden test

## Dependencies:

  1/ python packages in the import cell
  
  2/ bedtools
  
  3/ vcftools (vcf-subset, vcf-sort, vcf-merge)
  
  4/ htslib (bgzip, tabix)
  
  5/ exomiser 8.0.0 https://github.com/exomiser/Exomiser/releases/tag/8.0.0
  
  6/ rvtests https://github.com/zhanxw/rvtests

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

## Run Burden test

1/ Install all dependencies

2/ Download data here: https://mega.nz/#!ZXJkXSYR

      - Lithiasis_set_20082018-123-samples.vcf

3/ Uncompress data and put it in the input folder

4/ Specify parameters and run burdentest.ipynb

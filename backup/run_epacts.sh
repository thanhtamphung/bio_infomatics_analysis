#!/bin/bash

mkdir --p out/
NAME=annoexom
DIR=/usr/local/bin
VCF=/data/EPACTS/tamtest/${NAME}.vcf.gz
PED=/data/EPACTS/tamtest/pheno_all1.ped
OUT=/data/EPACTS/tamtest/out/test9
GRP=/data/EPACTS/tamtest/out/${NAME}.sites.anno.grp 


mkdir --p out/
echo "Performing EPACTS Binary Single Variant Score Test.."
${DIR}/epacts single --vcf ${VCF} --ped ${PED} --min-maf 0.001 --pheno PLq --test b.score --out ${OUT}.single.b.score --run 5 --anno --topzoom 1

echo "Creating Empirical Kinship matrix..."
${DIR}/epacts make-kin --vcf ${VCF} --ped ${PED} --min-maf 0.01 --out ${OUT}.single.q.emmax.kinf --run 5

echo "Running EMMAX single variant test..."
${DIR}/epacts single --vcf ${VCF} --ped ${PED} --min-maf 0.001 --pheno PLq --test q.emmax --out ${OUT}.single.q.emmax --kinf ${OUT}.single.q.emmax.kinf --run 5 --anno

echo "Taking the site list from the VCF file";
zcat ${VCF} | cut -f 1-8 | ${DIR}/bgzip -c > out/${NAME}.sites.vcf.gz

echo "Annotate the VCF file";
${DIR}/epacts anno --in out/${NAME}.sites.vcf.gz --out out/${NAME}.sites.anno.vcf.gz

echo "Create the input group file for burden test";
${DIR}/epacts make-group --vcf out/${NAME}.sites.anno.vcf.gz --out ${GRP} -format epacts -nonsyn

echo "Performing SKAT-O test";
${DIR}/epacts group --vcf ${VCF} --ped ${PED} --groupf ${GRP} --pheno PLq --cov age --cov sex --test skat -skat-adjust --out ${OUT}.gene.skat -field GT --max-maf 0.05 --run 4 

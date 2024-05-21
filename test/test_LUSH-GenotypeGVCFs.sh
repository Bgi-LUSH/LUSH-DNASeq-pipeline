mkdir -p ./outdir
export LD_LIBRARY_PATH=../bin/LUSH_toolkit-GenotypeGVCFs:$LD_LIBRARY_PATH
../bin/LUSH_toolkit-GenotypeGVCFs/lush_genotypegvcfs -I lush.g.vcf.gz -O lush.vcf.gz  --stand-call-conf 10 -t 52

mkdir -p ./outdir
export LD_LIBRARY_PATH=../bin/LUSH_toolkit-GenotypeGVCFs:$LD_LIBRARY_PATH
../bin/LUSH_toolkit-GenotypeGVCFs/lush_genotypegvcfs ../example_data/NA12878.g.vcf.gz ./outdir/NA12878.vcf.gz 10

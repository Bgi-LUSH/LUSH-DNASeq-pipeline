mkdir -p ./outdir
export LD_LIBRARY_PATH=../bin/LUSH_toolkit-HC:$LD_LIBRARY_PATH
../bin/LUSH_toolkit-HC/lush_hc HaplotypeCaller \
        --pcr-indel-model NONE \
        -I ../example_data/NA12878.sort.dup.bam \
        -R ../example_data/ref/chrM.fa \
        -O ./outdir/NA12878.vcf.gz

echo ALL done!

mkdir -p ./outdir
export LD_LIBRARY_PATH=../bin/LUSH_toolkit-HC:$LD_LIBRARY_PATH
../bin/LUSH_toolkit-HC/lush_hc \
        --pcr-indel-model NONE \
        --native-active-region-threads 3 \
        --native-main-spend-threads 16 \
        -I ../example_data/NA12878.sort.dup.BQSR.bam \
        -R ../example_data/ref/chrM.fa \
        -O ./outdir/NA12878.vcf.gz

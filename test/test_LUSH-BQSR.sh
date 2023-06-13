mkdir -p ./outdir
export LD_LIBRARY_PATH=../bin/LUSH_toolkit-BQSR:$LD_LIBRARY_PATH
../bin/LUSH_toolkit-BQSR/lush_bqsr \
 --bam_path ../example_data/NA12878.sort.dup.bam \
 --out_dir ./outdir  \
 --plugin_path ../bin/LUSH_toolkit-BQSR/libbqsr.so \
 --producer_number 2 --worker_number 21 \
 --fasta ../example_data/ref/chrM.fa \
 --known_site ../example_data/test.KnowSizes.vcf \
 --writer_thread 5 --pr_one_bam 1

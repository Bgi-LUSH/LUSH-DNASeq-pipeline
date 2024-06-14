mkdir -p ./outdir/ ./outdir/tem
../bin/LUSH_toolkit-Aligner/lush_aligner filter4mem \
        -6 ./outdir/ \
        -n 0.1 -J 0.5 -l 12 -g 2 -b 2 -t 20 -M \
        -r ../example_data/ref/chrM.fa \
        -o ./outdir/NA12878.sort.dup.bam \
        -Z ./outdir/tem \
        -i ../example_data/lush.config

echo ALL done!

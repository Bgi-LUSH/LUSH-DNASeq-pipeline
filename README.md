# Fast and accurate DNASeq Variant Calling workflow composed of LUSH-toolkits.
The LUSH pipeline reconstructs analysis tools SOAPnuke, BWA and GATK using C/C++, and employs a new parallel computing architecture. Several redundant intermediate steps of reading and writing the same file are eliminated in the LUSH pipeline, which greatly avoids unnecessary I/O throughput and improves CPU utilization. LUSH pipeline provides far superior computational speed to GATK/GATKSpark while maintaining a high level of accuracy comparable to that of GATK.

## USAGE for LUSH toolkits

###  USAGE for LUSH_Aligner 
LUSH_Aligner integrates four functional modules: preprocessing, alignment to genome, sort alignments and mark duplicates, corresponding to the software: SOAPnuke, Bwa mem, Samtools sort (Danecek et al. 2021) and GATK-MarkDuplicates (Picard), respectively. LUSH_Aligner reimplemented the original software algorithm using C and C++ and optimized the algorithm architecture to improve CPU utilization while reduce IO usage.  

LUSH_Aligner first needs to construct the index for the reference genome using the following command:
```
./bin/LUSH_toolkit-Aligner/lush_aligner1_share FilterBwaSortDup FilterBwaSortDup.yaml
```

LUSH_Aligner uses the filter4mem command to perform filtering, alignment, sorting alignments and marking duplicates.  
Common parameters for lush_aligner1_share FilterBwaSortDup in the yaml file：
```
--adapter1/-2	STR		3' adapter sequence
--adapter2/-3	STR 	5' adapter sequence [only for PE reads]
-4   Merge multiple pairs of clean FASTQ into one pair
--cleanfqOut/-6	STR Output directory. Processed fq files and statistical results would be output to here
--mixMatch/-b	INT misMatch: the max mismatch number when match the adapter (default: [1])
--lowQual/-l	FLOAT lowQual: low quality threshold (default: [5])
--qualRate/-J 	FLOAT qualRate: low quality rate (default: [0.5])
--nRate/-n 	FLOAT nRate: N rate threshold (default: [0.05])
--qualSys/-g 	INT qualSys: quality system 1:illumina, 2:sanger (default: [sanger])
--cutData/-C 	INT max reads number: limit the clean reads number(Mb)
-t  INT  number of threads [1]
-r  STR Reference sequence file
-I  STR config file for fq1 file and fq2 file
-o  STR the Bam file path after MarkDuplicates
-Y  use soft clipping for supplementary alignments
-M  mark shorter split hits as secondary
-Z  temporary file directory
```

config file for fq1 file and fq2 file should like this (tab-split):
```
./example_data/NA12878_l01_1.fq.gz  NA12878_l01_1       @RG\tID:NA12878.1\tLB:LibA\tSM:NA12878\tPL:COMPLETE\tCN:BGI
./example_data/NA12878_l01_2.fq.gz  NA12878_l01_2
./example_data/NA12878_l02_1.fq.gz  NA12878_l02_1       @RG\tID:NA12878.2\tLB:LibA\tSM:NA12878\tPL:COMPLETE\tCN:BGI
./example_data/NA12878_l02_2.fq.gz  NA12878_l02_2
```

Example:
```BASH
mkdir -p ./outdir/ ./outdir/tem
./bin/LUSH_toolkit-Aligner/lush_aligner1_share FilterBwaSortDup FilterBwaSortDup.yaml

```
We provide a quick run script in the test/ directory. You can directly navigate to the test directory by typing "cd ./test/" and then execute "sh ./test_LUSH-Aligner.sh".

###  USAGE for LUSH_BQSR
LUSH_BQSR is a C/C++ re-implementation of the GATK BaseRecalibrator and ApplyBQSR.  

Common parameters：
```
--bam_path	STR the input BAM file
--fasta	STR	Reference sequence file
--out_dir	STR	The output file to create
--plugin_path	STR	path for  libbqsr.so
--bed	STR	BED file for genomic intervals over which to operate	
--producer_number	INT	number of producer threads
--worker_number	INT	number of worker threads
--known_site	STR	One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis (should be decompression file).
--pr_one_bam	INT	When there are multiple producers, change the default split BAM output to merge BAM output
--writer_thread	INT	The number of threads used for writing
--ip	INT	Amount of padding (in bp) to add to each interval you are including
--bam_index	If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file
--quantize_qual	INT	Use static quantized quality scores to a given number of levels (with -bqsr)
```

Example:
```BASH
mkdir -p ./outdir
export LD_LIBRARY_PATH=./bin/LUSH_toolkit-BQSR:$LD_LIBRARY_PATH
./bin/LUSH_toolkit-BQSR/lush_bqsr \
 --bam_path /INPUT_PATH/NA12878.sort.dup.bam \
 --out_dir ./outdir/LUSH_BQSR  \
 --plugin_path ./bin/LUSH_toolkit-BQSR/libbqsr.so \
 --producer_number 2 --worker_number 21 \
 --fasta hg19.fa \
 --known_site Mills_and_1000G_gold_standard.indels.hg19.vcf \
 --writer_thread 5 --pr_one_bam 1
```
We provide a quick run script in the test/ directory. You can directly navigate to the test directory by typing "cd ./test/" and then execute "sh ./test_LUSH-BQSR.sh".

###  USAGE for LUSH_HC
LUSH_HC is a C/C++ re-implementation of the GATK HaplotypeCaller.

Common parameters：
```
Required:
  -I, --input <file1> <file2> ...             input file paths
  -O, --output <file>                         output file path
  -R, --reference <file>                      reference file path
Options:
  -H, --help                                  display help message
  -V, --version                               display version message
  -L, --interval <file>                       interval file path
  -P, --interval-padding <int>                interval padding size, must be a non-negative integer (default: 0)
  -Q, --base-quality-score-threshold <int>    base qualities threshold, must be in [6, 127](default: 18)
  -D, --max-reads-depth <int>                 maximum reads depth per alignment start position, must be a non-negative integer,set to 0 to disable(default: 50)
  -G, --gvcf-gq-bands <int1> <int2>...        specify GQ bands for merging non-variant sites in GVCF mode, must be in [1, 100], (default: 1, 2, 3,... 60, 70, 80, 90, 99)
      --nthreads <int>                        number of threads to use, must be in [1, 128] (default: 30)
      --pcr-indel-model <str>                 PCR indel model, available options: {NONE, HOSTILE, CONSERVATIVE, AGGRESSIVE} (default: CONSERVATIVE) 
      --emit-ref-confidence <str>             emit reference confidence score mode (default: NONE) available options: {NONE, BP_RESOLUTION, GVCF}
      --nstreampool <int> iostream pool size, must be in [1, 20] (default: 10)
      --inspect-reads strictly inspect input reads (default: false)
      --compression-level compression level, must be in [0, 9] (default: 6)

```

Example:
```BASH
mkdir -p ./outdir
export LD_LIBRARY_PATH=./bin/LUSH_toolkit-HC:$LD_LIBRARY_PATH
./bin/LUSH_toolkit-HC/lush_hc HaplotypeCaller \
        --pcr-indel-model NONE \
        -I /INPUT_PATH/NA12878.sort.dup.bam \
        -R hg19.fa \
        -O ./outdir/NA12878.vcf.gz
```
We provide a quick run script in the test/ directory. You can directly navigate to the test directory by typing "cd ./test/" and then execute "sh ./test_LUSH-HC.sh".

###  USAGE for LUSH_GenotypeGVCFs
LUSH_GenotypeGVCFs is a C/C++ re-implementation of the GATK GenotypeGVCFs.

**UASGE**: `LUSH_GenotypeGVCF inputGvcfFile  outputVcfFile  stand-call-conf`
```
inputGvcfFile   input VCF file
outputGvcfFile  output file name:/file/NA12878_PCR.vcf.gz
stand-call-conf The minimum phred-scaled confidence threshold at which variants should be called:10.0
```

Example:
```
mkdir -p ./outdir
export LD_LIBRARY_PATH=./bin/LUSH_toolkit-GenotypeGVCFs:$LD_LIBRARY_PATH
./bin/LUSH_toolkit-GenotypeGVCFs/lush_genotypegvcfs -I lush.g.vcf.gz  -O lush.vcf.gz  --stand-call-conf 10 -t 52
```
We provide a quick run script in the test/ directory. You can directly navigate to the test directory by typing "cd ./test/" and then execute "sh ./test_LUSH-GenotypeGVCFs.sh".

**If you encounter environmental issues when running the Lush toolkit on your system, you can try using the CentOS base image for execution, as Lush toolkit runs well in this environment. For example: `docker pull centos:centos7.6.1810`**

## Run the LUSH/GATK pipeline
###  Run the LUSH pipeline
```
Usage:
  LUSH_pipeline.sh [-i FQCONF] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX]
Description:
    FQCONF, the path of INPUT fastq file
    THREAD, the number of thread [10]
    OUTDIR, the path of outdir [./]
    PREFIX, the prefix of outputfile [LUSHtest]
    MODE, GVCF or not [Y/N]
```

Example:
```BASH
LUSH_pipeline.sh -i fq.config -t 40 -o ./  -m N -s samplename
```

fq.config should be like this:
```
/PATH/MGISEQ2000_PCR-free_NA12878_30X_1.fq.gz   NA12878_30X_1   @RG\tID:NA12878_30X.1\tLB:NA12878_30X\tSM:NA12878_30X\tPL:COMPLETE\tCN:BGI
/PATH/MGISEQ2000_PCR-free_NA12878_30X_2.fq.gz   NA12878_30X_2
```

### Run the GATK/GATKSpark pipeline

```
Usage:
  GATK_pipeline.sh [-i FQFile] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX] [-p SPARK]

Description:
    FQFile, the path of INPUT fastq file, should be like '/path/fastq1,/path/fatq2'
    THREAD, the number of thread [10]
    OUTDIR, the path of outdir [./]
    PREFIX, the prefix of outputfile [GATKtest]
    MODE, GVCF or not [Y/N]
    SPARK, Spark or not [Y/N]
```

Example:
```BASH
GATK_pipeline.sh -i /PATH/MGISEQ2000_PCR-free_NA12878_30X_1.fq.gz,/PATH/MGISEQ2000_PCR-free_NA12878_30X_2.fq.gz -t 40 -o ./  -m N -s samplename -p N

```

### Accuracy Validation
```
Usage:
  Haplotype_Comparison.sh [-i VCFFile] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX]
Description:
    VCFFile, the path of INPUT vcf file, should be like '/path/*1.vcf.gz,/path/*2.vcf.gz'
    THREAD, the number of thread [10]
    OUTDIR, the path of outdir [./]
    PREFIX, the prefix of outputfile [HAPtest]

```
Example:
```BASH
Haplotype_Comparison.sh -i LUSHtest.vcf.gz,GATKtest.vcf.gz -t 40 -o ./ -s sample
```

**Note that only sample scripts are provided here. Actual operation requires replacing the paths of the corresponding software and configuration files as required.**
#!/bin/bash

usage() {
    echo "Usage:"
    echo "  $0 [-i FQFile] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX] [-p SPARK]"
    echo "Description:"
    echo "    FQFile, the path of INPUT fastq file, should be like '/path/fastq1,/path/fatq2'"
    echo "    THREAD, the number of thread [10]"
    echo "    OUTDIR, the path of outdir [./]"
    echo "    PREFIX, the prefix of outputfile [GATKtest]"
	echo "    MODE, GVCF or not [Y/N]"
    echo "    SPARK, Spark or not [Y/N]"
    exit 1
}

while getopts i:t:o:s:m:p:h OPT; do
    case $OPT in
		i) FQFile="$OPTARG";;
		t) THREAD="$OPTARG";;
		o) OUTDIR="$OPTARG";;
		s) PREFIX="$OPTARG";;
		m) MODE="$OPTARG";;
        p) SPARK="$OPTARG";;
		h) usage;;
		?) usage;;
    esac
done; if [ ! $FQFile ]; then  usage;fi

infq1=`echo $FQFile|cut -f 1 -d ','`
infq2=`echo $FQFile|cut -f 2 -d ','`
outdir=${OUTDIR:-'./'}
samname=${PREFIX:-'GATKtest'}
thread=${THREAD:-'10'}

if [ ! -d "$outdir" ]  
then  
    mkdir -p "$outdir"  
    echo "Directory created: $outdir"  
else  
    echo "Directory already exists: $outdir"  
fi  

logfile=$outdir/LUSH_pipeline_$(date +"%Y%m%d%H%M").log
echo -e "Input configuration:\n\tFQ1: $infq1\n\tFQ2: $infq2\n\tOUTDIR: $outdir\n\tPREFIX: $samname\n\tTHREAD: $thread\nSee more detail in $logfile \n"
# Log file
exec >$logfile 2>&1

#Some parameters to be provided by the user
#1. software  
bin=/USER_HOME  # Fill it according to the user's real directory
SOAPnuke=$bin/SOAPnuke-SOAPnuke2.1.5/SOAPnuke
java=$bin/java8/jdk1.8.0_60/bin/java
bwa=$bin/bwa
gatk=$bin/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
samtools=$bin/samtools
#2. Reference files
db=/USER_HOME  #Fill it according to the user's real directory
fasta=$db/hg19/hg19.fasta
db_bqsr_mills=$db/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz
db_bqsr_1000g=$db/hg19/gatk/1000G_phase1.indels.hg19.vcf.gz
db_bqsr_dbsnp=$db/hg19/gatk/dbsnp_151.hg19.vcf.gz


# Some functions
check_path() 
{
	for arg in "$@"
	do
		if [ ! -f $arg ]
		then
			echo "ERROR: $arg not exists!!!"
			exit 1
		fi
	done
}

check_error() 
{
    if [ "$1" = "0" ]
    then
        echo "$2 successfully."
        return
    fi
    echo $2 failed.
    exit $1
}

delta_time()
{
    start=$1
    end=$2
    dt=$((end-start))
    dh=`echo $dt|awk '{print int($1/3600)}'`
    dt2=`echo "$dt $dh" |awk '{print int($1-3600*$2)}'`
    dm=`echo "$dt2"|awk '{print int($1/60)}' `
    ds=`echo "$dt $dm" |awk '{print int($1-60*$2)}'`
    echo "$dh:$dm:$ds"
}

run()
{
    cmd=$1
    step=$2
    echo "Starting to run $step: $cmd"
    start=`date +"%D %T"`
    start_s=`date +%s`
    echo "$step start time: $start"
    eval "$cmd"
    check_error $? "$step"
    end=`date +"%D %T"`
    end_s=`date +%s`
    echo "$step end time: $end"
    runtime=$(delta_time $start_s $end_s)
    echo "$step runtime: $runtime"
}

#check input file
check_path $SOAPnuke $java $bwa $gatk $samtools $fasta $db_bqsr_mills $db_bqsr_1000g $db_bqsr_dbsnp


### SOAPnuke 
mkdir -p $outdir
mkdir -p $outdir/filter
#cut='-w 500000'
para="-n 0.1 -q 0.5 -T $thread -l 12 "
filter_cmd="$SOAPnuke filter \
  $para \
 -1 $infq1 -2 $infq2 \
 -C $samname.clean_1.fq.gz -D $samname.clean_2.fq.gz \
 -o $outdir/filter"
run "$filter_cmd" SOAPnuke_filter

### bwa mem alignment
mkdir -p  $outdir/bwaMem
bwa_cmd="$bwa mem \
 -t $thread -M -Y \
 -R \"@RG\tID:$samname\tLB:$samname\tSM:$samname\tPL:COMPLETE\tCN:BGI\" \
 $fasta \
 $outdir/filter/$samname.clean_1.fq.gz  $outdir/filter/$samname.clean_2.fq.gz | $samtools view -Sb \
 -o  $outdir/bwaMem/$samname.raw.bam -"
run "$bwa_cmd" BWA-alignment

### samtools
sort_cmd="$samtools sort -@ $thread $outdir/bwaMem/$samname.raw.bam -o $outdir/bwaMem/$samname.sorted.bam"
run "$sort_cmd" Samtools_Sort


if [ "${SPARK}" = "Y" ] ;then
    MARKDUP="MarkDuplicatesSpark --conf spark.local.dir=$outdir/MarkDuplicates/tmp"
    BQSR="BaseRecalibratorSpark --conf spark.local.dir=$outdir/BQSR/tmp"
    APPBQSR="ApplyBQSRSpark --conf spark.local.dir=$outdir/BQSR/tmp"
	HC="HaplotypeCallerSpark --conf spark.local.dir=$outdir/HC/tmp"
else
    MARKDUP="MarkDuplicates"
    BQSR="BaseRecalibrator"
    APPBQSR="ApplyBQSR"
	HC="HaplotypeCaller"
fi

### MarkDuplicates
mkdir -p $outdir/MarkDuplicates
mkdir -p $outdir/MarkDuplicates/tmp
markdup_cmd1="$java -Xmx6g -XX:ParallelGCThreads=4 -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $gatk $MARKDUP -I $outdir/bwaMem/$samname.sorted.bam  \
 -M $outdir/MarkDuplicates/$samname.sorted.bam.mat \
 -O $outdir/MarkDuplicates/$samname.sorted.rmdup.bam"
markdup_cmd2="$samtools index -@ $thread $outdir/MarkDuplicates/$samname.sorted.rmdup.bam"
run "$markdup_cmd1 && $markdup_cmd2" GATK_MarkDuplicates

###BQSR
mkdir -p $outdir/BQSR
mkdir -p $outdir/BQSR/tmp
bqsr_cmd1="$java -Xmx6g -XX:ParallelGCThreads=4 -jar $gatk $BQSR \
 -R $fasta  \
 -I $outdir/MarkDuplicates/$samname.sorted.rmdup.bam \
 --known-sites $db_bqsr_mills \
 --known-sites $db_bqsr_1000g  \
 --known-sites $db_bqsr_dbsnp \
 -O $outdir/BQSR/$samname.recal.table"

bqsr_cmd2="$java -Xmx6g -XX:ParallelGCThreads=4 -jar $gatk $APPBQSR \
 -R $fasta  \
 -I $outdir/MarkDuplicates/$samname.sorted.rmdup.bam --bqsr-recal-file $outdir/BQSR/$samname.recal.table \
 -O $outdir/BQSR/$samname.sorted.rmdup.bqsr.bam"
run "$bqsr_cmd1 && $bqsr_cmd2" GATK_BQSR


### HaplotypeCaller
mkdir -p $outdir/HC
mkdir -p $outdir/HC/tmp

hc_cmd="$java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xss4m -jar $gatk $HC \
 --pcr-indel-model NONE \
 -R $fasta \
 -I $outdir/BQSR/$samname.sorted.rmdup.bqsr.bam\
 -O $outdir/HC/$samname.g.vcf.gz"
hcg_cmd="$java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xss4m -jar $gatk $HC \
 -ERC GVCF \
 --pcr-indel-model NONE \
 -R $fasta \
 -I $outdir/BQSR/$samname.sorted.rmdup.bqsr.bam\
 -O $outdir/HC/$samname.g.vcf.gz"
gvcf_cmd="$java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -jar $gatk GenotypeGVCFs \
 -R $fasta \
 --variant $outdir/HC/$samname.g.vcf.gz \
 -O $outdir/HC/$samname.vcf.gz \
 -stand-call-conf 10"

if [ "${MODE}" = "Y" ] ;then
	# echo "$hcg_cmd" LUSH_HC-GVCF-mode
	# echo "$gvcf_cmd" LUSH_GenotypeGVCFs
	run "$hcg_cmd" GATK_HC-GVCF-mode
	run "$gvcf_cmd" GATK_GenotypeGVCFs
else
	# echo "$hc_cmd" LUSH_HC
	run "$hc_cmd" GATK_HC
fi

#!/bin/bash

usage() {
    echo "Usage:"
    echo "  $0 [-i FQCONF] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX]"
    echo "Description:"
    echo "    FQCONF, the path of INPUT fastq file"
    echo "    THREAD, the number of thread [10]"
    echo "    OUTDIR, the path of outdir [./]"
    echo "    PREFIX, the prefix of outputfile [LUSHtest]"
	echo "    use GVCF mode or not [Y/N]"
    exit 1
}

while getopts i:t:o:s:m:h OPT; do
    case $OPT in
		i) FQCONF="$OPTARG";;
		t) THREAD="$OPTARG";;
		o) OUTDIR="$OPTARG";;
		s) PREFIX="$OPTARG";;
		m) MODE="$OPTARG";;
		h) usage;;
		?) usage;;
    esac
done; if [ ! $FQCONF ]; then  usage;fi

fqconfig=$FQCONF
outdir=${OUTDIR:-'./outdir_for_LUSH_pipeline'}
samname=${PREFIX:-'LUSHtest'}
thread=${THREAD:-'10'}

if [ ! -d "$outdir" ]  
then  
    mkdir -p "$outdir"  
    echo "Directory created: $outdir"  
else  
    echo "Directory already exists: $outdir"  
fi  

logfile=$outdir/LUSH_pipeline_$(date +"%Y%m%d%H%M").log
echo -e "Input configuration:\n\tFQCONF: $fqconfig\n\tOUTDIR: $outdir\n\tPREFIX: $samname\n\tTHREAD: $thread\nSee more detail in $logfile \n"
# Log file
exec >$logfile 2>&1

#Some parameters to be provided by the user
#1. software  
bin=../bin  # Fill it according to the user's real directory
lush_aligner=$bin/LUSH_toolkit-Aligner/lush_aligner
lush_sdk=$bin/LUSH_toolkit-BQSR/lush_bqsr
lush_hc=$bin/LUSH_toolkit-HC/lush_hc
LUSH_GenotypeGVCFs=$bin/LUSH_toolkit-GenotypeGVCFs/lush_genotypegvcfs
samtools=$bin/samtools
#2. Reference files
db=/USER_HOME  #Fill it according to the user's real directory
fasta=../example_data/ref/chrM.fa  # use example  ref
db_knowsize=../example_data/test.KnowSizes.vcf

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
check_path $fqconfig $lush_aligner $lush_sdk $lush_hc $LUSH_GenotypeGVCFs $fasta $db_bqsr_mills $db_bqsr_1000g $db_bqsr_dbsnp $samtools


###########RUN LUSH pipeline

###Lush_aligner
mkdir -p $outdir
mkdir -p $outdir/Lush_aligner
mkdir -p $outdir/tem
#cut='-C 0.01'
aligner_opts="-n 0.1 -J 0.5 -l 12 -g 2 -b 2 -M -Y"
align_cmd1="$lush_aligner filter4mem \
	-6 $outdir/Lush_aligner/ \
	$aligner_opts \
	-t $thread\
	-r $fasta \
	-o $outdir/Lush_aligner/$samname.sort.dup.bam \
	-Z $outdir/tem \
	-i $fqconfig"
align_cmd2="$samtools index -@ $thread  $outdir/Lush_aligner/$samname.sort.dup.bam"
#echo "$align_cmd1 && $align_cmd2" LUSH_Aligner
run "$align_cmd1 && $align_cmd2" LUSH_Aligner

### LUSH BQSR
mkdir -p $outdir/LUSH_BQSR
export LD_LIBRARY_PATH=$bin/LUSH_toolkit-BQSR:$LD_LIBRARY_PATH
bqsr_cmd="$lush_sdk \
 --bam_path $outdir/Lush_aligner/$samname.sort.dup.bam \
 --out_dir $outdir/LUSH_BQSR  \
 --plugin_path $bin/LUSH_toolkit-BQSR/libbqsr.so \
 --producer_number 2 --worker_number 21 \
 --fasta $fasta \
 --known_site $db_knowsize \
 --writer_thread 5 --pr_one_bam 1"
# echo "$bqsr_cmd" LUSH_BQSR
run "$bqsr_cmd" LUSH_BQSR

### LUSH HC
mkdir -p $outdir/LUSH_HC

hc_cmd="$lush_hc HaplotypeCaller \
	--pcr-indel-model NONE \
	--nthreads $thread \
	-I $outdir/LUSH_BQSR/$samname.sort.dup.BQSR.bam \
	-R $fasta -O $outdir/LUSH_HC/$samname.g.vcf.gz"
hcg_cmd="$lush_hc HaplotypeCaller \
	--pcr-indel-model NONE \
	--emit-ref-confidence GVCF \
	--nthreads $thread \
	-I $outdir/LUSH_BQSR/$samname.sort.dup.BQSR.bam \
	-R $fasta -O $outdir/LUSH_HC/$samname.g.vcf.gz"

gvcf_cmd="$LUSH_GenotypeGVCFs $outdir/LUSH_HC/$samname.g.vcf.gz $outdir/LUSH_HC/$samname.vcf.gz 10"

if [ "${MODE}" == "Y" ] ;then
	# echo "$hcg_cmd" LUSH_HC-GVCF-mode
	# echo "$gvcf_cmd" LUSH_GenotypeGVCFs
	export LD_LIBRARY_PATH=$bin/LUSH_toolkit-HC:$LD_LIBRARY_PATH
	run "$hcg_cmd" LUSH_HC-GVCF-mode
	export LD_LIBRARY_PATH=$bin/LUSH_toolkit-GenotypeGVCFs:$LD_LIBRARY_PATH
	run "$gvcf_cmd" LUSH_GenotypeGVCFs
else
	# echo "$hc_cmd" LUSH_HC
	export LD_LIBRARY_PATH=$bin/LUSH_toolkit-HC:$LD_LIBRARY_PATH
	run "$hc_cmd" LUSH_HC
fi


echo "all DONE!"

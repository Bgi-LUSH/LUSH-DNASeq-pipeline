#!/bin/bash

usage() {
    echo "Usage:"
    echo "  $0 [-i VCFFile] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX]"
    echo "Description:"
    echo "    VCFFile, the path of INPUT vcf file, should be like '/path/*1.vcf.gz,/path/*2.vcf.gz'"
    echo "    THREAD, the number of thread [10]"
    echo "    OUTDIR, the path of outdir [./]"
    echo "    PREFIX, the prefix of outputfile [HAPtest]"
    exit 1
}

while getopts i:t:o:s:h OPT; do
    case $OPT in
		i) VCFFile="$OPTARG";;
		t) THREAD="$OPTARG";;
		o) OUTDIR="$OPTARG";;
		s) PREFIX="$OPTARG";;
		h) usage;;
		?) usage;;
    esac
done; if [ ! $VCFFile ]; then  usage;fi

vcf1=`echo $VCFFile|cut -f 1 -d ','`
vcf2=`echo $VCFFile|cut -f 2 -d ','`
outdir=${OUTDIR:-'./'}
prefix=${PREFIX:-'HAPtest'}
thread=${THREAD:-'10'}
logfile=$outdir/_$(date +"%Y%m%d%H%M").log
echo -e "Input configuration:\n\tvcf1: $vcf1\n\tvcf2: $vcf2\n\tOUTDIR: $outdir\n\tPREFIX: $prefix\n\tTHREAD: $thread\nSee more detail in $logfile \n"
# Log file
exec >$logfile 2>&1

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

### The gold standard truth variant calls and high confidence genomic intervals 
dir=/USER_HOME  # Fill it according to the user's real directory
gold=$dir/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz #for NA12878
bed=$dir/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed #for NA12878
ref=$dir/db/hg19/hg19.fasta

#software
bin=/USER_HOME # Fill it according to the user's real directory
hap_py=$bin/hap.py

#check para path
check_path $vcf1  $vcf2 $hap_py $ref $bed 


mkdir -p $outdir
mkdir -p $outdir/tmp
export HGREF=$ref

echo "VCF1 vs TRUTH running ..."
# compare VCF1 vs TRUTH
python $hap_py $gold $vcf1 \
    -f  $bed \
    --preprocess-truth  --usefiltered-truth \
    --write-vcf \
    --scratch-prefix $outdir/tmp \
    --threads 10  -o $outdir/$prefix-VCF1vsTRUTH > $outdir/$prefix-VCF1vsTRUTH.log

# compare VCF2 vs TRUTH
echo "VCF2 vs TRUTH running ..."
python $hap_py $gold $vcf2 \
    -f  $bed \
    --preprocess-truth  --usefiltered-truth \
    --write-vcf \
    --scratch-prefix $outdir/tmp \
    --threads 10  -o $outdir/$prefix-VCF2vsTRUTH > $outdir/$prefix-VCF2vsTRUTH.log

# compare VCF1 vs VCF2
echo "VCF1 vs VCF2 running ..."
python $hap_py $vcf1 $vcf2 \
    --preprocess-truth  --usefiltered-truth \
    --write-vcf \
    --scratch-prefix $outdir/tmp \
    --threads 10  -o $outdir/$prefix-VCF1vsVCF2 > $outdir/$prefix-VCF1vsVCF2.log
# Fast and accurate DNASeq Variant Calling workflow composed of LUSH-toolkits.
The LUSH pipeline reconstructs analysis tools SOAPnuke, BWA and GATK using C/C++, and employs a new parallel computing architecture. Several redundant intermediate steps of reading and writing the same file are eliminated in the LUSH pipeline, which greatly avoids unnecessary I/O throughput and improves CPU utilization. LUSH pipeline provides far superior computational speed to GATK while maintaining a high level of accuracy comparable to that of GATK.

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
LUSH_pipeline.sh -i fq.config -t 40 -o ./  -m Y -s samplename
```

fq.config should be like this:
```
/PATH/MGISEQ2000_PCR-free_NA12878_30X_1.fq.gz   NA12878_30X_1   @RG\tID:NA12878_30X.1\tLB:NA12878_30X\tSM:NA12878_30X\tPL:COMPLETE\tCN:BGI
/PATH/MGISEQ2000_PCR-free_NA12878_30X_2.fq.gz   NA12878_30X_2

```

### Run the GATK pipeline

```
Usage:
  GATK_pipeline.sh [-i FQFile] [-t THREAD] [-o OUTDIR] [-m MODEL] [-s PREFIX]
Description:
    FQFile, the path of INPUT fastq file, should be like '/path/fastq1,/path/fatq2'
    THREAD, the number of thread [10]
    OUTDIR, the path of outdir [./]
    PREFIX, the prefix of outputfile [GATKtest]
    MODE, GVCF or not [Y/N]

```

Example:
```BASH
GATK_pipeline.sh -i /PATH/MGISEQ2000_PCR-free_NA12878_30X_1.fq.gz,/PATH/MGISEQ2000_PCR-free_NA12878_30X_2.fq.gz -t 40 -o ./  -m Y -s samplename
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


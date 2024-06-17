FROM centos:centos7.6.1810

COPY LUSH-DNASeq-pipeline/ /usr/LUSH-DNASeq-pipeline/

WORKDIR /usr/LUSH-DNASeq-pipeline/

ENV PATH="/usr/LUSH-DNASeq-pipeline/bin/LUSH_toolkit-Aligner/:/usr/LUSH-DNASeq-pipeline/bin/LUSH_toolkit-BQSR/:/usr/LUSH-DNASeq-pipeline/bin/LUSH_toolkit-GenotypeGVCFs/:/usr/LUSH-DNASeq-pipeline/bin/LUSH_toolkit-HC/:${PATH}"

CMD ["/bin/bash"]

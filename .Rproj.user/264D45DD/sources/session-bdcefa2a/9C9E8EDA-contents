#!/bin/bash

#code/trimSort.sh lib01

#input - folder in output/demux
lib=$1
libPath="output/demux/${lib}/"

#make folders for trimming output
outFungi="output/trim/${lib}/fungi/"
mkdir -p $outFungi

#make directory for scratch data
scratch="output/scratch/"
mkdir -p $scratch

#set gene primers for trimming
fwd_primer_fungi="TCCTCCGCTTATTGATATGC" #ITS4
rc_fwd_primer_fungi="GCATATCAATAAGCGGAGGA"
rev_primer_fungi="CAHCGATGAAGAACRYAG" #ITS3_kyo1
rc_rev_primer_fungi="CTRYGTTCTTCATCGDTG"

#Make file for trimming summary output
log="output/trim/${lib}.summary.tab"
echo -e "Sample\tRaw\tFungi.trim1\tFungi.trim2" > ${log}

#Loop over samples
for FWD in $(find $libPath -name "*R1.fastq.gz" | grep -v "undetermined"); do
    REV=`echo $FWD | sed 's/R1.fastq.gz/R2.fastq.gz/'`
    #current sample
    SAMP=`echo $FWD | awk -F'_' '{print $3 "_" $4}' | awk -F'.' '{print $1}'`  
    
    #############
    ### FUNGI ###
    #############

    #trim forward reads
    trim1FR1="${scratch}${SAMP}.f.R1.trim1.fq.gz"
    trim1FR2="${scratch}${SAMP}.f.R2.trim1.fq.gz"
    cutadapt --quiet -g $fwd_primer_fungi --discard-untrimmed --overlap 20 -e 0.17 --minimum-length 225 --maximum-length 228 -o $trim1FR1 -p $trim1FR2 $FWD $REV

    #trim reverse reads
    trim2FR1="${scratch}${SAMP}.f.R1.trim2.fq.gz"
    trim2FR2="${scratch}${SAMP}.f.R2.trim2.fq.gz"
    cutadapt --quiet -g $rev_primer_fungi --discard-untrimmed --overlap 18 -e 0.17 --minimum-length 227 --maximum-length 230 -o $trim2FR2 -p $trim2FR1 $trim1FR2 $trim1FR1

    #Count seqs after first trim step
    if [[ -f $trim2FR1 && -f $trim2FR2 ]]; then
        Fungi_mid=`gzip -cd $trim2FR1 | grep -c '^@M01498'`
    else
        Fungi_mid=0
    fi

    #trim overhang
    if [ $Fungi_mid -gt 0 ]; then
        trim3FR1="${scratch}${SAMP}.f.R1.trim3.fq.gz"
        trim3FR2="${scratch}${SAMP}.f.R2.trim3.fq.gz"
        SeqPurge -in1 $trim2FR1 -in2 $trim2FR2 -out1 $trim3FR1 -out2 $trim3FR2 -a1 $rc_rev_primer_fungi -a2 $rc_fwd_primer_fungi -qcut 0 -ncut 0 -min_len 100 -summary ${scratch}${SAMP}.f.log.txt
    fi

    #Count seqs after final step
    if [[ -f $trim3FR1 && -f $trim3FR2 ]]; then
        Fungi_final=`gzip -cd $trim3FR1 | grep -c '^@M01498'`
    else
        Fungi_final=0
    fi
    #Keep only samples that made it through the processing (also truncate reads)
    if [ $Fungi_final -gt 0 ]; then
        outFR1="${outFungi}${SAMP}.f.R1.fq.gz"
        outFR2="${outFungi}${SAMP}.f.R2.fq.gz"
        cutadapt --quiet -g XX -G XX --length 225 -o $outFR1 -p $outFR2 $trim3FR1 $trim3FR2
    fi
    
    #make summary and append to file
    RAW=`gzip -cd $FWD | grep -c '^@M01498'`
    echo -e $SAMP"\t"$RAW"\t"$Fungi_mid"\t"$Fungi_final| cat >> ${log}

    #remove tmp files
    rm ${scratch}${SAMP}.*
done

rm -r ${scratch}


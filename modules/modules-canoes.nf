#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref    = file(params.ref, type: 'file')
probes = file(params.probes, type: 'file')
outdir = file(params.outdir, type: 'dir')

process genReadCounts {
    tag { "read_counts" }
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    module 'bedtools/2.30.0'
    
    input:
    path(bam)

    output:
    tuple val("canoes_in"), path("canoes.reads.txt"), path("sample_list"), emit: canoes_reads
    
    """
    sort bam_list_unsorted.txt > bam_list_sorted.txt
    bedtools multicov -bams `cat bam_list_sorted.txt` -bed ${probes} -q 20 > canoes.reads.txt
    while read line
    do
        echo `basename \$line` | sed 's/.bam//'
    done < bam_list_sorted.txt > sample_list
    """
}

process calcGC_CANOES {
    tag {" calc_gc "}
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    output:
    path("gc.txt"), emit: gc_content 
    
    """
    java -Xmx2000m -Djava.io.tmpdir=TEMP -jar /home/phelelani/applications/gatk-2.1-9/GenomeAnalysisTK.jar \
        -T GCContentByInterval \
        -L ${probes} \
        -R ${ref} \
        -o gc.txt
    """
}

process runCANOES {
    tag { "run_canoes" }
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    input:
    tuple val(canoes_in), path(canoes_reads), path(sample_list)
    path(gc_content)

    output:
    path("*"), emit: canoes_out
    path("Sample_CNVs.csv"), emit: cnvs
    
    """
    sed 's/chr//g; /^X/d; /^Y/d' ${canoes_reads} > canoes.reads_new.txt
    sed 's/chr//g; /^X/d; /^Y/d' ${gc_content} > gc_new.txt 
    
    run_canoes.R gc_new.txt canoes.reads_new.txt ${sample_list}
    """
}

process filterCANOESCNVs {
    tag { 'filter_cnvs' }
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    input:
    path(cnvs)
    
    output:
    path("Sample_CNVs_filtered.csv"), emit: filtered_cnvs
    
    """
    awk '{ if( (\$10>=80) && (\$10!="NA") && (\$4>=100) ) { print } }' Sample_CNVs.csv > Sample_CNVs_filtered.csv
    """
}

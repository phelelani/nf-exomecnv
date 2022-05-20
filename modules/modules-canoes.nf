#!/usr/bin/env nextflow

ref    = file(params.ref, type: 'file')
probes = file(params.probes, type: 'file')
outdir = file(params.outdir, type: 'dir')

process genReadCounts {
    tag { "read_counts" }
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    module 'bedtools/2.30.0'
    
    input:
    path bams.collectFile () { item -> [ 'bamlist_unsorted.txt', "${item.get(1)[0]}" + '\n' ] }

    output:
    tuple val("canoes_in"), path("canoes.reads.txt"), path("sample_list"), emit canoes_reads
    
    """
    sort bamlist_unsorted.txt > bamlist.txt
    bedtools multicov -bams `cat bamlist.txt` -bed ${probes} -q 20 > canoes.reads.txt
    awk -F'/' '{ print \$8 }' bamlist.txt | sed 's/.bam//' > sample_list
    """
}

process calcGC {
    tag {" calc_gc "}
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    output:
    path("gc.txt"), emit gc_content 
    
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
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    input:
    tuple val(canoes_in), path(canoes_reads), path(sample_list)
    path(gc_content)

    output:
    path("*"), emit: canoes_out
    
    """
    sed 's/chr//g; s/^X/d; s/^Y/d' ${canoes_reads} > canoes.reads_new.txt
    sed 's/chr//g; s/^X/d; s/^Y/d' ${gc_content} > gc_new.txt 
    
    run_canoes.R gc_new.txt canoes.reads_new.txt ${sample_list}
    """
}

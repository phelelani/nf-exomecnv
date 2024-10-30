#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref    = file(params.ref, type: 'file')
probes = file(params.probes, type: 'file')
outdir = file(params.outdir, type: 'dir')

process calcGC_CANOES {
    tag { "${chr}" }
    label 'gatk'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
    
    input:
    each chr
    
    output:
    tuple val("${chr}"), path("${chr}_gc.txt"), emit: chr_gc_content 

    script:
    mem = task.memory.toGiga() - 4
    
    """
    grep ^"${chr}	" ${probes} > ${chr}_probes_sanger.bed

    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g" \
        AnnotateIntervals \
        -R ${ref} \
        -L ${chr}_probes_sanger.bed \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${chr}_gc.tmp

    grep ^"${chr}" ${chr}_gc.tmp | awk '{ print \$1":"\$2"-"\$3"\t"\$4 }' > ${chr}_gc.txt
    """
}

process genReadCounts {
    tag { "${chr}" }
    label 'bedtools'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
    // module 'bedtools/2.30.0'
    
    input:
    path(bam_list)
    each chr

    output:
    tuple val("${chr}"), path("${chr}_canoes_reads.txt"), path("${chr}_sample_list"), emit: chr_reads_cov
    
    """
    ## CREATE A SORTED LIST OF BAM FILES FOR BEDTOOLS
    sort bam_list_unsorted.txt > bam_list_sorted.txt

    ## SUBSET THE PROBES FOR A PARTICULAR CHROMOSOME
    grep ^"${chr}	" ${probes} > ${chr}_probes_sanger.bed

    ## CREATE SAMPLE LIST FILE FOR CANOES IN R
    while read line
    do
        echo `basename \$line` | sed 's/.bam//'
    done < bam_list_sorted.txt > ${chr}_sample_list

    ## RUN BEDTOOLS FOR THE CHROMOSOME
    bedtools multicov -bams `cat bam_list_sorted.txt` -bed ${chr}_probes_sanger.bed -q 20 > ${chr}_canoes_reads.txt
    """
}

process runCANOES {
    tag { "${chr}" }
    label 'canoes'
    memory '5 GB'
    publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
    
    input:
    tuple val(chr), path(canoes_reads), path(sample_list), path(gc_content)

    output:
    tuple val(chr), path("${chr}_CNVs_pass.csv"), emit: chr_cnvs_pass
    tuple val(chr), path("${chr}_CNVs_fail.csv"), emit: chr_cnvs_fail
    tuple val(chr), path("${chr}_xcnvs.RData"), emit: chr_cnvs_rdata
    tuple val(chr), path("${chr}_CNV_genotype"), emit: chr_cnvs_geno
    tuple val(chr), path("${chr}_CNV_plots"), emit: chr_cnvs_plots
    
    """
    sed 's/chr//g; /^X/d; /^Y/d' ${canoes_reads} > ${chr}_canoes_reads_new.txt
    sed 's/chr//g; /^X/d; /^Y/d' ${gc_content} > ${chr}_gc_new.txt 
    
    run_canoes.R ${chr}_gc_new.txt ${chr}_canoes_reads_new.txt ${sample_list}
    """
}

process filterCANOESCNVs {
    tag { "ALL_CNVs" }
    label 'canoes'
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true
    
    input:
    path(cnvs)
    
    output:
    path("Sample_CNVs.csv"), emit: all_cnvs
    path("Sample_CNVs_filtered.csv"), emit: filtered_cnvs
    
    """
    grep --no-filename SAMPLE *_CNVs_pass.csv | uniq > Sample_CNVs.csv
    grep -v --no-filename SAMPLE *_CNVs_pass.csv | sort -k1,1 -k3g,3 >> Sample_CNVs.csv
    awk '{ if( (\$10>=80) && (\$10!="NA") && (\$4>=100) ) { print } }' Sample_CNVs.csv > Sample_CNVs_filtered.csv
    """
}

// process splitBAMS {
//     tag { "${sample}:${chr}" }
//     label 'samtools'
//     module 'samtools/1.17'
//     cpus 12

//     input:
//     tuple val(sample), path(bam), path(index)
//     each chr

//     output:
//     tuple val("${chr}"), path("${chr}_${bam}"), path("${chr}_${index}"), emit: chr_bams
    
//     """
//     samtools view --bam ${bam} ${chr} --output ${chr}_${bam} --threads ${task.cpus}
//     samtools index ${chr}_${bam}
//     """
// }

// process genReadCounts {
//     tag { "${chr}" }
//     label 'canoes'
//     publishDir "${outdir}/out_CANOES/${chr}", mode: 'copy', overwrite: true
//     module 'bedtools/2.30.0'
    
//     input:
//     tuple val(sample), path(bam), path(index)

//     output:
//     tuple val("${sample}"), path("${sample}.cov"), emit: samplexcanoes_reads
    
//     """
//     ls *.bam | sort > ${chr}_bam_list_sorted.txt
//     bedtools multicov -bams `cat ${chr}_bam_list_sorted.txt` -bed ${probes} -q 20 | grep ^"${chr}	" > ${chr}_canoes.reads.txt
//     while read line
//     do
//         echo `basename \$line` | sed 's/.bam//; s/${chr}_//'
//     done < ${chr}_bam_list_sorted.txt > ${chr}_sample_list
//     """
// }

// java -Xmx2000m -Djava.io.tmpdir=TEMP -jar /home/phelelani/applications/gatk-2.1-9/GenomeAnalysisTK.jar \
//     -T GCContentByInterval \
//     -L ${chr}_probes_sanger.bed \
//     -R ${ref} \
//     -o ${chr}_gc.txt

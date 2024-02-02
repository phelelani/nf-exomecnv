#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// REQUIRED FILES
ref               = file(params.ref, type: 'file')
priors            = file(params.priors, type: 'file')
indelible_conf    = file(params.indelible_conf, type: 'file')
outdir            = file(params.outdir, type: 'dir')

// 1. THE FETCH COMMAND EXTRACTS THE READS FROM THE BAM FILE, IT TAKES 2 ARGUMENTS:
process run_Fetch {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.sc_reads"), path(bam), path(bai), emit: sc_reads
    
    """
    export REF_PATH=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=/dataB/aux/38/samtools_ref_cache/%2s/%2s/%s
    indelible.py fetch --config ${indelible_conf} --i ${bam} --o ${bam}.sc_reads
    """
}

// 2. THE AGGREGATE MERGES INFORMATION ACROSS READS TOWARDS A POSITION-LEVEL VIEW OF THE DATA:
process run_Aggregate {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(sc_read), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${bam}.counts"), emit: counts
    
    """
    indelible.py aggregate --i ${sc_read} --b ${bam} --o ${bam}.counts --r ${ref} --config ${indelible_conf}
    """
}

// 3. THE SCORE COMMAND SCORES POSITIONS BASED ON THE READ INFORMATION AND SEQUENCE CONTEXT:
process run_Score {
    tag { sample }
    label 'indelible'
    
    input:
    tuple val(sample), path(count)

    output:
    tuple val(sample), path("${count}.scored"), emit: scores
    path("${count}.scored"), emit: database_in
    
    """
    indelible.py score --i ${count} --o ${count}.scored --config ${indelible_conf}
    """
}

// 4. THE DATABASE COMMAND GENERATES THE ALLELE FREQUENCY AND BREAKPOINT DATABASE REQUIRED FOR THE NEXT STEP – ANNOTATE
process run_Database {
    tag { "Indel_DB" }
    label 'indelible'
    cpus 6
    publishDir "${outdir}/out_INDELIBLE/database", mode: 'copy', overwrite: true

    input:
    path(score)

    output:
    path("InDelible_db.tsv"), emit: indel_database
    
    """
    ls *.scored > scores.txt
    indelible.py database --f scores.txt --o InDelible_db.tsv --r ${ref} --priors ${priors} --config ${indelible_conf} --tb ${task.cpus}
    """
}

// 5. THE ANNOTATE COMMAND ENRICHES THE RESULT WITH GENE/EXON ANNOTATIONS AND MERGES THE DATABASE RESULTS WITH THE POSITION FILE:
process run_Annotate {
    tag { "${score}" }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/annotations", mode: 'copy', overwrite: true

    input:
    path(database)
    tuple val(sample), path(score)

    output:
    tuple val(sample), path("${score}.annotated"), emit: annotated
    
    """
    indelible.py annotate --i ${score} --o ${score}.annotated --d ${database} --config ${indelible_conf}
    """
}

// 6. ONE CAN THEN LOOK FOR DE NOVO MUTATION EVENTS USING THE DENOVO COMMAND:
// TRIO
process run_DenovoTrio {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_trio", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo

    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// MOM
process run_DenovoMom {
    tag { sample }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_mom", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(mom_bam), path(mom_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_mom

    """
    indelible.py denovo --c ${annotation} --m ${mom_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

// DAD
process run_DenovoDad {
    label 'indelible'
    tag { sample }

    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_dad", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_bam), path(child_bai), path(dad_bam), path(dad_bai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_dad

    """
    indelible.py denovo --c ${annotation} --p ${dad_bam} --o ${annotation}.denovo.tsv --config ${indelible_conf}
    """    
}

process filterINDELIBLE {
    tag { 'filter_cnvs' }
    label 'indelible'
    publishDir "${outdir}/out_INDELIBLE/filtered", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(annotated)
    
    output:
    path("${sample}.annotated.filtered.tsv"), emit: filtered_cnvs
    
    """
    awk '{ if ((\$39 < 2) && (\$40 < 2)) { print } }' ${annotated} > ${sample}.annotated.filtered.tsv
    """
}

#!/usr/bin/env nextflow

// REQUIRED FILES
ref               = file(params.ref, type: 'file')
priors            = file(params.priors, type: 'file')
indelible_config  = file(params.config, type: 'file')
outdir            = file(params.outdir, type: 'dir')

// 1. THE FETCH COMMAND EXTRACTS THE READS FROM THE BAM FILE, IT TAKES 2 ARGUMENTS:
process run_Fetch {
    tag { sample }
    errorStrategy 'ignore'
    
    input:
    tuple val(sample), path(cram), path(crai)
    
    output:
    tuple val(sample), path("${cram}.sc_reads"), path(cram), path(crai), emit: sc_reads
    
    """
    /bin/hostname
    indelible.py fetch --config ${config} --i ${cram} --o ${cram}.sc_reads
    """
}

// 2. THE AGGREGATE MERGES INFORMATION ACROSS READS TOWARDS A POSITION-LEVEL VIEW OF THE DATA:
process run_Aggregate {
    tag { sample }
    
    input:
    tuple val(sample), path(sc_read), path(cram), path(crai)
    
    output:
    tuple val(sample), path("${cram}.counts"), emit: counts
    
    """
    indelible.py aggregate --i ${sc_read} --b ${cram} --o ${cram}.counts --r ${ref} --config ${config}
    """
}

// 3. THE SCORE COMMAND SCORES POSITIONS BASED ON THE READ INFORMATION AND SEQUENCE CONTEXT:
process run_Score {
    tag { sample }

    input:
    tuple val(sample), path(count)

    output:
    tuple val(sample), path("${count}.scored"), emit: scores
    path("${count}.scored"), emit: database_in
    
    """
    indelible.py score --i ${count} --o ${count}.scored --config ${config}
    """
}

// 4. THE DATABASE COMMAND GENERATES THE ALLELE FREQUENCY AND BREAKPOINT DATABASE REQUIRED FOR THE NEXT STEP – ANNOTATE
process run_Database {
    tag { "Indel_DB" }
    cpus 6
    publishDir "${outdir}/out_INDELIBLE/database", mode: 'copy', overwrite: true

    input:
    path(score)

    output:
    path("InDelible_db.tsv"), emit: indel_database
    
    """
    ls *.scored > scores.txt
    indelible.py database --f scores.txt --o InDelible_db.tsv --r ${ref} --priors ${priors} --config ${config} --tb 6
    """
}

// 5. THE ANNOTATE COMMAND ENRICHES THE RESULT WITH GENE/EXON ANNOTATIONS AND MERGES THE DATABASE RESULTS WITH THE POSITION FILE:
process run_Annotate {
    tag { "${score}" }
    publishDir "${outdir}/out_INDELIBLE/annotations", mode: 'copy', overwrite: true

    input:
    path(database)
    tuple val(sample), path(score)

    output:
    tuple val(sample), path("${score}.annotated"), emit: annotated
    
    """
    indelible.py annotate --i ${score} --o ${score}.annotated --d ${database} --config ${config}
    """
}

// 6. ONE CAN THEN LOOK FOR DE NOVO MUTATION EVENTS USING THE DENOVO COMMAND:
// TRIO
process run_DenovoTrio {
    tag { sample }
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_trio", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_cram), path(child_crai), path(mom_cram), path(mom_crai), path(dad_cram), path(dad_crai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo

    """
    indelible.py denovo --c ${annotation} --m ${mom_cram} --p ${dad_cram} --o ${annotation}.denovo.tsv --config ${config}
    """    
}

// MOM
process run_DenovoMom {
    tag { sample }
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_mom", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_cram), path(child_crai), path(mom_cram), path(mom_crai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_mom

    """
    indelible.py denovo --c ${annotation} --m ${mom_cram} --o ${annotation}.denovo.tsv --config ${config}
    """    
}

// DAD
process run_DenovoDad {
    tag { sample }
    publishDir "${outdir}/out_INDELIBLE/denovo_annotation_dad", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(child_cram), path(child_crai), path(dad_cram), path(dad_crai), path(annotation)
    
    output:
    tuple val(sample), path("${annotation}.denovo.tsv"), emit: indelible_denovo_dad

    """
    indelible.py denovo --c ${annotation} --p ${dad_cram} --o ${annotation}.denovo.tsv --config ${config}
    """    
}

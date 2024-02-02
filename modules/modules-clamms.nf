#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref         = file(params.ref, type: 'file')
probes      = file(params.probes, type: 'file')
mappability = file(params.mappability, type: 'file')
special_reg = file(params.special_reg, type: 'file')
outdir      = file(params.outdir, type: 'dir')
sexinfo     = file(params.sexinfo, type: 'file')

process generateWindows {
    tag { 'generate_windows' }
    memory '11 GB'
    module 'bedtools/2.30.0'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true
    
    output:
    path("windows.bed"), emit: windows
    
    """
    export INSERT_SIZE=200
    export CLAMMS_DIR=/home/phelelani/applications/clamms

    sort -k1,1 -k2,2n ${probes} > targets_sorted.bed
    
    \$CLAMMS_DIR/annotate_windows.sh targets_sorted.bed ${ref} ${mappability} \
        \$INSERT_SIZE ${special_reg} > windows.bed
    """
}


process samtoolsDOC {
    tag { sample }
    maxForks 20
    memory '11 GB'
    module 'samtools/1.15'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam), path(bai)
    path(windows)
    
    output:
    tuple val(sample), path("${sample}.coverage.bed"), emit: coverage
    
    """
    samtools bedcov -Q 30 ${windows} ${bam} | awk '{ printf \"%s\\t%d\\t%d\\t%.6g\\n\", \$1, \$2, \$3, \$NF/(\$3-\$2); }' > ${sample}.coverage.bed
    """
}

process normalizeDOC {
    tag { sample }
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(coverage)
    path(windows)
    
    output:
    path("${sample}.norm.cov.bed"), emit: coverage_norm
    tuple val(sample), path("${sample}.norm.cov.bed"), emit: coverage_norm_set

    """
    export CLAMMS_DIR=/home/phelelani/applications/clamms

    $CLAMMS_DIR/normalize_coverage ${sample}.coverage.bed ${windows} | sed 's/^chr//g' > ${sample}.norm.cov.bed
    """
}

process trainModels {
    tag { 'train_models' }
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    path(coverage_norm)
    path(windows)
    
    output:
    path("models.bed"), emit: models
    
    """
    export CLAMMS_DIR=/home/phelelani/applications/clamms
    sed 's/^chr//g' ${windows} > windows_new
    $CLAMMS_DIR/fit_models ${sexinfo} windows_new | sed 's/^chr//g' > models.bed
    """
}

// grep "^Y" \$FILE | awk '{ x += \$4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'

process callCNVs {
    tag { sample }
    memory '11 GB'
    errorStrategy 'ignore'
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(norm_cov)
    path(models)
    
    output:
    tuple val(sample), path("${sample}.cnv.bed"), emit: cnvs

    """
    export CLAMMS_DIR=/home/phelelani/applications/clamms
    sex=`grep ${sample} ${sexinfo} | cut -f 2`
    $CLAMMS_DIR/call_cnv ${norm_cov} ${models} --sex \$sex > ${sample}.cnv.bed
    """
}

process filterCLAMMSCNVs {
    tag { 'filter_cnvs' }
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true
    
    input:
    path(cnvs)
    
    output:
    tuple path("samples.cnv.bed"), path("samples.cnv.filtered.bed"), emit: filtered_cnvs
    
    """
    cat *.cnv.bed > samples.cnv.bed
    awk '{ if( (\$9>=500) && (\$10>0) ) { print } }' samples.cnv.bed > samples.cnv.filtered.bed
    """
}

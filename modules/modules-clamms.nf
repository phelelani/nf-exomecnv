#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref           = file(params.ref, type: 'file')
probes        = file(params.probes, type: 'file')
interval_list = file(params.interval_list, type: 'file')
mappability   = file(params.mappability, type: 'file')
special_reg   = file(params.special_reg, type: 'file')
outdir        = file(params.outdir, type: 'dir')
sexinfo       = file(params.sexinfo, type: 'file')

process generateWindows {
    tag { 'generate_windows' }
    label 'clamms'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/sample_coverage", mode: 'copy', overwrite: false
    
    output:
    path("windows.bed"), emit: windows
    
    """
    export INSERT_SIZE=200

    sort -k1,1 -k2,2n ${probes} > targets_sorted.bed
    
    \$CLAMMS_DIR/annotate_windows.sh targets_sorted.bed ${ref} ${mappability} \$INSERT_SIZE ${special_reg} > windows.bed
    """
}


process samtoolsDOC {
    tag { sample }
    label 'clamms'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/sample_coverage", mode: 'copy', overwrite: false

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
    label 'clamms'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/sample_coverage", mode: 'copy', overwrite: false

    input:
    tuple val(sample), path(coverage)
    path(windows)
    
    output:
    tuple val(sample), path("${sample}.norm.cov.bed"), emit: norm_coverage

    """
    \$CLAMMS_DIR/normalize_coverage ${sample}.coverage.bed ${windows} | sed 's/^chr//g' > ${sample}.norm.cov.bed
    wc -l ${sample}.norm.cov.bed 
    """
}

process getPicardQCMetrics {
    tag { sample }
    label 'picard'
    publishDir "${outdir}/out_CLAMMS/qc_metrics", mode: 'copy', overwrite: false

    input:
    tuple val(sample), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${sample}.hs_metrics.txt"), emit: qc_metrics

    script:
    mem = task.memory.toGiga() - 4
    
    """
    java -XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g -jar /usr/picard/picard.jar \
        CollectHsMetrics \
        --INPUT ${bam} \
        --OUTPUT ${sample}.hs_metrics.txt \
        --TARGET_INTERVALS	${interval_list} \
        --BAIT_INTERVALS ${interval_list} \
        --REFERENCE_SEQUENCE ${ref}
    """
}

process getPicardMeanInsertSize {
    tag { sample }
    label 'picard'
    publishDir "${outdir}/out_CLAMMS/qc_metrics", mode: 'copy', overwrite: false

    input:
    tuple val(sample), path(bam), path(bai)
    
    output:
    tuple val(sample), path("${sample}.insert_size_metrics.txt"), emit: ins_size_metrics
    tuple val(sample), path("${sample}.insert_size_histogram.pdf"), emit: ins_size_hist

    script:
    mem = task.memory.toGiga() - 4
    
    """
    java -XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g -jar /usr/picard/picard.jar \
        CollectInsertSizeMetrics \
        --INPUT ${bam} \
        --OUTPUT ${sample}.insert_size_metrics.txt \
        --Histogram_FILE ${sample}.insert_size_histogram.pdf
    """
}

process combinePicardQCMetrics {
    tag { 'combine_qc_metrics' }
    label 'picard'
    publishDir "${outdir}/out_CLAMMS/qc_metrics", mode: 'copy', overwrite: false

    input:
    path(qc_metrics)
    
    output:
    path("qcs_metrics"), emit: qcs_metrics
    
    """
    combine_picard_qc_metrics.sh
    """
}

process createPCAData {
    tag { 'create_pca' }
    label 'clamms'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/custom_reference_panels", mode: 'copy', overwrite: false

    input:
    path(coverage_norm)
    
    output:
    path("pca.coordinates.txt"), emit: pca_data
    
    """
    custom_ref_panel.sh
    """
}

process createCustomRefPanel {
    tag { 'train_models' }
    label 'R'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/custom_reference_panels", mode: 'copy', overwrite: false

    input:
    path(coverage_norm)
    path(pca_data)
    path(sample_metrics)
    
    output:
    path("*.pdf"), emit: plots
    path("*.ref.panel.files.txt"), emit: ref_panel
    
    """
    custom_ref_panel.R ${pca_data} ${sample_metrics} ${sexinfo}
    """
}

process trainModels {
    tag { sample }
    label 'clamms'
    memory '11 GB'
    publishDir "${outdir}/out_CLAMMS/models", mode: 'copy', overwrite: false

    input:
    tuple val(sample), path(ref_panel)
    path(windows)
    path(coverage_norm)
    
    output:
    tuple val(sample), path("${sample}.models.bed"), emit: sample_models
    
    """
    sed 's/chr//g' ${windows} > windows.new.bed

   \$CLAMMS_DIR/fit_models ${ref_panel} windows.new.bed > ${sample}.models.bed
    """
}

//find . -type f -iname "*coverage.bed" -exec sh -c 'size=`cut -d'\t' -f 1 "$1" | sort | uniq | wc -l`; if [[ $size == 196768 ]]; then : #echo -e "$1\tPASS: $size"; else echo -e "$1\tFAIL: $size"; fi' - {} \; 
// CHECK WHAT IS GOING ON HERE!
process callCNVs {
    tag { sample }
    label 'clamms'
    memory '11 GB'
    errorStrategy 'ignore'
    publishDir "${outdir}/out_CLAMMS/cnv_calls", mode: 'copy', overwrite: false

    input:
    tuple val(sample), path(models), path(norm_cov)
    
    output:
    tuple val(sample), path("${sample}.cnv.bed"), emit: cnvs

    """
    sex=`grep ${sample} ${sexinfo} | cut -f 2`
    \$CLAMMS_DIR/call_cnv ${norm_cov} ${models} --sex \$sex > ${sample}.cnv.bed
    """
}

process filterCLAMMSCNVs {
    tag { 'filter_cnvs' }
    publishDir "${outdir}/out_CLAMMS/cnv_calls_filtered", mode: 'copy', overwrite: false
    
    input:
    path(cnvs)
    
    output:
    tuple path("samples.cnv.bed"), path("samples.cnv.filtered.bed"), emit: filtered_cnvs
    
    """
    cat *.cnv.bed > samples.cnv.bed
    awk '{ if( (\$9>=500) && (\$10>0) ) { print } }' samples.cnv.bed > samples.cnv.filtered.bed
    """
}

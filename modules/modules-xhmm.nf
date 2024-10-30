#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref       = file(params.ref, type: 'file')
probes    = file(params.probes, type: 'file')
xhmm_conf = file(params.xhmm_conf, type: 'file')
outdir    = file(params.outdir, type: 'dir')

// SPLIT BAM FILES TO GROUPS
process groupBAMs {
    tag { 'group_bams' }

    input:
    path(bams)

    output:
    path("bam_group_*"), emit: bam_groups
    
    """
    sort bam_list_unsorted.txt > bam_list_sorted.txt
    split -l 5 bam_list_sorted.txt --numeric-suffixes=000 --suffix-length=3 --additional-suffix=.list bam_group_
    """
}

// RUN GATK FOR DEPTH OF COVERAGE (FOR SAMPLES IN EACH GROUP):
process gatkDOC {
    tag { "${group}" }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    tuple val(group), path(list)
    
    output:
    tuple val(group), path("${group}.DATA.sample_interval_summary"), emit: bam_group_doc

    script:
    mem = task.memory.toGiga() - 4
    
    """
    gatk --java-options "-XX:+UseSerialGC -Xss456k -Xms2g -Xmx${mem}g -Djava.io.tmpdir=\$PWD" \
        DepthOfCoverage \
        -I ${list} \
        -L ${probes} \
        -R ${ref} \
        --max-depth-per-sample 5000 \
        --verbosity INFO \
        --omit-depth-output-at-each-base true \
        --omit-locus-table true \
        --min-base-quality 0 \
        --read-filter MappingQualityReadFilter \
        --minimum-mapping-quality 20 \
        --start 1 --stop 5000 --nBins 200 \
        --include-ref-n-sites true \
        --count-type COUNT_READS \
        --output ${group}.DATA
    """
}

// COMBINES GATK DEPTH-OF-COVERAGE OUTPUTS FOR MULTIPLE SAMPLES (AT SAME LOCI):
process combineDOC {
    tag { 'combine_doc' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(list)
    
    output:
    path("DATA.RD.txt"), emit: combined_doc
    
    """
    for i in *.sample_interval_summary; do sed 's/,/	/g' \$i > \${i%.sample_interval_summary}.fixed_sample_interval_summary ; done
    ls *.fixed_sample_interval_summary > the_list
    xhmm --mergeGATKdepths -o DATA.RD.txt --GATKdepthsList the_list
    """
}

// OPTIONALLY, RUN GATK TO CALCULATE THE PER-TARGET GC CONTENT AND CREATE A LIST OF THE TARGETS WITH EXTREME GC CONTENT:
process calcGC_XHMM {
    tag { 'calc_dc' }
    label 'gatk'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    output:
    path("DATA.locus_GC.txt"), emit: gc_targets
    path("extreme_gc_targets.txt"), emit: extreme_gc_targets
    
    """
    gatk AnnotateIntervals \
        -R ${ref} \
        -L ${probes} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O gc.tmp

    grep "^chr" gc.tmp | awk '{ print \$1":"\$2"-"\$3"\t"\$4 }' > DATA.locus_GC.txt
    cat DATA.locus_GC.txt | awk '{ if (\$2 < 0.1 || \$2 > 0.9) print \$1 }' > extreme_gc_targets.txt
    """
}

// java -Xmx3072m -Djava.io.tmpdir=TEMP -jar /home/phelelani/applications/gatk-2.1-9/GenomeAnalysisTK.jar \
//     -T GCContentByInterval \
//     -L ${probes} \
//     -R ${ref} \
//     -o DATA.locus_GC.txt

// process calc_PLINKSeq {
//     tag {}
//     publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

//     output:

//     input:

//     """
//     awk 'BEGIN{OFS="\t"; print "#CHR\tBP1\tBP2\tID"} {chr=$1; if (match(chr,"chr")==0) {chr="chr"chr} bp1=$2; bp2=bp1; if (length($2-$3) > 1) {bp2=$3}print chr,bp1,bp2,NR}' ${probes} > EXOME.targets.reg



//     ./sources/scripts/interval_list_to_pseq_reg EXOME.interval_list > ./EXOME.targets.reg
//     pseq . loc-load --locdb ./EXOME.targets.LOCDB --file ./EXOME.targets.reg --group targets --out ./EXOME.targets.LOCDB.loc-load
    
//     pseq . loc-stats --locdb ./EXOME.targets.LOCDB --group targets --seqdb ./seqdb | \
//     awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste ./EXOME.interval_list - | \
//     awk '{print $1"\t"$2}' > ./DATA.locus_complexity.txt

//     cat ./DATA.locus_complexity.txt | awk '{if ($2 > 0.25) print $1}' > ./low_complexity_targets.txt
//     """
// }

// FILTERS SAMPLES AND TARGETS AND THEN MEAN-CENTERS THE TARGETS:
process filterSamples {
    tag { 'filter_samples' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true    

    input:
    path(combined_doc)
    path(extreme_gc_targets)
    
    output:
    path("DATA.filtered_centered.RD.txt"), emit: filtered_centered
    path("DATA.filtered_centered.RD.txt.filtered_targets.txt"), emit: excluded_filtered_targets
    path("DATA.filtered_centered.RD.txt.filtered_samples.txt"), emit: excluded_filtered_samples
    
    """
    xhmm --matrix -r DATA.RD.txt --centerData --centerType target \
        -o DATA.filtered_centered.RD.txt \
        --outputExcludedTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeTargets extreme_gc_targets.txt \
        --minTargetSize 10 --maxTargetSize 10000 \
        --minMeanTargetRD 10 --maxMeanTargetRD 500 \
        --minMeanSampleRD 25 --maxMeanSampleRD 200 \
        --maxSdSampleRD 150
   """
}

// RUNS PCA ON MEAN-CENTERED DATA:
process runPCA {
    tag { 'run_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    
    output:
    tuple val("pca_data"), path("DATA.RD_PCA*"), emit: pca_data

    """
    xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA
    """
}

// NORMALIZES MEAN-CENTERED DATA USING PCA INFORMATION:
process normalisePCA {
    tag { 'norm_pca' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(filtered_centered)
    tuple val(pca), path(pca_data)
    
    output:
    path("DATA.PCA_normalized.txt"), emit: data_pca_norm

    """
    xhmm --normalize -r DATA.filtered_centered.RD.txt \
        --PCAfiles DATA.RD_PCA \
        --normalizeOutput DATA.PCA_normalized.txt \
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
    """
}

// FILTERS AND Z-SCORE CENTERS (BY SAMPLE) THE PCA-NORMALIZED DATA:
process filterZScore {
    tag { 'filter_zscore' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(data_pca_norm)
    
    output:
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt"), emit: pca_norm_zscore
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"), emit: excluded_zscore_targets
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"), emit: excluded_zscore_samples
    
    """
    xhmm --matrix -r DATA.PCA_normalized.txt \
        --centerData --centerType sample --zScoreData \
        -o DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        --outputExcludedTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --outputExcludedSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        --maxSdTargetRD 30
    """
}

// FILTERS ORIGINAL READ-DEPTH DATA TO BE THE SAME AS FILTERED, NORMALIZED DATA:
process filterRD {
    tag { 'filter_rd' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(combined_doc)
    path(excluded_filtered_targets)
    path(excluded_zscore_targets)
    path(excluded_filtered_samples)
    path(excluded_zscore_samples)
    
    output:
    path("DATA.same_filtered.RD.txt"), emit: orig_filtered

    """
    xhmm --matrix -r DATA.RD.txt \
        --excludeTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --excludeTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
        --excludeSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
        -o DATA.same_filtered.RD.txt
    """
}

// DISCOVERS CNVS IN NORMALIZED DATA:
process discoverCNVs {
    tag { 'discover_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    
    output:
    path("*"), emit: cnvs
    
    """
    xhmm --discover -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -c DATA.xcnv -a DATA.aux_xcnv -s DATA
    """
}

// GENOTYPES DISCOVERED CNVS IN ALL SAMPLES:
process genotypeCNVs {
    tag { 'genotype_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    path(cnvs)
    
    output:
    path("*"), emit: genotypes
    
    """
    xhmm --genotype -p ${xhmm_conf} \
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        -R DATA.same_filtered.RD.txt \
        -g DATA.xcnv -F ${ref} \
        -v DATA.vcf
    """
}

process filterXHMMCNVs {
    tag { 'filter_cnvs' }
    label 'xhmm'
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true
    
    input:
    path(cnvs)
    
    output:
    path("DATA_filtered.xcnv"), emit: filtered_cnvs
    
    """
    awk '{ if( (\$9>=60) && (\$10>=60) && (\$11>=60) ) { print } }' DATA.xcnv > DATA_filtered.xcnv
    """
}

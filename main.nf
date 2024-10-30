#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// INPUT FILES FOR INDELIBLE
Channel.fromFilePairs([params.crams + '/*{.cram,.cram.crai}'])
    .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }
    .filter { it -> it[1] =~ '_01_1' }
    .set { crams }

Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: 6)
    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5] ] }
    .set { cram_trios }

Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
    .filter { it -> it[1].size() == 4 }
    .filter { it -> it[1][3] =~ '02_2.cram' }
    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
    .set { cram_mom }

Channel.fromFilePairs([params.crams + '/*{_01_1,_02_2,_03_3}*'], size: -1)
    .filter { it -> it[1].size() == 4 }
    .filter { it -> it[1][3] =~ '03_3.cram' }
    .map { it -> [ it[0], it[1][0], it[1][1], it[1][2], it[1][3] ] }
    .set { cram_dad }

// OUTPUT DIR
outdir    = file(params.outdir, type: 'dir')
workflow  = params.workflow
outdir.mkdir()

include { run_Fetch; run_Aggregate; run_Score; run_Database;
         run_Annotate; run_DenovoTrio; run_DenovoMom; run_DenovoDad; filterINDELIBLE } from './modules/modules-indelible.nf'
include { genReadCounts; calcGC_CANOES; runCANOES; filterCANOESCNVs } from './modules/modules-canoes.nf'
include { groupBAMs; gatkDOC; combineDOC; calcGC_XHMM; filterSamples; runPCA;
         normalisePCA; filterZScore; filterRD; discoverCNVs; genotypeCNVs; filterXHMMCNVs } from './modules/modules-xhmm.nf'
include { generateWindows; samtoolsDOC; normalizeDOC; createPCAData; getPicardQCMetrics; getPicardMeanInsertSize; combinePicardQCMetrics;
         createCustomRefPanel; trainModels; callCNVs; filterCLAMMSCNVs } from './modules/modules-clamms.nf'

// INDELIBLE WORKFLOW
workflow RUN_INDELIBLE {
    take:
    crams
    
    main:
    run_Fetch(crams)
    run_Aggregate(run_Fetch.out.sc_reads)
    run_Score(run_Aggregate.out.counts)
    run_Database(run_Score.out.database_in.collect())
    run_Annotate(run_Database.out.indel_database, run_Score.out.scores)
    run_DenovoTrio(cram_trios.join(run_Annotate.out.annotated))
    run_DenovoMom(cram_mom.join(run_Annotate.out.annotated))
    run_DenovoDad(cram_dad.join(run_Annotate.out.annotated))
    filterINDELIBLE(run_Annotate.out.annotated)
}

// CANOES WORKFLOW
workflow RUN_CANOES {
    take:
    bams
    chroms

    main:
    calcGC_CANOES(chroms)
    bams
        .collectFile() { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] }
        .set { bam_list }
    genReadCounts(bam_list,chroms)
    genReadCounts.out.chr_reads_cov
        .join(calcGC_CANOES.out.chr_gc_content)
        .set { chr_canoes_input }
    runCANOES(chr_canoes_input)
    runCANOES.out.chr_cnvs_pass
        .map { it -> it[1] }
        .collect()
        .set { all_cnvs_pass }
    filterCANOESCNVs(all_cnvs_pass)
    filterCANOESCNVs.out.filtered_cnvs.view()
}

// XHMM WORKFLOW
workflow RUN_XHMM {
    take:
    bams

    main:
    groupBAMs(bams.collectFile () { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] })
    gatkDOC(groupBAMs.out.bam_groups.flatMap().map { it -> [it.name[0..-6], it] })
    combineDOC(gatkDOC.out.bam_group_doc.collect { it -> it[1] })
    calcGC_XHMM()
    filterSamples(combineDOC.out.combined_doc, calcGC_XHMM.out.extreme_gc_targets)
    runPCA(filterSamples.out.filtered_centered)
    normalisePCA(filterSamples.out.filtered_centered, runPCA.out.pca_data)
    filterZScore(normalisePCA.out.data_pca_norm)
    filterRD(combineDOC.out.combined_doc,
             filterSamples.out.excluded_filtered_targets,filterSamples.out.excluded_filtered_samples,
             filterZScore.out.excluded_zscore_targets, filterZScore.out.excluded_zscore_samples )
    discoverCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore)
    genotypeCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore, discoverCNVs.out.cnvs)
    filterXHMMCNVs(discoverCNVs.out.cnvs.collect())
}

// CLAMMS WORKFLOW
workflow RUN_CLAMMS {
    take:
    bams

    main:
    generateWindows()
    samtoolsDOC(bams, generateWindows.out.windows)
    normalizeDOC(samtoolsDOC.out.coverage, generateWindows.out.windows)
    normalizeDOC.out.norm_coverage
        .map { it -> it[1] }
        .collect()
        .set { norm_coverage_files }
    createPCAData(norm_coverage_files)
    getPicardQCMetrics(bams)
    getPicardMeanInsertSize(bams)
    getPicardQCMetrics.out.qc_metrics
        .join(getPicardMeanInsertSize.out.ins_size_metrics, by:0)
        .map { it -> [ it[1], it[2] ] }
        .flatten()
        .collect()
        .set { picard_metrics }
    combinePicardQCMetrics(picard_metrics)
    createCustomRefPanel(norm_coverage_files,createPCAData.out.pca_data,combinePicardQCMetrics.out.qcs_metrics)
    createCustomRefPanel.out.ref_panel
        .flatten()
        .map { it -> [ "${it.baseName.replaceAll('.ref.panel.files','')}", it ] }
        .set { for_training }
    trainModels(for_training, generateWindows.out.windows, norm_coverage_files)
    trainModels.out.sample_models
        .join(normalizeDOC.out.norm_coverage)
        .set { cllin }
    callCNVs(cllin)
    filterCLAMMSCNVs(callCNVs.out.cnvs.map { it -> it[1] }.collect())
}

// PICK AND CHOOSE
workflow {
    switch (workflow) {
        case['indelible']:
            RUN_INDELIBLE(crams)
            break
            // =====
        case['canoes']:
            chroms = (1..22).toList().collect { 'chr' + "${it}" } //+ ["chrX", "chrY","chrM"]
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { bams }
            RUN_CANOES(bams,chroms)
            break
            // =====
        case['xhmm']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { bams }
            RUN_XHMM(bams)
            break
            // =====
        case['clamms']:
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { bams }
            RUN_CLAMMS(bams)
            break
            // =====
        case['all']:
            chroms = (1..22).toList().collect { 'chr' + "${it}" } //+ ["chrX", "chrY","chrM"]
            Channel.fromPath(params.samplesheet_bams)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.BAM}", "${row.BAM}".replaceAll("\\b.bam\\b",".bam.bai") ] }
                .set { bams }
            //RUN_INDELIBLE(crams)
            RUN_CANOES(bams,chroms)
            RUN_XHMM(bams)
            RUN_CLAMMS(bams)
            break
            // =====}
        default:
            exit 1, """
OOOPS!! SEEMS LIE WE HAVE A WORFLOW ERROR!

No workflow \'mode\' give! Please use one of the following options for workflows:
    --mode indelible    // To run the INDELIBLE WORKFLOW
    --mode canoes       // To run the CANOES
    --mode canoes       // To run the CANOES workflow
    --mode all          // To run all workflows
"""
            break
    }
}

#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// INPUT FILES FOR CANOES & XHMM
Channel.fromFilePairs([params.bams + '/*{.bam,.bam.bai}'])
    .map { it -> [ it[0], it[1].find { it =~ '.bam$' }, it[1].find { it =~ '.bai$' }] }
    .set { bams }

// bams.collectFile () { item -> [ 'bam_list_unsorted.txt', "${item.get(1)[0]}" + '\n' ] }
//     .set { bam_list }

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
outdir = file(params.outdir, type: 'dir')
outdir.mkdir()
workflow                  = params.workflow

include { run_Fetch; run_Aggregate; run_Score; run_Database;
         run_Annotate; run_DenovoTrio; run_DenovoMom; run_DenovoDad; filterINDELIBLE } from './modules/modules-indelible.nf'
include { genReadCounts; calcGC_CANOES; runCANOES; filterCANOESCNVs } from './modules/modules-canoes.nf'
include { groupBAMs; gatkDOC; combineDOC; calcGC_XHMM; filterSamples; runPCA;
         normalisePCA; filterZScore; filterRD; discoverCNVs; genotypeCNVs; filterXHMMCNVs } from './modules/modules-xhmm.nf'
include { generateWindows; samtoolsDOC; normalizeDOC; trainModels; callCNVs; filterCLAMMSCNVs } from './modules/modules-clamms.nf'

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
    
    main:
    genReadCounts(bams.collectFile () { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] })
    calcGC_CANOES()
    runCANOES(genReadCounts.out.canoes_reads, calcGC_CANOES.out.gc_content)
    filterCANOESCNVs(runCANOES.out.cnvs)
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
    trainModels(normalizeDOC.out.coverage_norm.collect(), generateWindows.out.windows)
    callCNVs(normalizeDOC.out.coverage_norm_set, trainModels.out.models)
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
            RUN_CANOES(bams)
            break
            // =====
        case['xhmm']:
            RUN_XHMM(bams)
            break
            // =====
        case['clamms']:
            RUN_CLAMMS(bams)
            break
            // =====
        case['all']:
            // RUN_INDELIBLE(crams)
            RUN_CANOES(bams)
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

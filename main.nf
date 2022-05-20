#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.crams          = '/dataG/ddd/data/batch_1/cram_ddd_id'
params.bams           = '/external/diskC/ddd/data/batch_1'
params.outdir         = '/home/phelelani/projects/nadja/test_exomecall'
params.ref            = '/home/phelelani/projects/nadja/data/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.priors         = '/home/phelelani/projects/nadja/data/indelible/Indelible_db_10k.hg38.bed'
params.indelible_conf = '/home/phelelani/projects/nadja/data/indelible/config.yml'
params.xhmm_conf      = '/home/phelelani/projects/nadja/data/xhmm/params.txt'
params.probes         = '/home/phelelani/projects/nadja/data/canoes/probes_sanger.bed'

// INPUT FILES FOR CANOES & XHMM
Channel.fromFilePairs([params.bams + '/*{.bam,.bam.bai}'])
    .set { bams }

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

include { run_Fetch; run_Aggregate; run_Score; run_Database;
         run_Annotate; run_DenovoTrio; run_DenovoMom; run_DenovoDad } from './modules/modules-indelible.nf'

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
    }


workflow RUN_CANOES {
    take:
        bams
    main:
        genReadCounts(bams)
        calcGC()
        runCANOES(genReadCounts.out.canoes_reads, calcGC.out.)
}

workflow RUN_XHMM {
    take:
        bams
    main:
        groupBAMs(bams.collectFile () { item -> [ 'bamlist_unsorted.txt', "${item.get(1)[0]}" + '\n' ] })
        calcDOC(groupBAMs.out.bam_groups.flatMap().map { it -> [it.name[0..-6], it] })
        combineDOC(calcDOC.out.bam_group_doc.collect { it -> it[1] })
        calcGC()
        filterSamples(combineDOC.out.combined_doc, calcGC.out.extreme_gc_targets)
        runPCA(filterSamples.out.filtered_centered)
        normalisePCA(filterSamples.out.filtered_centered, runPCA.out.pca_data)
        filterZScore(normalisePCA.out.data_pca_norm)
        filterRD(combineDOC.out.combined_doc,
                 filterSamples.out.excluded_filtered_targets,filterSamples.out.excluded_filtered_samples,
                 filterZScore.out.excluded_zscore_targets, filterZScore.out.excluded_zscore_samples )
        discoverCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore)
        genotypeCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore, discoverCNVs.out.cnvs)

}

workflow {
    switch (mode) {
        case['indelible']:
            RUN_INDELIBLE(crams)
            break
            // =====
        case['canoes']:
            RUN_CANOES(bams)
            break
            // =====
        case['xhmm']:
            RUN_XHMM
            break
            // =====
        default:
            
            break
            // =====}
    }
}

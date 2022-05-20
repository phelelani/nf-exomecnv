#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.crams   = '/dataG/ddd/data/batch_1/cram_ddd_id'
params.bams    = ''
params.outdir  = '/home/phelelani/projects/nadja/test_exomecall'
params.ref     = '/home/phelelani/projects/nadja/data/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.priors  = '/home/phelelani/projects/nadja/data/hg38/Indelible_db_10k.hg38.bed'
params.config  = '/home/phelelani/projects/nadja/config.yml'

// INPUT FILES FOR CANOES
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

workflow {
    RUN_INDELIBLE(crams)
}

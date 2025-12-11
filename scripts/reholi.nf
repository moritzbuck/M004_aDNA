#!/usr/bin/env nextflow

include { unpigz } from "./processes.nf"
include { merge_fastq } from "./processes.nf"
include { remove_low_complex } from "./processes.nf"
include { bowtie2 } from "./processes.nf"
include { sam2bam } from "./processes.nf"
include { bam2lca } from "./processes.nf"
include { make_tables } from "./processes.nf"
include { merge_bams } from "./processes.nf"
include { sort_bam } from "./processes.nf"
include { index_bam } from "./processes.nf"
include { check_databases } from "./processes.nf"
include { bt2_index } from "./processes.nf"
include { validate } from "./processes.nf"



workflow make_dbs{

    main :
    db_pipe = Channel.from(file(params.fasta_db_dir))
    db_pipe
        | check_databases
        | splitCsv
        | map { x -> x[0]}
        | bt2_index
        | flatten
        | map { x -> params.db_dir + "/" + x.getFileName().toString() }
        | set{ indexes }

    // indexes | map( x -> "to_join : " + x )
    //         | view

    emit :
    indexes
}

workflow mapping{

    take :
        extra_dbs
    main : 

    Channel.fromFilePairs(file(params.fastq_dir).resolve('*_{R1,R2}*.fastq.gz')).set { ch_fastq_gzs }
    Channel.fromPath(file(params.db_dir).resolve("*.rev.1.bt2l")).concat(extra_dbs).map{ x -> x.toString().split(".rev")[0] }.unique().set { made_ch_blast_dbs }
    
//    made_ch_blast_dbs.view()

    ch_fastq_gzs 
        | unpigz 
        | merge_fastq 
        | remove_low_complex
        | combine( made_ch_blast_dbs ) 
        | bowtie2 
        | sam2bam
        | groupTuple()
        | merge_bams
        | sort_bam
        | index_bam
        | bam2lca
        | map { v -> tuple(v[1], v[2], v[3]) }
        | flatten
        | branch { v -> 
            bests : v.toString().contains(".lca.best.json")
            consensuses : v.toString().contains(".lca.consensus.json")
            annotations : v.toString().contains(".read_annotation.json")
        }
        | set{ jsons }


    jsons.bests.toList()
        | concat(jsons.consensuses.toList())
        | concat(jsons.annotations.toList())
        | toList
        | make_tables
}

workflow {

new_dbs = make_dbs()
new_dbs.view()
mapping(new_dbs)

}

// End of file
            
    

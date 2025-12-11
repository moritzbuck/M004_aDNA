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




workflow mapping{

    main : 



    Channel.fromPath(file(params.bam_dir).resolve("*.bam")).map(x -> tuple( x.toString().split("/")[-1].split(".sorted")[0]  , x ) ).set { bams }
     
    bams.view()

    bams
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


mapping()

}

// End of file
            
    

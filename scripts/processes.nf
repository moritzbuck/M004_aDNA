process validate{
    label "low_cpu"
    input:
    val(check_dbs)
    output:
    val([])
    script:
    """
    """
}

process unpigz{
    label "low_cpu"
    cache 'true'
    input:
    tuple val(sampleID), path(fastq_pair)
    output:
    tuple val(sampleID), file("*_R[12].fastq")
    script:
    """
    unpigz -c ${fastq_pair[0]} > ${sampleID}_R1.fastq
    unpigz -c ${fastq_pair[1]} > ${sampleID}_R2.fastq
    """
}

process check_databases{
    label "no_cpu"
    cache 'false'
    input:
    path(fasta_db_dir)
    output :
    file("to_process_dbs.txt")
    script :
    """
    python ${params.script_folder}/check_databases.py --dir ${fasta_db_dir} --accession2taxid ${params.accession2taxid_path} --output to_process_dbs.txt --db_dir ${params.db_dir}
    """
}

process bt2_index{
    publishDir "${params.db_dir}", mode: 'copy', overwrite: true
    input:
    val(genome)
    output :
    file("${genome}*")
    script :
    """
    bowtie2-build --large-index --threads ${task.cpus} ${params.fasta_db_dir}/${genome}/${genome}*_genomic.fna ${genome} 
    """
}

process merge_fastq {
    input:
    tuple val(sampleID), file(fastq_pair)
    output:
    tuple val(sampleID), file("${sampleID}.assembled.fastq")
    script:
    """
    echo $sampleID
    pear -f ${fastq_pair[0]} -r ${fastq_pair[1]} -o  ${sampleID} --threads ${task.cpus}
    """
}

process remove_low_complex{
    label "low_cpu"
    input:
    tuple val(sampleID), file(assembled_fastq)
    output:
    tuple val(sampleID), file("${sampleID}.cleaned.fastq")
    script:
    """
    sga preprocess --min-length ${params.min_length} --dust-threshold ${params.dust_threshold} -p 0  --dust -o ${sampleID}.cleaned.fastq ${assembled_fastq}
    """
}

process sam2bam {
    input:
    tuple val(sampleID), val(db_id), file(sam_file)
    output:
    tuple val(sampleID), val(db_id), file("${sampleID}_x_${db_id}.sorted.bam")
    script:
    """
    samtools view -b ${sam_file} | samtools sort -n -@ ${task.cpus} -m 6G -o ${sampleID}_x_${db_id}.sorted.bam -
    """
} 

process bowtie2 {
    cache 'lenient'
    input:
    tuple val(sampleID), file(assembled_fastq), val(db)
    output:
    tuple val("${sampleID}"), val("${db_id}"), file("${sampleID}_x_${db_id}.sam")

    script:
    if (params.bowtie_preset == "slow") {
        bowtie_setting_string = "-t -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -k 10 "
    } else if (params.bowtie_preset == "less_slow") {
        // without the mismatch
        bowtie_setting_string = "-t -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -k 10"
    } else if (params.bowtie_preset == "new_holi") {
        bowtie_setting_string = "-t -k 1000"
    } else if (params.bowtie_preset == "bowtie2_default") {
        bowtie_setting_string = ""
    }

    db_id = file(db).getFileName()

    """
    bowtie2 -x ${db} -U ${assembled_fastq} ${bowtie_setting_string} --end-to-end --no-unal  --threads ${task.cpus}  -S ${sampleID}_x_${db_id}.sam
    """
}

process filter_sam {
    label "low_cpu"
    cache 'lenient'
    input:
    tuple val(sampleID), val(db_id), file(sam_file)
    output:
    tuple val(sampleID), val(db_id), file("${sampleID}_x_${db_id}.filtered")
    script:
    """
    #!/usr/bin/env python
    import os, sys, shutil
    from subprocess import call
    import re
    from tqdm import tqdm
    import gzip
    from ete3 import NCBITaxa
    from Bio import Entrez
    from speedict import Rdict, Options, AccessType
    from random import choice


    print("loading taxdb")
    db = Rdict(f"${params.accession2taxid_path}", access_type=AccessType.read_only())
    seq2tax = {}

    kranks = ['superkingdom',
    'kingdom',
    'phylum',
    'subphylum',
    'class',
    'subclass',
    'infraclass',
    'cohort',
    'order',
    'suborder',
    'infraorder',
    'superfamily',
    'family',
    'subfamily',
    'genus',
    'species']

    sranks = ['superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species']

    ncbi = NCBITaxa()
    _taxid2str = {}
    def taxid2str(taxid):
        taxid = int(taxid)
        if taxid not in _taxid2str:
            if ncbi.get_taxid_translator([taxid]):
                lineage = ncbi.get_lineage(taxid)
                taxa = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(lineage)
                common_names = ncbi.get_common_names(lineage)
                common_name = ";".join([common_names.get(l, "") for l in lineage if ranks[l] in kranks])
                ranks = {v : k for k,v in ranks.items() if v in kranks}
                strs = [ "" if k not in ranks else taxa[ranks[k]] for k in kranks]
                taxstr = ";".join(strs)
                _taxid2str[taxid] = (taxstr, common_name)
            else :
                _taxid2str[taxid] = (None, None)
        return _taxid2str[taxid]

    def match_counts(md):
        str_ = md.split(":")[-1]
        assert all([s in "^ATGC012345678N9" for s in str_])
        matches = sum([int(s) for s in re.split("[ATCGN^]", str_) if s])
        return matches

    def processing_sam(parsed_sam, identity) :
        parsed_sam = [k for k in parsed_sam if k['identity'] > identity]
        for p in tqdm(parsed_sam):
            taxid = db.get(str(p['target'].split(".")[0]), None)
            if taxid :
                if taxid not in seq2tax:
                    seq2tax[taxid] = taxid2str(taxid)
                taxstrs = seq2tax[taxid]
                p['taxid'] = taxid
                p['taxstr_sci'] = taxstrs[0]
                p['taxstr_com'] = taxstrs[1]
        return parsed_sam

    def parse_samline(l):
        saml = l.split()
        target = saml[2]
        read = saml[0]
        read_l = len(saml[9])
        md = [ll for ll in saml if ll.startswith("MD:")][0]
        matches = match_counts(md)
        return {'read'    : read,
                'target'  : target,
                'identity': matches/read_l}

    print("parsing ${sam_file}")
    with open("${sam_file}") as handle:
        parsed_sam = [parse_samline(l) for l in tqdm(handle) if not l.startswith("@")]
    print(f"filtering sam with identity > ${params.identity_cutoff}")
    processed_sam = processing_sam(parsed_sam, ${params.identity_cutoff})
    print(f"writing filtered sam to ${sampleID}_x_${db_id}.filtered")
    with open("${sampleID}_x_${db_id}.filtered", "w") as out:
        header = processed_sam[0].keys()
        out.write("\\t".join(header) + "\\n")
        out.writelines(["\\t".join([str(p[k]) for k in header]) + "\\n" for p in processed_sam])
    sys.exit(0)
    """
}

process merge_hits {
    label "no_cpu"
    input:
    tuple val(sampleID), val(dbs), file(sam_files)
    output:
    tuple val(sampleID), file("${sampleID}.hits")
    script:
    """
    head -n1 ${sam_files[0]} > ${sampleID}.hits
    tail -n +2 ${sam_files} >> ${sampleID}.hits
    """
}

process merge_bams {
    label "high_mem"
    input:
    tuple val(sampleID), val(dbs), file(bam_files)
    output:
    tuple val(sampleID), file("${sampleID}.merged.bam")
    script:
    """
    samtools merge -@ ${task.cpus} ${bam_files} -o ${sampleID}.merged.bam
    """
}

process filter_bam {
    input:
    tuple val(sampleID), file(bam_file)
    output:
    tuple val(sampleID), file("${sampleID}.filtered.bam")
    script:
    """
    filterBAM reassign --bam ${bam_file} -t ${task.cpus} -i 0 -A 92 -M 30G -m 5G -n 10 -s 0.75 -o ${sampleID}.filtered.bam
    """
}

process sort_bam {
    label "high_mem"
    input:
    tuple val(sampleID), file(bam_file)
    output:
    tuple val(sampleID), file("${sampleID}.sorted.bam")
    script:
    """
    samtools sort -@ ${task.cpus} -m 500M -o ${sampleID}.sorted.bam ${bam_file}
    """
}

process index_bam {
//    label "high_mem"
    input:
    tuple val(sampleID), file(bam_file)
    output:
    tuple val(sampleID), file("${sampleID}.sorted.bam")
    script:
    """
    samtools index ${bam_file}
    """
}

process bam2lca {
    label "no_cpu"
    input:
    tuple val(sampleID), file(bam_file)
    output:
    tuple val(sampleID), file("${sampleID}.lca.best.json"), file("${sampleID}.lca.consensus.json"), file("${sampleID}.read_annotation.json")
    script:
    """
    python ${params.script_folder}/bam2lca.py --bam ${bam_file} --accession2taxid ${params.accession2taxid_path} --min-identity ${params.identity_cutoff} --prefix ${sampleID}

    """
}

process make_tables {
    label "no_cpu"
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true
    input: tuple(file(lca_bests) ,file(lca_alls), file(read_annotations))
    output:
    tuple file("lca.best.tsv"), file("lca.consensus.tsv"), file("reads.fastq")
    script:
    """
    python ${params.script_folder}/make_tables.py --best ${lca_bests} --consensus ${lca_alls} --annotation ${read_annotations}
    """
}
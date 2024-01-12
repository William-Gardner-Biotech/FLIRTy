#!/opt/homebrew/Caskroom/miniconda/base/envs/aegis/bin/nextflow

nextflow.enable.dsl=2

workflow {

    ch_primers = Channel
    .fromPath( params.primers )
    .splitText()

    MATCH_AMP ( ch_primers )

    REMOVE_AMP_REGION (
        MATCH_AMP.out
    )

    BACK_BLAST ( 
        REMOVE_AMP_REGION.out 
    )

}

// process that will find the matching amplicon for the primer
process MATCH_AMP {

    tag "${primer_name}"

    publishDir "primer_bed", mode: 'symlink'

    input:
    val primer

    output:
    tuple val(primer), val( primer_name ), path( amplicon_bed )

    script:
    String[] parts = primer.split("\t")
    // naming the amplicons based on primer name for uniquenss as each primer shares an amp
    // and would cause collision
    primer_name = parts[3]
    amplicon_bed = primer_name + ".bed"
    """
    python3 ${projectDir}/bin/quick_match.py "${primer}" ${params.amplicons} ${amplicon_bed}
    """
}

// 
// process that will take the bam file and remove any reads found in the amplicon region
process REMOVE_AMP_REGION {
    
    tag "${primer_name}"

    publishDir "masked_beds", mode: 'Copy'

    input:
    tuple val(primer), val(primer_name), path( amplicon_bed )

    output:
    tuple val( primer ), val( primer_name ), path( masked_reads )

    script:
    masked_reads = primer_name + ".fasta"
    """
    # samtools faidx ${params.ref_genome}

    # cut -f1,2 ${params.ref_genome}.fai > genome.txt

    # Input BAM is from the sequencing reads we care about
    bedtools bamtobed -i ${params.inputBAM} > ${primer_name}_raw.bed

    bedtools subtract -a ${primer_name}_raw.bed -b ${amplicon_bed} > ${primer_name}_masked.bed

    bedtools getfasta -fi ${params.ref_genome} -bed ${primer_name}_masked.bed -fo ${masked_reads}

    # bedToBam -i ${primer_name}_masked.bed -g genome.txt > ${primer_name}.bam

    """
}

/*process REMAP {
    tag "${amplicon}"

    input:
    tuple val(primer), path(amplicon_bam), path(amplicon_fa)

    output:
    path "${amplicon_bam.getName().replace('.bam', '.fq')}"

    publishDir "back2reads", mode: 'symlink'

    script:
    println (primer)

    def amp = amplicon_bam.getName().replace(".bam", "")
    """
    samtools sort -o sorted_${amplicon_bam} ${amplicon_bam}

    samtools fastq sorted_${amplicon_bam} > ${amplicon_bam.getName().replace('.bam', '.fq')}
    """
}*/

// This process will now map the original primer to 
process BACK_BLAST {
    tag "${primer_name}"

    input:
    tuple val(primer), val( primer_name ), path( masked_reads )

    output:
    tuple val( primer_name ), path( primer_fasta ), path( blast_results ), path( rc_blast_results )

    publishDir "BLASTED", mode: 'Copy'

    script:
    primer_fasta = primer_name + "_primer.fa"
    blast_results = primer_name + "_results.txt"
    rc_blast_results = "RC_" + primer_name + "_results.txt"
    """
    # echo stored value into a bed file for use by bedtools
    echo "${primer}" > ${primer_name}_primer.bed

    # Extract the fasta from the bed
    bedtools getfasta -fi ${params.ref_genome} -bed ${primer_name}_primer.bed -fo ${primer_fasta}

    makeblastdb -in ${masked_reads} -dbtype nucl -out ${primer_name}_db

    blastn -query ${primer_fasta} -db ${primer_name}_db -out ${blast_results} -evalue 15 -task blastn-short

    # awk command to make a reverse complement and BLAST it as well
    awk '/^>/ {print ">RC_"\$0; getline seq} {print seq | "rev | tr ATCG TAGC"; close("rev | tr ATCG TAGC")}' ${primer_fasta} > RC_${primer_fasta}

    blastn -query RC_${primer_fasta} -db ${primer_name}_db -out ${rc_blast_results} -evalue 15 -task blastn-short

    """
}
// Add a slider feature that detects if there are results at x E-val then slides differently to find more or less
// Build a new blastdb for each different primer as diff amplicons remove diff reads
process BUILD_BLAST_DB {

    tag "${primer_name}"

    publishDir "BLASTDBs", mode: 'symlink', pattern: "${primer_name}_db.*"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    input:
    tuple val( primer ), val( primer_name ), path( masked_reads )

    output:
    tuple val( primer_name ), path( "${primer_name}_db.*" )

    script:
    """
    makeblastdb -in ${masked_reads} -dbtype nucl -out ${primer_name}_db
    """
}

process PERFORM_BLAST {

    tag "${primer_name}"

    publishDir "BLAST_RESULTS", mode: 'symlink'

    input:
    tuple val( blast_db_name ), path( blast_db_files )
    tuple val( primer_name ), path( query_primer )

    output:
    tuple val( primer_name ), path( blast_results )

    script:

    blast_results = primer_name + "_blast_results.txt"

    """
    blastn -query ${query_primer} -db ${blast_db_name} -out ${blast_results} -evalue 0.1
    """
}
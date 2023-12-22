#!/opt/homebrew/Caskroom/miniconda/base/envs/aegis/bin/nextflow

nextflow.enable.dsl=2

workflow {

    ch_primers = Channel
    .fromPath( params.primers)
    .splitText()

    MATCH_AMP ( ch_primers )

    REMOVE_AMP_REGION (
        MATCH_AMP.out
    )

}

// process that will find the matching amplicon for the primer
process MATCH_AMP {

    publishDir "primer_bed", mode: 'symlink'

    input:
    val primer

    output:
    tuple val(primer), path( amplicon )

    script:
    String[] parts = primer.split("\t")
    // naming the amplicons based on primer name for uniquenss as each primer shares an amp
    // and would cause collision
    amplicon = parts[3]

    """
    python3 ${projectDir}/bin/quick_match.py "${primer}" ${params.amplicons} ${amplicon}
    """
}

// process that will take the bam file and remove any reads found in the amplicon region
process REMOVE_AMP_REGION {
    
    tag "${amplicon}"

    publishDir "masked_beds", mode: 'Copy'

    input:
    tuple val(primer), path(amplicon)

    output:
    tuple val(primer), path("${amplicon}.sam")

    script:
    """
    samtools faidx ${params.ref_genome}

    cut -f1,2 ${params.ref_genome}.fai > genome.txt

    bedtools bamtobed -i ${params.inputBAM} > ${amplicon}_raw.bed

    bedtools subtract -a ${amplicon}_raw.bed -b ${amplicon} > ${amplicon}_masked.bed

    bedToBam -i ${amplicon}_masked.bed -g genome.txt > ${amplicon}.sam

    """
}




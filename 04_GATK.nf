#! /usr/bin/env nextflow
/* Name: GATK4 Pipeline
 * Auth: Jennifer Chang
 * Date: 2021/04/01
 * Desc:
 * Test:
 */

nextflow.enable.dsl=2

params.genome="00_Raw-Data/test-data/ref/*.fasta"
params.reads="00_Raw-Data/test-data/fastq/*_{R1,R2}.fastq.gz"
params.outdir="04_GATK"
params.reads_file=false
params.window=100000

// === Define Processes
process bwamem2_index {
    tag "${genome_fasta.simpleName}"
    label 'bwamem'
    executor 'slurm'
    clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(genome_fasta)

    output: // [ genome.fasta, [genome index files] ]
    tuple path("$genome_fasta"), path("${genome_fasta}*")

    script:
    """
    #! /usr/bin/env bash
    bwa-mem2 index ${genome_fasta}
    """
}

process bwamem2_mem {
    tag "${readname}"
    label 'bwamem'
    executor 'slurm'
    clusterOptions '-N 1 -n 16 -t 04:00:00 --account=isu_gif_vrsc'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(genome_fasta), path(genome_index_files), \
      val(readname), path(readpairs), val(i_readname)

    output:
    path("${i_readname}_mapped.bam")

    script:
    """
    #! /usr/bin/env bash
    bwa-mem2 mem -t 16 ${genome_fasta} ${readpairs} |\
     samtools view --threads 16 -bS - > ${i_readname}_mapped.bam
    """
}

process FastqToSam {
  tag "$readname"
  label 'gatk'
  executor 'slurm'
  clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

  publishDir "${params.outdir}"

  input:  // [readgroup, [left.fq.gz, right.fq.gz], increment_readgroup]
  tuple val(readname), path(readpairs), val(i_readname)

  output: // increment_readgroup.bam since one readgroup can have multiple lanes
  path("${i_readname}.bam")

  """
  #! /usr/bin/env bash
  gatk FastqToSam \
  --FASTQ ${readpairs.get(0)} \
  --FASTQ2 ${readpairs.get(1)} \
  --OUTPUT ${i_readname}.bam \
  --READ_GROUP_NAME ${readname} \
  --SAMPLE_NAME ${readname}_name \
  --LIBRARY_NAME ${readname}_lib \
  --PLATFORM ILLUMINA \
  --SEQUENCING_CENTER ISU
  """
}
// USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

process MarkIlluminaAdapters {
  tag "${bam.fileName}"
  label 'gatk'
  executor 'slurm'
  clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

  publishDir "${params.outdir}"

  input:  // reads.bam
  path(bam)

  output: // reads_marked.bam
  path("${bam.simpleName}_marked.bam")

  script:
  """
  #! /usr/bin/env bash
  gatk MarkIlluminaAdapters \
  --INPUT ${bam} \
  --OUTPUT ${bam.simpleName}_marked.bam \
  --METRICS ${bam.simpleName}_marked_metrics.txt
  """
}

process SamToFastq {
  tag "${bam.fileName}"
  label 'gatk'
  executor 'slurm'
  clusterOptions '-N 1 -n 16 -t 02:00:00 --account=isu_gif_vrsc'

  publishDir "${params.outdir}"

  input:  // reads.bam
  path(bam)

  output: // reads_interleaved.fq
  path("${bam.simpleName}_interleaved.fq")

  script:
  """
  #! /usr/bin/env bash
  gatk SamToFastq \
  --INPUT ${bam} \
  --FASTQ ${bam.simpleName}_interleaved.fq \
  --CLIPPING_ATTRIBUTE XT \
  --CLIPPING_ACTION 2 \
  --INTERLEAVE true \
  --INCLUDE_NON_PF_READS true
  """
}

process CreateSequenceDictionary {
  tag "$fasta"
  label 'gatk'
  publishDir "${params.outdir}"

  input:  // genome.fasta
  path(fasta)

  output: // genome.dict
  path("${fasta.simpleName}.dict")

  script:
  """
  #! /usr/bin/env bash
  gatk CreateSequenceDictionary \
  -R ${fasta} \
  -O ${fasta.simpleName}.dict
  """
}

process samtools_faidx {
  tag "$fasta"
  label 'samtools'
  publishDir "${params.outdir}"

  input:  // genome.fasta
  path(fasta)

  output:
  path("*.fai")

  """
  #! /usr/bin/env bash
  samtools faidx $fasta
  """
}

process MergeBamAlignment {
  tag "$i_readname"
  label 'gatk'
  publishDir "${params.outdir}"

  input:  // [readgroup, unmapped reads, mapped reads]
  tuple val(i_readname), path(read_unmapped), path(read_mapped), path(genome_fasta), path(genome_dict) //, path(genome_fai)
  //, path(genome_index), path(genome_fai), path(genome_dict)

  output: // merged bam and bai files
  tuple path("${i_readname}_merged.bam"), path("${i_readname}_merged.bai")

  script:
  """
  #! /usr/bin/env bash
  gatk MergeBamAlignment \
  --REFERENCE_SEQUENCE $genome_fasta \
  --UNMAPPED_BAM ${read_unmapped} \
  --ALIGNED_BAM ${read_mapped} \
  --OUTPUT ${i_readname}_merged.bam \
  --CREATE_INDEX true \
  --ADD_MATE_CIGAR true \
  --CLIP_ADAPTERS false \
  --CLIP_OVERLAPPING_READS true \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --ATTRIBUTES_TO_RETAIN XS
  """
}

process makewindows {
  tag "$fasta"
  label 'samtools'
  publishDir "${params.outdir}"

  input:  // genome.fasta
  path fasta

  output:
  path("${fasta.simpleName}_coords.bed")

  """
  #! /usr/bin/env bash
  samtools faidx $fasta
  awk -F'\t' '{print \$1"\t"\$2}' ${fasta}.fai > genome_length.txt
  bedtools makewindows -w $params.window -g genome_length.txt |\
    awk '{print \$1"\t"\$2+1"\t"\$3}' |\
    sed \$'s/\t/:/1' |\
    sed \$'s/\t/-/g' > ${fasta.simpleName}_coords.bed
  """
}

// === Main Workflow
workflow {
  genome_ch = channel.fromPath(params.genome, checkIfExists:true)
  if (params.reads) {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
  } else {
    reads_ch = channel.fromPath(params.reads_file, checkIfExists:true) |
    splitCsv(sep:'\t') |
    map { n -> [ n.get(0), [n.get(1), n.get(2)]] }
  }

  // == Since one sample may be run on multiple lanes
  i = 1
  ireads_ch = reads_ch | map { n -> [n.get(0), n.get(1), "${i++}_"+n.get(0)] }

  genome_ch | bwamem2_index | combine(ireads_ch) | bwamem2_mem
  ireads_ch.take(3) | FastqToSam | MarkIlluminaAdapters | SamToFastq

  index_pattern = ~/^\d+_/
  mapped_ch = bwamem2_mem.out |
    map { n -> [n.simpleName.replaceFirst("_mapped",""), n] }

  // Might be FastqToSam.out
  unmapped_ch = MarkIlluminaAdapters.out |
    map { n -> [n.simpleName.replaceFirst("_marked",""), n] }

  genome_ch | ( CreateSequenceDictionary & samtools_faidx & makewindows )
  both_ch = unmapped_ch | join(mapped_ch) | combine(genome_ch) | combine(CreateSequenceDictionary.out) | MergeBamAlignment

  allbam_ch = MergeBamAlignment.out | map { n -> n.get(0)} | collect
  allbai_ch = MergeBamAlignment.out | map { n -> n.get(1)} | collect
  allbam_ch | combine(allbai_ch) | view

/*
  makewindows.out | splitText( by:1) |
      map { n -> n.replaceFirst("\n","") } |
      combine(MergeBamAlignment.out) | view
*/

  
}

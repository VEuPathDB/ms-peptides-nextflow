#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { indexResults as indexProteins } from './modules/utils.nf'
include { indexResults as indexGenome } from './modules/utils.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  if(params.inputDirectory) {
    samples = channel.fromPath( params.inputDirectory + "/*" )
      .map { file ->
        def basename = file.getBaseName()
        return tuple(basename, file)
      }
  }
  else {
    error "Missing params.inputDirectory"
  }

  res = massSpecPeptides(samples);

  indexProteins(res.protein_gff.collectFile(), params.proteinGffFileName)
  indexGenome(res.genome_gff.collectFile(), params.genomeGffFileName)
}

process massSpecPeptides {
  container 'veupathdb/bioperl:1.0.0'

  input:
  tuple val(sampleName), path(sampleFile)

  output:
  path 'peptides_protein_align.gff', emit: protein_gff
  path 'peptides_genome_align.gff', emit: genome_gff

  script:
  """
  massSpecPeptides.pl \
   --sampleFile $sampleFile \
   --outputProteinGffFile peptides_protein_align.gff \
   --outputGenomicGffFile peptides_genome_align.gff \
   --proteinFastaFile $params.proteinFastaFile \
   --recordMinPeptidePct 50 \
   --inputAnnotationGff $params.annotationGff \
   --sampleName $sampleName \
   --runName ${workflow.runName}
  """
}



#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(params.inputDirectory) {
  samples = Channel.fromPath( params.inputDirectory + "/*" )
    .map { file ->
      def basename = file.getBaseName()
      return tuple(basename, file)
    }
}
else {
  throw new Exception("Missing params.inputDirectory")
}

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  massSpecPeptides(samples);
}

process massSpecPeptides {
  container = 'bioperl/bioperl:stable'

  input:
  tuple val(sampleName), path(sampleFile)

  script:
  """
  massSpecPeptides.pl \
   --sampleFile $sampleFile \
   --outputProteinGffFile peptides_protein_align.gff \
   --outputGenomicGffFile peptides_genome_align.gff \
   --proteinFastaFile $params.proteinFastaFile \
   --recordMinPeptidePct 50 \
   --inputAnnotationGff $params.annotationGff
   --sampleName $sampleName
  """








}

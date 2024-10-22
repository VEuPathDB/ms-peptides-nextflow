#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(params.inputDirectory) {
  samples = Channel.fromPath( params.inputDirectory + "/life_cycle_nature01107-Gmt.tab" )
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
    path sample

  script:
  """
  massSpecPeptides.pl \
   --sampleFile $sample \
   --outputFile temp \
   --proteinFastaFile $params.proteinFastaFile \
   --recordMinPeptidePct 50 \
   --outputProteinGffFile peptides.gff
  """
}

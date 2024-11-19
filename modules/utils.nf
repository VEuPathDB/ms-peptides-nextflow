#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process indexResults {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'
  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gff.gz'
    path '*.tbi'

  script:
  """
  sort -k1,1 -k4,4n $gff > $outputFileName
  bgzip $outputFileName
  tabix -p gff ${outputFileName}.gz
  """
}

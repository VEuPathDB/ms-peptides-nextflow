# ms-peptides-nextflow

Nextflow pipeline that maps legacy mass spectrometry peptide evidence onto a current
proteome and genome annotation.

## Overview

This pipeline is part of VEuPathDB's genomic data workflows. Genome annotations are
periodically revised, but the mass spectrometry (MS) peptide identifications that
support gene models are often generated against older, now-stale versions of the
proteome. This pipeline re-anchors that historical MS evidence onto the current
annotation: for each sample's search results, it locates every identified peptide
(and any annotated post-translational modification residues) within the current
protein sequences, then projects those peptide coordinates onto the genome using the
transcript's exon/CDS structure from the current annotation GFF. The result is a set
of protein- and genome-coordinate GFF3 tracks of MS peptide support, ready to be
indexed and loaded into a genome browser.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine: Docker (default) or Singularity/Apptainer

Containers used:
- `veupathdb/bioperl:1.0.0` — peptide-to-protein/genome mapping
- `biocontainers/tabix:v1.9-11-deb_cv1` — sorting, bgzip, and tabix indexing of output GFFs

## Usage

```
nextflow run VEuPathDB/ms-peptides-nextflow -r main -resume \
  --inputDirectory /path/to/sample_files \
  --proteinFastaFile /path/to/AnnotatedProteins.fsa \
  --annotationGff /path/to/annotation.gff \
  --outputDir /path/to/output \
  -C conf/docker.config
```

To run under Singularity/Apptainer (e.g. on an HPC cluster via LSF):

```
nextflow run VEuPathDB/ms-peptides-nextflow -r main -resume \
  --inputDirectory /path/to/sample_files \
  --proteinFastaFile /path/to/AnnotatedProteins.fsa \
  --annotationGff /path/to/annotation.gff \
  --outputDir /path/to/output \
  -C conf/lsf.config
```

There is a single, unnamed entry point — no `-entry` flag is needed. Each file in
`inputDirectory` is treated as one sample's MS search results, and samples are
processed as separate parallel tasks (the pipeline runs one task at a time by
default, via `process.maxForks = 1`).

### Sample file format

Each sample file (see `bin/massSpecPeptides.pl`) is a flat text export of an MS
protein/peptide/residue search result: a `#`-prefixed line starts a protein record
(source ID, description, molecular weight, pI, score, coverage, sequence/spectrum
counts), a `## start` line begins a peptide block for that protein (peptide
coordinates, observed/expected/calculated mass, sequence, ions score, spectrum
count), and a `## relative_position` line begins a modified-residue block for that
peptide (position, modification type).

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `inputDirectory` | `data/test` | Directory of per-sample MS search result files to process |
| `proteinFastaFile` | `data/AnnotatedProteins.fsa` | Current proteome FASTA that peptides are mapped onto |
| `annotationGff` | `data/pfal3D7.gff` | Current genome annotation GFF3 (CDS/exon structure) used to project peptides onto genomic coordinates |
| `proteinGffFileName` | `ms_peptides_protein_align.gff` | Name of the merged, protein-coordinate output GFF |
| `genomeGffFileName` | `ms_peptides_genome_align.gff` | Name of the merged, genome-coordinate output GFF |
| `outputDir` | `$launchDir/output` | Directory where output files are published |

The minimum percentage of a protein record's peptides that must match the target
sequence for a peptide-to-protein assignment to be accepted is fixed in `main.nf`
(`recordMinPeptidePct 50`).

## Output

Written to `outputDir`, both bgzip-compressed and tabix-indexed:

- **`ms_peptides_protein_align.gff.gz`** (+ `.tbi`) — `ms_peptide` and `ms_residue`
  features in protein-relative coordinates, one track entry per identified peptide
  (and its modified residues) per sample, tagged with ions score, spectrum count,
  peptide sequence, and sample name.
- **`ms_peptides_genome_align.gff.gz`** (+ `.tbi`) — the same peptide evidence
  projected onto genomic coordinates as `ms_peptide` parent features with `align`
  child features spanning exon boundaries, suitable for display alongside the gene
  model in a genome browser.

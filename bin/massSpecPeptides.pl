#!/usr/bin/perl

use strict;

use Getopt::Long;

use Bio::SeqIO;

use Bio::Tools::GFF;

use Bio::Coordinate::GeneMapper;

use Data::Dumper;

my ($sampleFile, $genomicGff, $proteinFastaFile, $minPeptidePct, $proteinGff, $inputAnnotationGff, $sampleName);

&GetOptions ("sampleFile=s" => \$sampleFile,
             "proteinFastaFile=s" => \$proteinFastaFile,
             "recordMinPeptidePct=i" => \$minPeptidePct,
             "outputProteinGffFile=s" => \$proteinGff,
             "outputGenomicGffFile=s" => \$genomicGff,
             "inputAnnotationGff=s" => \$inputAnnotationGff,
             "sampleName=s" => \$sampleName
    );


$minPeptidePct = $minPeptidePct ? $minPeptidePct : 50;

my $GFF_SOURCE = "veupathdb";

my ($transcriptLocations, $proteinToTranscriptMap) = &makeProteinGenomeCoordinatesHash($inputAnnotationGff);

open (my $proteinGffFh, ">$proteinGff") or die "Cannot open $proteinGff file for writing: $!";

open (my $genomicGffFh, ">$genomicGff") or die "Cannot open $genomicGff file for writing: $!";

my ($recordSet, $peptideSequencesToRecords) = &parseSampleFile($sampleFile);

my $seqio = Bio::SeqIO->new(-file => $proteinFastaFile, -format => 'fasta');

my $proteinCount;
my $seenGenomicCoords = {}; # keep track of unique genmic coords for a protein

while (my $seq = $seqio->next_seq) {

  my $msRecordWithPeptideLocations = &searchProteinSeq($seq, $recordSet, $peptideSequencesToRecords);

  &writeProteinGffOutput($proteinGffFh, $msRecordWithPeptideLocations, $seq, $sampleName);

  my $proteinId = $seq->id();
  my $transcriptId = $proteinToTranscriptMap->{$proteinId};
  my $locations = $transcriptLocations->{$transcriptId};

  my ($genomicSequenceSourceId, $mapper) = &getProteinToGenomicCoordMapper($locations);

  &writeGenomeGffOutput($genomicGffFh, $msRecordWithPeptideLocations, $genomicSequenceSourceId, $mapper, $sampleName);

  if($proteinCount++ % 500 == 0) {
    print STDERR "Processed $proteinCount Proteins", "\n";
  }

}

sub makeProteinGenomeCoordinatesHash {
  my ($gff) = @_;

  my %transcriptLocations;
  my %proteinToTranscriptMap;

  my $parser = Bio::Tools::GFF->new(-gff_version => 3,
                                    -file        => $gff);

  while(my $feature = $parser->next_feature()) {
    my $primaryTag = $feature->primary_tag();

    next unless($primaryTag eq 'CDS' || $primaryTag eq 'exon');

    # parent for exon and cds rows is the transcript
    # which we can use to join them
    my ($parentId) = $feature->get_tag_values("Parent");
    my $start = $feature->start();
    my $end = $feature->end();
    my $strand = $feature->strand();
    my $seqId = $feature->seq_id();

    my $loc = [$start, $end, $strand];

    if($primaryTag eq 'CDS') {
      my ($proteinId) = $feature->get_tag_values("protein_source_id");

      $proteinToTranscriptMap{$proteinId} = $parentId;

      push @{$transcriptLocations{$parentId}->{cds}}, $loc;
      $transcriptLocations{$parentId}->{seq_id} = $seqId;
    }

    if($primaryTag eq 'exon') {
      push @{$transcriptLocations{$parentId}->{exon}}, $loc;
    }
  }

  return \%transcriptLocations, \%proteinToTranscriptMap;
}

sub writeGenomeGffOutput {
  my ($genomicGffFh, $record, $genomicSequenceSourceId, $mapper, $sampleName) = @_;

  foreach my $peptide (@{$record->{peptides}}) {
    my $peptideSequence = $peptide->get("sequence");
    my $ionsScore = $peptide->get("ions_score");
    my $spectrumCount = $peptide->get("spectrum_count");

    my $peptideLocations = $record->{peptideLocations}->{$peptideSequence};
    foreach my $location (@{$peptideLocations}) {

      my $peptideCoords = Bio::Location::Simple->new (-start => $location->[0],
                                                      -end   => $location->[1]
          );

      my $map = $mapper->map($peptideCoords);
      return undef if ! $map;

      foreach (sort { $a->start <=> $b->start } $map->each_Location ) {

        my $genomicPeptide = new Bio::SeqFeature::Generic(
          -start      => $_->start,
          -end        => $_->end,
          -strand     => $_->strand,
          -primary    => 'ms_peptide',
          -source_tag => $GFF_SOURCE,
          -seq_id     => $genomicSequenceSourceId,
          -tag        => {
            ions_score => $ionsScore,
            spectrum_count => $spectrumCount,
            sample_name => $sampleName,
            peptide => $peptideSequence,
          });

        $genomicPeptide->gff_format(Bio::Tools::GFF->new(-gff_version => 3));

        print $genomicGffFh $genomicPeptide->gff_string(), "\n";
      }
    }
  }

}

sub writeProteinGffOutput {
  my ($proteinGffFh, $record, $seq, $sampleName) = @_;

  my $protein = $seq->id();

  foreach my $peptide (@{$record->{peptides}}) {
    my $peptideSequence = $peptide->get("sequence");

    my $ionsScore = $peptide->get("ions_score");
    my $spectrumCount = $peptide->get("spectrum_count");

    my $peptideLocations = $record->{peptideLocations}->{$peptideSequence};
    foreach my $location (@{$peptideLocations}) {
      my $pepStart = $location->[0];
      my $pepEnd = $location->[1];
      my $peptideId = "${sampleName}_${protein}_${pepStart}_${pepEnd}";

      my $peptideFeature = new Bio::SeqFeature::Generic(
        -start      => $pepStart,
        -end        => $pepEnd,
        -strand     => ".",
        -primary    => 'ms_peptide',
        -source_tag => $GFF_SOURCE,
        -seq_id     => $protein,
        -tag        => {
          ions_score => $ionsScore,
          spectrum_count => $spectrumCount,
          peptide => $peptideSequence,
          sample_name => $sampleName,
          ID => $peptideId
        });

      $peptideFeature->gff_format(Bio::Tools::GFF->new(-gff_version => 3));

      print $proteinGffFh $peptideFeature->gff_string(), "\n";

      foreach my $residue (@{$peptide->{residues}}) {
        my $relativePosition = $residue->get("relative_position");
        my $modificationType = $residue->get("modification_type");

        my $residueLocation = $pepStart + $relativePosition;

        my $residue = substr($peptideSequence, $relativePosition, 1);

        my $residueId =  "${sampleName}_${protein}_${pepStart}_${pepEnd}_r${residueLocation}";

        my $residueFeature = new Bio::SeqFeature::Generic(
          -start      => $residueLocation,
          -end        => $residueLocation,
          -strand     => ".",
          -primary    => 'ms_residue',
          -source_tag => $GFF_SOURCE,
          -seq_id     => $protein,
          -tag        => {
            sample_name => $sampleName,
            Parent => $peptideId,
            ID => $residueId,
            residue => $residue,
            modification_type => $modificationType,
            relative_position => $relativePosition
          });

        $residueFeature->gff_format(Bio::Tools::GFF->new(-gff_version => 3));

        print $proteinGffFh $residueFeature->gff_string(), "\n";
      }
    }
  }
}


sub searchProteinSeq {
  my ($seq, $recordSet, $peptideSequencesToRecords) = @_;

  my %proteinAtts;
  my $desc = $seq->desc();
  foreach my $pair (split(/\s/, $desc)) {
    my ($k, $v) = split(/=/, $pair);
    $proteinAtts{$k} = $v;
  }

  my $geneId = $proteinAtts{gene};
  my $seqString = $seq->seq();

  my ($matchingRecords, $peptideLocations, $bestMatch);

  # NOTE: if the id from the record matches this gene then we only need to deal
  # with the peptides from the record.
  if(my $records = $recordSet->{$geneId}) {
    my $peptides = &getPeptidesFromRecords($records);
    $peptideLocations = &mapPeptidesToProtein($seqString, $peptides);
    $matchingRecords = $records;

    $bestMatch = &findBestRecord($matchingRecords, $peptideLocations, $minPeptidePct);
  }

  # NOTE: if we were able to find a match using the id mapping then we're done
  # otherwise we need to try to map all peptides to this protein
  unless($bestMatch) {

    $peptideLocations = &mapPeptidesToProtein($seqString, $peptideSequencesToRecords);
    $matchingRecords = &getRecordsFromPeptideSubset($peptideLocations, $peptideSequencesToRecords);

    # NOTE: here we require 100% peptide matches
    $bestMatch = &findBestRecord($matchingRecords, $peptideLocations, 100);
  }

  return &cloneAndAddPeptideLocations($bestMatch, $peptideLocations);
}


sub cloneAndAddPeptideLocations {
  my ($record, $peptideLocations) = @_;

  return unless($record);

  # make a copy of the MSRecord
  my $clone = $record->clone();

  my $peptideCount = scalar @{$clone->{peptides}};

  foreach my $peptide (@{$clone->{peptides}}) {
    my $peptideSequence = $peptide->get("sequence");

    foreach my $location (@{$peptideLocations->{$peptideSequence}}) {
      push @{$clone->{peptideLocations}->{$peptideSequence}}, $location;
    }
  }
  return $clone;
}


sub findBestRecord {
  my ($matchingRecords, $peptideLocations, $minPeptidePct) = @_;

  my @results;

  foreach my $record (@$matchingRecords) {
    my $count = 0;
    my $peptideCount = scalar @{$record->{peptides}};

    foreach my $peptide (@{$record->{peptides}}) {
      my $peptideSequence = $peptide->get("sequence");

      $count++ if($peptideLocations->{$peptideSequence});
    }

    my $matchPct = ($count / $peptideCount) * 100;

    if($matchPct >= $minPeptidePct) {

      # NOTE: Use count here, not the pct
      push @results, [$count, $record];
    }
  }

  return unless scalar @results > 0;

  my @sorted = sort { $b->[0] <=> $a->[0]} @results;


  # NOTE:  we are returning the record which has the highest
  # number of peptides matching above our minPeptidePct threshold
  return $sorted[0]->[1];
}


sub getRecordsFromPeptideSubset {
  my ($peptideLocations, $peptideSequencesToRecords) = @_;

  my (@res, %seen);

  foreach my $peptide (keys %$peptideLocations) {
    foreach my $record (@{$peptideSequencesToRecords->{$peptide}}) {

      my $uid = $record->{_unique_id};
      
      # NOTE: multiple peptides will return the same record reference
      # we only want it one time in our array
      push @res, $record unless($seen{$uid});

      $seen{$uid} = 1;

    }
  }
  return \@res;
}

sub mapPeptidesToProtein {
  my ($seqString, $peptides) = @_;

  my %res;

  foreach my $peptide (keys %$peptides) {
    my $pattern = quotemeta($peptide);

    my $pepLength = length $peptide;

    while ($seqString =~ /$pattern/g) {
      my $end = pos($seqString) + 0;
      my $start = $end - $pepLength + 1;

      push @{$res{$peptide}}, [$start, $end];
    }
  }

  return \%res;
}

sub getPeptidesFromRecords {
  my ($records) = @_;

  my %peptides;

  foreach my $record (@$records) {
    foreach my $peptide (@{$record->{peptides}}) {

      my $peptideSequence = $peptide->get("sequence");
      $peptides{$peptideSequence}++;
    }
  }
  return \%peptides;
}

sub parseSampleFile {
  my ($sampleFile) = @_;

  warn "now reading input file $sampleFile \n";

  open(F, $sampleFile) or die "Could not open $sampleFile: $!\n";

  my $recordSet = {};
  my ($record, $peptide, $residue);

  my $recordUniqueId = 1;

  my $state;

  # NOTE: This Variable is mapping unique peptide sequences
  # to all of the records which contain it
  my %peptideSequencesToRecords;

  while (<F>) {
    chomp;

    next if /^\s*$/;
    if (/^# /) {
      $state = 'record';
      next;
    } elsif(/^## start/) {
      $state = 'peptide';
      undef $peptide;
      next;
    } elsif(/^## relative_position/) {
      $state = 'residue';
      next;
    } elsif($state eq 'record') {
      $record = MSRecord->new($_);
      $record->{_unique_id} = $recordUniqueId++;
      $record->set("file", $sampleFile);
      push @{$recordSet->{$record->get("sourceId")}}, $record;
    } elsif($state eq 'peptide') {
      $peptide = MSPeptide->new($_);
      my $peptideSequence = $peptide->get("sequence");
      if(!$peptideSequence || $peptideSequence =~ /\d/) {
        next;# TODO:  why not die here instead??
      }
      $record->addPeptide($peptide);
      # NOTE:  I'm allowing for multiple records having the same id here
      # This shouldnt' happen but some of the input files used alias mappings so it might?
      push @{$peptideSequencesToRecords{$peptideSequence}}, $record;
    } elsif($state eq 'residue') {
      $residue = MSResidue->new($_);
      $peptide->addResidue($residue);
    } else {
      die "Something wrong in the file on this line: $_\n";
    }
  }
  close F;

  return $recordSet, \%peptideSequencesToRecords;
}


sub getProteinToGenomicCoordMapper {
  my ($locations) = @_;

  my ($sequenceSourceId, $exonLocs, $cdsRange) = &getExonLocsAndCdsRangeFromLocations($locations);

  my $mapper = Bio::Coordinate::GeneMapper->new(
    -in    => 'peptide',
    -out   => 'chr',
    -exons => $exonLocs,
    -cds => $cdsRange
      );

  return $sequenceSourceId, $mapper;
}


sub getExonLocsAndCdsRangeFromLocations {
  my ($locations) = @_;

  my $naSequenceSourceId = $locations->{seq_id};

  my ($minCds, $maxCds, $strand);

  my @exonLocs;

  foreach my $exon (@{$locations->{exon}}) {
    my $exonLoc = Bio::Location::Simple->new( -seq_id => $naSequenceSourceId,
                                              -start => $exon->[0],
                                              -end => $exon->[1],
                                              -strand => $exon->[2]);

    push @exonLocs, $exonLoc;
  }

  foreach my $cds (@{$locations->{cds}}) {
    # init min and max w/ first value
    $minCds = $cds->[0] unless($minCds);
    $maxCds = $cds->[0] unless($maxCds);

    my ($min, $max) = sort {$a <=> $b} ($cds->[0], $cds->[1]);

    $minCds = $min if($min < $minCds);
    $maxCds = $max if($max > $maxCds);

    #used for exon and cds
    $strand = $cds->[2]

  }

  my $cdsRange = Bio::Location::Simple->new( -seq_id => $naSequenceSourceId,
                                             -start => $minCds,
                                             -end => $maxCds,
                                             -strand => $strand);


  return($naSequenceSourceId, \@exonLocs, $cdsRange);
}


1;

#------------------------------------------------------------
# Helper packages start here
#------------------------------------------------------------
package MSBase;
use strict;

sub new {
  my ($class, $ln) = @_;

  my $self = bless({}, $class);

  $self->initRecord($ln);

  return $self;
}

sub get {
  my ($self, $k) = @_;

  return $self->{$k};
}

sub set {
  my ($self, $k, $v) = @_;

  return $self->{$k} = $v;
}


1;

#------------------------------------------------------------
# Helper for MS Protein Record
#------------------------------------------------------------

package MSRecord;
use base qw/MSBase/;
use strict;

sub addPeptide {
  my ($self, $pep) = @_;

  push @{$self->{peptides}}, $pep;
}


sub initRecord {
  my ($self, $ln) = @_;

  ( $self->{sourceId},
    $self->{description},
    $self->{seqMolWt},
    $self->{seqPI},
    $self->{score},
    $self->{percentCoverage},
    $self->{sequenceCount},
    $self->{spectrumCount},
    $self->{sourcefile},
  ) = split "\t", $ln;

  return $self;
}


# NOTE:  the peptide objects will NOT be copies. they will be the existing references
sub clone {
    my $self = shift;
    return bless { %$self }, ref $self;
}

1;

#------------------------------------------------------------
# Helper for MS Peptide
#------------------------------------------------------------

package  MSPeptide;
use base qw/MSBase/;
use strict;

sub addResidue {
  my ($self, $residue) = @_;

  push @{$self->{residues}}, $residue;
}

sub initRecord {
  my ($self, $ln) = @_;

  ( $self->{start},
    $self->{end},

    $self->{observed},
    $self->{mr_expect},
    $self->{mr_calc},
    $self->{delta},
    $self->{miss},
    $self->{sequence},
    $self->{modification},
    $self->{query},
    $self->{hit},

    $self->{ions_score},

    $self->{spectrum_count}
  ) = split "\t", $ln;

  $self->{residues} = [];

  die "missing pep spectrum_count" unless($self->{spectrum_count});
  return $self;
}

1;

#------------------------------------------------------------
# Helper for MS Residue
#------------------------------------------------------------

package  MSResidue;
use base qw/MSBase/;
use strict;

sub initRecord {
  my ($self, $ln) = @_;

  ( $self->{relative_position},
    $self->{order},
    $self->{description},
    $self->{modification_type}
  ) = split "\t", $ln;

  return $self;
}

1;

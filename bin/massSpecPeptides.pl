#!/usr/bin/perl

use strict;

use Getopt::Long;

use Bio::SeqIO;

use Bio::Tools::GFF;

use Data::Dumper;

my ($sampleFile, $outputFile, $proteinFastaFile, $minPeptidePct, $proteinGff);

&GetOptions ("sampleFile=s" => \$sampleFile,
             "proteinFastaFile=s" => \$proteinFastaFile,
             "outputFile=s" => \$outputFile,
             "recordMinPeptidePct=i" => \$minPeptidePct,
             "outputProteinGffFile=s" => \$proteinGff,
    );


$minPeptidePct = $minPeptidePct ? $minPeptidePct : 50;

my $GFF_SOURCE = "veupathdb";

open (my $proteinGffFh, ">$proteinGff") or die "Cannot open $proteinGff file for writing: $!";

my ($recordSet, $peptideSequencesToRecords) = &parseSampleFile($sampleFile);

my $seqio = Bio::SeqIO->new(-file => $proteinFastaFile, -format => 'fasta');

my $proteinCount;

while (my $seq = $seqio->next_seq) {

  my $msRecordWithPeptideLocations = &searchProteinSeq($seq, $recordSet, $peptideSequencesToRecords);



  &writeOutput($proteinGffFh, $msRecordWithPeptideLocations, $seq);
  #&genomeGff($msRecordWithPeptideLocations)

  #&writeMSTabOutput($msRecordWithPeptideLocations)
  #&writeMSResiduesTabOutput($msRecordWithPeptideLocations)

  if($proteinCount++ % 500 == 0) {
    print STDERR "Processed $proteinCount Proteins", "\n";
  }

}

sub writeOutput {
  my ($proteinGffFh, $record, $seq) = @_;

  foreach my $peptide (@{$record->{peptides}}) {
    my $peptideSequence = $peptide->get("sequence");

    my $ionsScore = $peptide->get("ions_score");
    my $spectrumCount = $peptide->get("spectrum_count");

    my $peptideLocations = $record->{peptideLocations}->{$peptideSequence};
    foreach my $location (@{$peptideLocations}) {

      my $peptide = new Bio::SeqFeature::Generic(
        -start      => $location->[0],
        -end        => $location->[1],
        -strand     => ".",
        -primary    => 'ms_peptide',
        -source_tag => $GFF_SOURCE,
        -seq_id     => $seq->id(),
        -tag        => {
          ions_score => $ionsScore,
          spectrum_count => $spectrumCount,
        });

      $peptide->gff_format(Bio::Tools::GFF->new(-gff_version => 3));

      print $proteinGffFh $peptide->gff_string(), "\n";
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
    die "";
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

    my $pepLength = scalar $peptide;

    while ($seqString =~ /$pattern/g) {
      my $start = pos($seqString) + 1;
      my $end = $start + $pepLength + 0;

      push @{$res{$peptide}}, [$start, $end];
    }
  }

  return \%res;
}

sub getPeptidesFromRecords {
  my ($records) = @_;

  my %peptides;

  foreach my $record (@$records) {
    foreach my $peptide ({$record->{peptides}}) {

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

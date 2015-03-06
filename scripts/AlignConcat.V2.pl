#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;
use Getopt::Std;

#OPTIONS
my %opts;
getopts('f:r:', \%opts);
#DEFAULT OPTIONS
my $format = $opts{f} || 'clustalw';
my $ref_genome = $opts{r} || 1; #REFERENCE GENOME IS 1 BY DEFAULT

#REFERENCE NUMBER FOR EACH ORTHOLOG FASTA FILE
my $ref = 0;
my (%gid_seen, %results, @ref_orf);

foreach my $fh (@ARGV)
  {
    #INPUT FILE
    my $in = Bio::AlignIO->new(-file => $fh, -format => $format);
    my $aln = $in->next_aln;

    #EXTRACT DATA
    foreach my $seq ($aln->each_seq() )
      {
	my $id = $seq->id();
	if ($id =~/^gid_(\d+)\|(\S+)$/)
	  {
	    my ($gid, $locus_tag) = ($1, $2);
	    if (!(exists $gid_seen{$gid}))
	      {
		$gid_seen{$gid}++;
	      }
	    #STORE SEQUENCE DATA
	    $results{"$ref-$gid"} = $seq->seq;
	    #CREATE A REFERENCE ORDER FOR ORTHOLOG FILES
	    #ARRAY IS NOT PROPERLY ORDERED AT THIS POINT
	    if ($gid == $ref_genome)
	      {
		my $orf = {
                       'locus_tag' => $locus_tag,
                       'ref' => $ref
                      };
            push(@ref_orf, $orf);
	      }
	  }
	else
	  {
	    warn "Wrong id!\n";
	  }
      }
    $ref++;
  }

#PUT THE ORTHOLOGS IN ALPHABETICAL ORDER
#THIS IS ASSUMED THAT THE ORTHOLOGS ARE NAMED BY THE ORDER IN WHICH THEY START AND STOP
my @ordered_orf = sort { $a->{locus_tag} cmp $b->{locus_tag} } @ref_orf;

my %aligned;
#CONCATENTATE BY REFERENCE ORDER
foreach my $orf (@ordered_orf)
  {
    #GET REFERENCE
    my $ref = $orf->{ref};
    foreach my $seen (keys %gid_seen)
      {
	$aligned{$seen} .= $results{"$ref-$seen"};
      }
  }


#PRINT OUT ALIGNED FASTA
foreach my $seen (sort {$a <=> $b} keys %gid_seen)
  {
    print ">gid_", $seen, "\n", $aligned{$seen}, "\n";
  }

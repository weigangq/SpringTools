#!/usr/bin/env perl

use strict;
use warnings;
use DBI;
use Bio::SeqIO;

die "Usage: $0 <fasta file>\n" unless @ARGV > 0;

my $filename = shift @ARGV;

#########################
# PARSE UNKNOWN ORFS
#########################
my $in = Bio::SeqIO->new( -file => $filename, -format => 'fasta');

while(my $seq = $in->next_seq())
{
    my $id = $seq->id();
    #Ex: >orf00085|scf7180000000018|gid_13
    my ($locus_tag, $contig_id, $gid, $start, $stop) = split('\|', $id);
    # Insure that the start position is strictly before the stop.
    my $strand = "t";
    if ($start > $stop)
    {
	my $temp = $start;
	$start = $stop;
	$stop = $temp;
	$strand = "f";
    }
    $locus_tag =~ s/orf/$contig_id\_/;
    $gid =~ s/gid_//;
    print $gid, "\t", $locus_tag, "\t", $start, "\t", $stop, "\t", $strand, "\t", $seq->seq(), "\n";
}

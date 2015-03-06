#!/usr/bin/env perl

use strict;
use warnings;
use DBI;

die "Usage: $0 <Reciprocol Top Hits File>\n" unless @ARGV > 0;

my $filename = shift @ARGV;

open (FILE, "<", $filename);
while (my $line = <FILE>)
{
    chomp($line);
    #Found reciprocol top hits:      orf00003|scf7180000000008|gid_14        PA1392|orth_id:1357|gid_1
    my ($junk, $orf_info, $id) = split('\t', $line);
    my ($orf_locus_tag, $contig_id, $orf_gid) = split('\|', $orf_info);
    $orf_locus_tag =~ s/orf/$contig_id\_/;
    my ($bh_locus_tag, $orth_orf_id, $bh_gid) = split('\|', $id);
    $orth_orf_id =~ s/orth_id://;
    if ($orth_orf_id =~ /NO_GROUP/)
      {
	next;
      }
    $orf_gid =~ s/gid_//;
    print $orth_orf_id, "\t", $orf_locus_tag, "\t", $orf_gid, "\t", 5, "\n";
}


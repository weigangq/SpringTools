#!/usr/bin/perl

use strict;
use warnings;

my $file = shift @ARGV;

open(FILE, $file);
open(OUT, ">$file.parsed");
my $cluster_id;

while(<FILE>) {
    chomp;
    # KEEP TRACK OF CLUSTER ID
    if (/^>Cluster (\d+)/) {
	$cluster_id = $1;
    }
    elsif (/>(.+)\.\.\./) {
      my $unparsed = $1;
      my @fields = split('\|', $unparsed);
      # Adjust the cluster id by 1 since cdhit cluster numbers start at 0.
      print OUT join("\t", @fields), "\t", $cluster_id+1, "\n";
    }
}

close(FILE);
close(OUT);

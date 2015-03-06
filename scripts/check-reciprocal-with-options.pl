#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std; 

my %opts; 
getopts ('l:p:', \%opts) ;

die "$0 -l <coverage cuttoff> -p <percent id cutoff> <forward_blast_m6> <reverse_blast_m6>\n" unless @ARGV == 2;
my ($fwd, $rev) = @ARGV;
my $cutoff_perc_id = $opts{p} || 80;
my $cutoff_perc_length = $opts{l} || 0.80;
my (@fwd_top_hits, @rev_top_hits);
open FWD, "<" . $fwd;
my (%fwd_top_hits, %rev_top_hits, @query);
my $query_ct=0;
while (<FWD>) {
    chomp;
    my @data = split;
    next if $fwd_top_hits{$data[0]};
    $fwd_top_hits{$data[0]} = 
	{ 'hit' => $data[1],
	  'perc_id' => $data[2],
	  'perc_length' => sprintf "%.4f", $data[3]/$data[4],
	  'evalue' => $data[5],
	}; 
    push @query, $data[0];
    $query_ct++;
}
close FWD;
warn "Total query having hits:" . $query_ct . "\n";

open REV, "<" . $rev;
while (<REV>) {
    chomp;
    my @data = split;
    next if $rev_top_hits{$data[0]};
    $rev_top_hits{$data[0]} = $data[1]; 
}
close REV;

foreach my $q (@query) { # e.g., BafPKo_0002
    my $top = $fwd_top_hits{$q}; # e.g. BB_0002
    if ( $q eq $rev_top_hits{$top->{hit}}) {
	if ($top->{perc_id} >= $cutoff_perc_id && $top->{perc_length} >= $cutoff_perc_length ) {
	    print "Found reciprocol top hits:\t", $q, "\t", $top->{hit}, "\n";
	} else {
	    warn "Skipped for ", $q, ": reciprocal hits low identify or length coverage: ", $top->{hit}, "\t", $top->{perc_id}, "\t", $top->{perc_length},  "\n";
	}
    } else {
	warn "Not reciprocol top hits:\t", $q, "\t", $top->{hit}, "\t", $rev_top_hits{$top->{hit}}, "\n";
    }
}

exit;

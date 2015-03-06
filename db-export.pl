#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;

my %opts;
GetOptions(
    \%opts,
    "help|h",
    "man|m",
    "exclude_genome|G=s",
    "exclude_annot|A=s",
    "genome|g=s",
    "annot|a=s"
    ) or pod2usage(2);
pod2usage(1) if $opts{"help"};
pod2usage( -exitstatus => 0, -verbose => 2) if $opts{"man"};

my $dbh = DBI->connect("dbi:Pg:dbname=paerug;host=borreliabase.org", "lab", "homology", {RaiseError => 1, AutoCommit => 1}) || die "Could not connect to database.";

# EXCLUDE GENOME_ID
#SELECT cdhit_id, genome_id, count(genome_id) FROM orf WHERE exclude != TRUE GROUP BY cdhit_id, genome_id
my $statement = "SELECT cdhit_id, genome_id, count(genome_id) FROM orf WHERE exclude != TRUE GROUP BY cdhit_id, genome_id";

if (defined($opts{"exclude_genome"})) {
    my $excluded = $opts{"exclude_genome"};
    my @exclude_gid = split(",", $excluded);
    my @exclude_query;
    foreach my $gid (@exclude_gid) {
	if ($gid =~ /(\d+)\-(\d+)/) {
	    for (my $i = $1; $i < 1 + $2; $i++) {
		push @exclude_query, $i;
	    }
	}
	elsif ($gid =~ /\d+/) {
	    push @exclude_query, $gid;
	}
	else {
	    die "Error with excluded genome(s) syntax.\n";
	}
    }
    my $exclude = join(",", @exclude_query);
    $statement =~ s/WHERE/WHERE genome_id NOT IN \($exclude\) AND/;
#    print $statement, "\n";
}

#INCLUDE GENOME_ID
if (defined($opts{"genome"})) {
    my $included = $opts{"genome"};
    my @include_gid = split(",", $included);
    my @include_query;
    foreach my $gid (@include_gid) {
	if ($gid =~ /(\d+)\-(\d+)/) {
	    for (my $i = $1; $i < 1 + $2; $i++) {
		push @include_query, $i;
	    }
	}
	elsif ($gid =~ /\d+/) {
	    push @include_query, $gid;
	}
	else {
	    die "Error with included genome(s) syntax.\n";
	}
    }
    my $include = join(",", @include_query);
    $statement =~ s/WHERE/WHERE genome_id IN \($include\) AND/;
}

#EXCLUDE ANNOTATION
if (defined($opts{"exclude_annot"})) {
    my $excluded = $opts{"exclude_annot"};
    my @exclude_annot = split(",", $excluded);
    my @exclude_query;
    foreach my $annot (@exclude_annot) {
	$annot =~ s/^\s*(\S+)\s*$/$1/;
	$annot = "lower(product_name) NOT LIKE '%" . $annot . "%'";
	push(@exclude_query, $annot);
    }
    my $exclude = join(" AND ", @exclude_query);
    $statement =~ s/WHERE/WHERE \($exclude\) AND/;
#    print $statement, "\n";
}

#INCLUDE ANNOTATION
if (defined($opts{"annot"})) {
    my $included = $opts{"annot"};
    my @include_annot = split(",", $included);
    my @include_query;
    foreach my $annot (@include_annot) {
	$annot =~ s/^\s*(\S+)\s*/$1/;
	$annot = "lower(product_name) LIKE '%" . $annot . "%'";
	push(@include_query, $annot);
    }
    my $include = join(" OR ", @include_query);
    $statement =~ s/WHERE/WHERE \($include\) AND/;
#    print $statement, "\n";
}

print $statement, "\n";


#QUERY
my $query = $dbh->prepare($statement);
$query->execute();

my (@results, %seen_genome, %seen_cdhit_id);

while (my @data = $query->fetchrow_array()) {
    my ($cdhit_id, $genome_id, $count) = ($data[0], $data[1], $data[2]);
    $seen_genome{$genome_id}++;
    $seen_cdhit_id{$cdhit_id}++;
    push @results, {
	'cdhit_id' => $cdhit_id,
	'genome_id' => $genome_id,
	'cdhit_cts' => $count
    };
}

#FILTER
my $num_genomes = scalar keys %seen_genome;
my @filtered = &filter(\@results);
my %filtered_seen_cdhit;
foreach my $cdhit_id (@filtered) {
    $filtered_seen_cdhit{$cdhit_id->{cdhit_id}}++;
}

#PRINT MATRIX
print "gid\t", join("\t", sort {$a <=> $b} keys %filtered_seen_cdhit), "\n";

$query = $dbh->prepare('SELECT strain_name FROM genome WHERE genome_id = ?');

foreach my $gid (sort {$a <=> $b} keys %seen_genome) {
    $query->execute($gid);
    my @temp = $query->fetchrow_array();
    print $temp[0];
    foreach my $cdhit (sort {$a <=>$b} keys %filtered_seen_cdhit) {
	my @cts = grep {$_->{genome_id} eq $gid && $_->{cdhit_id} eq $cdhit} @filtered;
	if (@cts) {
	    my $ct = shift @cts;
	    print "\t", $ct->{cdhit_cts};
	}
	else {
	    print "\t0";
	}
    }
    print "\n";
}

sub filter() {
    my @Filtered;
    my @Nonfiltered = @{$_[0]};
    foreach my $cdhit_id (sort {$a <=> $b} keys %seen_cdhit_id) {
	my @cts = grep {$_->{cdhit_id} eq $cdhit_id} @Nonfiltered;
	if (scalar @cts < $num_genomes-1 && @cts > 1) {
	    push @Filtered, @cts;
	}
    }
    return @Filtered;
}

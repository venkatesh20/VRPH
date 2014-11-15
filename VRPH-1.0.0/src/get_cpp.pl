#!/usr/bin/perl

use Data::Dumper;

open(FILE,$ARGV[0]) or die $!;
my @cpp=();

while (my $line = <FILE>) {
  if ($line=~ m/^(\S+)\.cpp/i) {
     if (! ($1 ~~ @cpp) ) {
       push(@cpp,$1);
     }
  }
}


close(FILE);

print Dumper @cpp;
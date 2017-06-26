#!/usr/bin/perl -w
use strict;
use warnings;

my $pos = 0;
open(my $fd, $ARGV[0]) or die("can't open file\n");
while(my $l = <$fd>)
{
    chomp $l;
	my @sp = split('\t', $l);
	$pos += $sp[6];
}

close $fd;

my $self = 0;
open(my $fe, $ARGV[1]) or die("can't open file\n");
while(my $m = <$fe>)
{
    chomp $m;
	my @sq = split('\t', $m);
	$self += $sq[6];
}

close $fe;

my $TP = 0;
open(my $ff, $ARGV[2]) or die("can't open file\n");
while(my $n = <$ff>)
{
    chomp $n;
	my @sr = split('\t', $n);
	$TP += $sr[6];
}

close $ff;

if($TP == 0 || $self == 0 || $pos == 0)
{
    print("\nprecision: 0\nrecall: 0\n");
    
    print("F-measure with beta weight ",$ARGV[3],": 0\n\n");
}
else
{
  my $precision = $TP/($self);
  
  my $recall = $TP/$pos;
  
  my $F_measure = ((($ARGV[3])*($ARGV[3])+1)*$precision*$recall)/((($ARGV[3])*($ARGV[3])*$precision)+$recall);
  
  print("\nprecision: ",$precision,"\nrecall: ",$recall,"\n");
  
  print("F-measure with beta weight ",$ARGV[3],": ",$F_measure,"\n\n");
}

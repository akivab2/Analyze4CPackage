#!/usr/bin/perl -w
use strict;
use warnings;

open(my $fd, $ARGV[0]) or die("can't open file\n");
while(my $l = <$fd>)
{
    chomp $l;
	my @sp = split('\t', $l);
	my $sec = $sp[1] + 1;
	print"$sp[0]\t$sp[1]\t$sec\n";
}

close $fd;

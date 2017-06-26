#!/usr/bin/perl -w
use strict;
use warnings;

open(my $fd, $ARGV[0]) or die("can't open file\n");
while(my $l = <$fd>)
{
    chomp $l;
	my @sp = split('\t', $l);
  $sp[2] =~ s/\D//g; #the code removes '\n' from the end of $sp[2]. this is so the output won't have an extra \n.
	print"$sp[0]\t$sp[1]\t$sp[2]\t\n";
}

close $fd;

#!/usr/bin/perl -w
use strict;
use warnings;

my $counter = 1;
open(my $fd, $ARGV[0]) or die("can't open file\n");
while(my $l = <$fd>)
{
	if($counter % 2 == 1)
	{
		chomp $l;
		my @data1 = split('>',$l);
		my @data2 = split(':',$data1[1]);
		my @data3 = split('-',$data2[1]);
		print("$data2[0]\t$data3[0]\t$data3[1]\n");
	}
 $counter = $counter+1;
}
close $fd;

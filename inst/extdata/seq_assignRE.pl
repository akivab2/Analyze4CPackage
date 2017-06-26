#!/usr/bin/perl -w
use strict;
use warnings;

##############
#This script get- ARGV[1]:RE sites on a genome (strand    chr  spot), ARGV[0]:BED file with reads
#calculate the number of reads that are assign to each RE site.
#Output to STDOUT- chr  spot    numberOfReads
##############


#open BED file. fill Hash that the key is chr and the value is hash that the key is the start on this chr if the strand is + or the end if the strand is - (end+6).
open(my $fd, $ARGV[0]) or die("can't open file $ARGV[0]");
my %reads;
while(<$fd>)
{
    chomp;
    my $line=$_;
    my @data = split(' ',$line);
    if ($data[5] eq "+")
    {
        if (exists $reads{$data[0]}{$data[1]})
        {
            $reads{$data[0]}{$data[1]}+=1;
        }
        else
        {
            $reads{$data[0]}{$data[1]}=1;
        }
    }
    else
    {
        if (exists $reads{$data[0]}{$data[2]-4})
        {
            $reads{$data[0]}{$data[2]-4}+=1;
        }
        else
        {
            $reads{$data[0]}{$data[2]-4}=1;
        } 
    }
}
###print the proccesd BED file
#foreach my $chr (sort keys %reads)
#{
#    foreach my $pos (sort keys %{$reads{$chr}})
#    {
#        print "$chr\t$pos\t$reads{$chr}{$pos}\n"
#    }
#}
close ($fd);
#open file to write to
#open (my $wr, ">RE_".$ARGV[0]) or die ("can't open file for writing");
#open RE file
open($fd, $ARGV[1]) or die("can't open file $ARGV[1]");
while(<$fd>)
{
    chomp;
    my $line=$_;
    my @data = split(' ',$line);
    print "$data[1]\t$data[2]\t";
    if (exists $reads{$data[1]}{$data[2]})	
    {
          print "$reads{$data[1]}{$data[2]}\n";
    }
    else {print "0\n"};
}
close ($fd);

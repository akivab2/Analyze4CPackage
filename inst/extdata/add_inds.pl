#this code gets a bed file of RE sites and a number of bps to add and subtract
#the code adds that number ($ARGV[1]) to the second index and subtracts from the first in order
#to create the dpn and its surroundings on the right and to the left (since we could read it also from the 
#opposite strand and align backwards)

open(my $fd, $ARGV[0]) or die("can't open file\n");
while(my $l = <$fd>)
{
  chomp $l;
	my @data = split('\t',$l);
	my $first = $data[1]-$ARGV[1];
	my $second = $data[2]+$ARGV[1];
	print("$data[0]\t$data[1]\t$second\n$data[0]\t$first\t$data[2]\n");
}
close $fd;

#' @export

#this function creates a file in the folder "genomes", with all the chromosomes of an organism and gets from the user the sizes of each chromosome
#if the file already exists in "genomes" it will not create it

#the point of the flag is in case we want to create the file after we already asked the user for the name of the organism
#since we don't want the user to have to enter the name multiple times, the function already knows the name and uses it from the input
#when the flag equals 0 then the user hasn't been asked yet, when its equal to 1 the user has been asked
#when the input is 1 the other inputs should also be the name of the organism and the whole name of the chromosome size file
#the default is just "null" for each, will be changed when the user writes the actual names
chrom_sizes_fileCreator <- function(flag=0,name1="null1",name2="null2")
{
	setTOroot()
	if(flag==0)
	{
		name1 <- readline(prompt="what is the name of the genome that you want to create a size file for?\n")
		name2 <- paste(name1,"_chromosome_sizes.txt",sep="")
		ls_files <- system("ls ~/Analyze4C/genomes",intern=TRUE)
		match <- pmatch(name2,ls_files,nomatch=0)
	}	
	else
	{
		match <- 0
	}
	
	if(match == 0)
	{
		len <- readline(prompt="please enter the number of chromosomes in the genome:\n")
		len <- as.integer(len)
		all <- matrix(0,len,2)
		for(i in 1:len)
		{
			
			siz <- readline(prompt=paste("what is the size of chromosome",i,"?\n"))
			ch <- paste("chr",i,sep="")
			all[i,1] <- ch
			all[i,2] <- siz
		}
		setwd("genomes")
		write.table(all, name2 , sep="\t",  row.names = FALSE,col.names = FALSE,quote=FALSE)
		setTOroot()
	}
	else
	{
		cat("\n",name1,"'s chromosome sizes file already exists\n\n",sep="")
	}
}

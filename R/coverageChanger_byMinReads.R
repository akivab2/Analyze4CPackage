#' @export

#this function changes the coverage by changing the number of reads to 0 on any RE site that has a number of reads less than what the user requested

coverageChanger_byMinReads <- function(dat,min_reads,cis,choice=0,ch=0)
{
	ln <- length(unique(dat[,1])) #getting the number of chromosomes
	new <- dat
	if(choice == 0)
	{	
		choice <- as.integer(readline(prompt="What chromosomes would you like to change the minimum reads for:\n\n1) All\n2) trans\n3) cis\n4) by chromosome\n\nPlease enter your choice\n\n"))
	}	
	switch(choice,
	#all
	{
		new[new[,3]<min_reads,3] <- 0
		type <- "all chromosomes minimum reads changed"
	},
	#trans
	{
		new[new[,1] != cis & new[,3]<min_reads,3] <- 0
		type <- "trans chromosomes minimum reads changed"
	},
	#cis
	{
		new[new[,1] == cis & new[,3]<min_reads,3] <- 0
		type <- "cis chromosome minimum reads changed"
	},
	#by chromosome
	{
		if(ch==0)
		{
			ch <- as.integer(readline(prompt="Enter the chromosome number that you would like to change change the minimum reads for:\n"));
		}
		new[new[,1] == ch & new[,3]<min_reads,3] <- 0
		type <- paste("chromosome",ch,"minimum reads changed")
	}
	)
	out <- list(new,type)
	return(out)
}	

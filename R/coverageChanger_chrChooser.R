#' @export

coverageChanger_chrChooser <- function(filename,remove,cis,choice=0,ch=0)
{
	ln <- length(unique(filename[,1])); #getting the number of chromosomes
	if(choice == 0)
	{	
		choice <- as.integer(readline(prompt="Which coverage would you like to change:\n\n1) All\n2) trans\n3) cis\n4) by chromosome\n\nPlease enter your choice\n\n"));
	}	
	switch(choice,
	#all
	{
		data <- filename;
		new <- cov.remover(data,remove);
		type <- "all chromosomes coverage changed"
	},
	#trans
	{
		data <- filename[filename[,1]!=cis,]; #getting all the trans data
		temp <- cov.remover(data,remove);
		new <- c();
		for(i in 1:ln)
		{
			if(i!=cis)
			{
				new <- rbind(new,temp[temp[,1]==i,]);
			}
			else
			{
				new <- rbind(new,filename[filename[,1]==i,]);
			}
		}
		type <- "trans chromosomes coverage changed"
	},
	#cis
	{
		data <- filename[filename[,1]==cis,];
		temp <- cov.remover(data,remove);
		new <- c();
		for(i in 1:ln)
		{
			if(i!=cis)
			{
				new <- rbind(new,filename[filename[,1]==i,]);
			}
			else
			{
				new <- rbind(new,temp[temp[,1]==i,]);
			}
		}
		type <- "cis chromosome coverage changed"
	},
	#by chromosome
	{
		if(ch==0)
		{
			ch <- as.integer(readline(prompt="Enter the chromosome number that you would like to change its coverage:\n"));
		}	
		data <- filename[filename[,1]==ch,];
		temp <- cov.remover(data,remove);
		new <- c();
		for(i in 1:ln)
		{
			if(i!=ch)
			{
				new <- rbind(new,filename[filename[,1]==i,]);
			}
			else
			{
				new <- rbind(new,temp[temp[,1]==i,]);
			}
		}
		type <- paste("chromosome",ch,"coverage removed")
	}
	)
	out <- list(new,type)
	return(out);
}	

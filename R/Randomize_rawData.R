#' @export

#this function recieves a raw data file and randomizes the reads on the RE sites per chromosome
#the purpose is to have random data that can be compared to the original data, usually proving that the original datas results are not random

Randomize_rawData <- function()
{
	#get the data file
	repeat
	{
		choice <- as.integer(readline(prompt=cat("\nenter the number of the folder from which you would like to take the raw data file from:\n\n1) original\n2) rearranged\n3) coverage removed\n\n")))		
		if(choice == 2)
		{
			ls_files <- system("ls ~/Analyze4C/rawData/rearranged",intern=TRUE)
			pth <- "~/Analyze4C/rawData/rearranged/"
		}
		else if(choice == 3)
		{
			ls_files <- system("ls ~/Analyze4C/rawData/coverage_removed",intern=TRUE)
			pth <- "~/Analyze4C/rawData/coverage_removed/"
		}
		else
		{
			ls_files <- system("ls ~/Analyze4C/rawData/original",intern=TRUE)
			pth <- "~/Analyze4C/rawData/original/"
		}
		
		if(length(ls_files) != 0)
		{	
			file.name <- readline(prompt=cat("\nThese are the raw data files available\nwhich would you like to use?\n",ls_files,"\n",sep="\n"))
			ind_files <- pmatch(file.name,ls_files,nomatch=0)
			if(ind_files == 0)
			{
				cat("no such file exists.\nplease try again.\n\n")
			}
			else
			{
				break
			}
		}
		else
		{
			cat("\nthis folder is empty. please choose a different folder from which you would like to use a file\n\n")
		}
	}
	
	#importing the data from the file
	if(choice == 2)
	{
		dat <- read.table(paste("~/Analyze4C/rawData/rearranged/",file.name,sep=""))
	}
	else if(choice == 3)
	{
		dat <- read.table(paste("~/Analyze4C/rawData/coverage_removed/",file.name,sep=""))
	}
	else
	{
		dat <- read.table(paste("~/Analyze4C/rawData/original/",file.name,sep=""))
	}
	
	#randomize the data (per chromosome)
	cat("\nrandomizing the data\nplease wait...\n\n")
	randAll <- c();
	for (r in unique(dat[,1]))
	{
		rand <- sample(dat[dat[,1]==r,3]); 
		randAll	<- c(randAll,rand)
	}
	rand_final <- cbind(dat[,1],dat[,2],randAll);
			
	#saving the file
	DandT1 <- toString(Sys.time())
	DandT2 <- gsub(" ","_",DandT1)
	DandT2 <- gsub(":","",DandT2)	
	spl1 <- unlist(strsplit(file.name,"[.]"))[1]
	write.table(rand_final,paste(pth,spl1,"_randomized_",DandT2,".sgr",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)	
}

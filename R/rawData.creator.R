#' @export

#inputs: 'Experiments_4C', this is in order to present to the user the names of the experiments that they could choose from in order to create the raw data file
#this function first gets the RE site name and organism from user
#then it will check to see if RE site files are created already for these, and if not it will create them using 'RE.file.creator' function
#the user is then asked which of the RE files they would like to use, if the reads file of that experiment exists it will create a raw data sgr file
#using 'seq_assignRE.pl' and the inputs: the aligned reads bed file of the experiment ([ARGV[0]), and the RE site file txt file (ARGV[1])

rawData.creator <- function(Experiments_4C)
{
	#getting the name of RE site and organism
	RE_name <- readline(prompt="\nplease provide the RE sites name:\n\n")
	RE_org <- readline(prompt="\nplease provide the organisms name:\n\n")
	
	RE.file.creator(RE_name,RE_org)
	
	##################################################################################################################################################	
	
	#asking which of the files in sgr_orgs folder the user wants to use
	ls_REs <- system("ls ~/Analyze4C/sgr_orgs",intern=TRUE)
	RE_file <- readline(prompt=cat("which RE file would you like to use?\n",ls_REs,"\n",sep="\n"))
	
	#here we create the actual sgr file (raw data) for the experiment
	#for this we need the bed file of the alignments (not reads like what is here, this needs to be changed) and the txt file of the RE sites
	
	cat("the 4C experiments that are available are:\n\n")
	print(Experiments_4C)
	cat("\n")
	
	inp2_bait <- readline(prompt="enter the name of the bait:\n")
	inp2_tissue <- readline(prompt="enter the name of the tissue:\n")
	inp2_lane <- readline(prompt="enter the name of the lane:\n")
	
	reads_name <- paste(inp2_bait,"_",inp2_tissue,"_",inp2_lane,"_AlignedReads.bed",sep="")
	ls_reads <- system("ls ~/Analyze4C/Reads",intern=TRUE)
	ind_reads <- pmatch(reads_name,ls_reads,nomatch=0)
	if(ind_reads == 0)
	{
		cat("\nthis experiments read files don't exist\nfirst extract read files for this experiment, align them, and create a bed file, and then try again\n\n")
	}
	else
	{
		#creating the file name according to the RE file and information given above
		spl1 <- strsplit(RE_file,"[.]")
		if(length(spl1[[1]])==3)
		{
			name_sgr <- paste(inp2_bait,"_",inp2_tissue,"_",inp2_lane,".sgr",sep="")
		}
		else if(length(spl1[[1]])==4)
		{
			name_sgr <- paste(inp2_bait,"_",inp2_tissue,"_",inp2_lane,"_",spl1[[1]][3],".sgr",sep="")				
		}
		else if(length(spl1[[1]])==5)
		{
			name_sgr <- paste(inp2_bait,"_",inp2_tissue,"_",inp2_lane,"_",spl1[[1]][3],"_",spl1[[1]][4],".sgr",sep="")
		}
		#creating the raw data file
		system(paste("perl ~/Analyze4C/proxy/seq_assignRE.pl ~/Analyze4C/Reads/",reads_name," ~/Analyze4C/sgr_orgs/",RE_file," > ~/Analyze4C/rawData/original/",name_sgr,sep=""))
	}
}

#' @export

#1)extract reads from sequencing

#explanation:
#this function gets the name of the bait and tissue from user (details of a specific experiment),
#then creates a text file and setting its name accordingly (if a file with that name already exists then it will not create a new file),
#then gets the primer for that experiment, and finally extracts the reads from the sequencing fastq files (as requested by user)

extract_reads <- function(Experiments_4C)
{
	#making sure the sequencing files are in the right folder
	repeat
	{
		inp1_start <- readline(prompt="Please make sure the sequencing files are in the folder called 'Seq'\n\npress 1 to continue\n")
		inp1_start <- as.integer(inp1_start)
		if(inp1_start == 1)
		{break}
	}
	#might need to add here a condition that if the files are moved the user should press 1 or something like that and only then the code will work

	#####################################################################################################################
	
	#creating the output file (inp1 means that it is the input from user in option 1), this is how i should do it for all the other options as well)
	#if the output file already exists then there is no need to create it
	inp1_bait <- readline(prompt="enter the name of the bait:\n")
	inp1_tissue <- readline(prompt="enter the name of the tissue:\n")
	inp1_lane <- readline(prompt="enter the name of the lane:\n")
	
	#creating name of file and output file (if doesn't exist yet)
	name1 <- paste(inp1_bait,"_",inp1_tissue,"_",inp1_lane,"_","reads.txt",sep="")
	ls_reads <- system("ls ~/Analyze4C/Reads",intern=TRUE)
	ind_reads <- pmatch(name1,ls_reads,nomatch=0)
	if(ind_reads == 0)
	{
		system(paste("touch ~/Analyze4C/Reads/",name1,sep=""))
	}
	
	#####################################################################################################################
		
	#here we check if this experiment is in 'Experiments_4C', if not we need to add it to the list
	#the function returns the updated Experiments_4C and the primer that was added (it's returned in a list so make sure to separate them later)
	
	out1 <- findIn.Experiments_4C(inp1_bait,inp1_tissue,inp1_lane,Experiments_4C)
	Experiments_4C <- out1[[1]]
	inp1_primer <- out1[[4]]
	
	#if i want to also get an excel file with the list i could try this
	#first print out the list of primers and then tell user to choose one from the list and copy it
	#try using this https://www.datacamp.com/community/tutorials/r-tutorial-read-excel-into-r

	#i think that here i need to update the Experiments_4C file
	
	#####################################################################################################################
		
	#extracting the reads
	repeat
	{
		cat("\n")
		system("ls ~/Analyze4C/Seq | grep '.fastq$'")
		inp1_seqFile <- readline(prompt="\nchoose a sequencing file from the list to extract from.\ncopy the name of the file and enter it into the prompt line\n\n")
		system(paste("grep -A 2 -B 1 ^",inp1_primer," ~/Analyze4C/Seq/",inp1_seqFile," >> ~/Analyze4C/Reads/",name1,sep=""))
		inp1_continue <- readline(prompt="would you like to extract more reads? y/n\n")
		if(inp1_continue == "n")
		{break}
	}
	
	#removing '--' from file
	system(paste("sed --in-place '/--/d' ~/Analyze4C/Reads/",name1,sep=""))
	
	#returning Experiments_4C
	return(Experiments_4C)
}

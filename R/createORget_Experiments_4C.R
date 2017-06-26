#' @export

#creating the file "Experiments_4C" if it doesn't exist
#the file contains the list of experiments which contains the bait, tissue, experiment name or lane
#the list will be filled in by user while they use the program, every time they enter a new experiment it will add a new row

createORget_Experiments_4C <- function()
{
	ls_files2 <- system("ls ~/Analyze4C",intern=TRUE)
	ind2 <- pmatch("Experiments_4C.txt",ls_files2,nomatch=0)
	if(ind2 != 0) #if the file exists we will get the file and it will be sorted
	{
		Experiments_4C <- read.table("Experiments_4C.txt",header=TRUE,stringsAsFactors=FALSE)
		#sorting the list of experiments by bait alphabetically (and sorting the row indices)
		Experiments_4C <- Experiments_4C[order(Experiments_4C$Bait),]
		rownames(Experiments_4C) <- seq(length=nrow(Experiments_4C))
	}
	else #if the file does not exist it will be created
	{
		Experiments_4C <- data.frame(Bait=character(),Tissue=character(),Lane_Experiment=character(),Cis=character(),Bait_position=character(),Primer=character(),Primer_Length=character(),stringsAsFactors=FALSE)
		write.table(Experiments_4C, "Experiments_4C.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
	}
	return(Experiments_4C)
}

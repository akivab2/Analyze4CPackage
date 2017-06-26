#' @export
#aligning reads

aligner <- function(Experiments_4C)
{
	inp_bait <- readline(prompt="enter the name of the bait:\n")
	inp_tissue <- readline(prompt="enter the name of the tissue:\n")
	inp_lane <- readline(prompt="enter the name of the lane:\n")
	
	name <- paste(inp_bait,"_",inp_tissue,"_",inp_lane,"_reads",sep="")
	
	#check to see that there is no bed file existing already
	
	out <- findIn.Experiments_4C(inp_bait,inp_tissue,inp_lane,Experiments_4C)
	primer <- out[[4]]
	
	cat("the primer is",primer,"\n")
	
	#prim_length <- readline(prompt="\nplease enter the length of the primer not including the RE site\n")
	prim_length <- out[[5]]
	
	link <- readline(prompt="\nplease provide a link to the whole genomes bowties index file and add '/genome' at the end, this should be in the bowtie index folder of the genomes sequence folder\ne.g.: Arabidopsis/Arabidopsis_thaliana_TAIR10/Sequence/BowtieIndex/genome\n\n")

	#aligning the reads using bowtie
	cat("\naligning the reads\nplease wait...\n\n")
	system(paste("bowtie -m 1 -q -S --trim5 ",prim_length," ",link," ~/Analyze4C/Reads/",name,".txt"," ~/Analyze4C/Reads/",name,".sam",sep=""))
	
	#converting the sam file to a bed file
	name2 <- paste(inp_bait,"_",inp_tissue,"_",inp_lane,"_AlignedReads",sep="")
	cat("\nconverting the sam file to a bed file\nplease wait...\n\n")
	system(paste("sam2bed"," ~/Analyze4C/Reads/",name,".sam"," ~/Analyze4C/Reads/",name2,".bed",sep=""))
	
	#remove the sam file
	system(paste("rm"," ~/Analyze4C/Reads/",name,".sam",sep=""))
}	


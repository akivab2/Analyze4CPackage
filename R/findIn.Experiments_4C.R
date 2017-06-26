#' @export

#here we check if this experiment is in 'Experiments_4C', if not we need to add it to the list
#the function returns the updated Experiments_4C and the primer that was added (it's returned in a list so make sure to separate them later)

findIn.Experiments_4C <- function(bait,tissue,lane,Experiments_4C)
{
	exp_ind <- which(Experiments_4C$Bait %in% bait & Experiments_4C$Tissue %in% tissue & Experiments_4C$Lane_Experiment %in% lane)
	if(length(exp_ind) == 0)
	{
		#asking if to add this experiment to the Experiments_4C file
		inp <- readline(prompt="\nthe details of this experiment do not exist in the file.\n\nwould you like to add them to the file?\ny/n\n\n")
		if(inp == "y")
		{
			#getting the cis chromosome number from user
			cis <- readline(prompt="enter the cis chromosome of the experiment:\n")
			#getting the bait position from user
			bait_pos <- readline(prompt="enter the bait position:\n")
			#getting the primer from user
			primer <- readline(prompt="enter the primer of the experiment:\n")
			#getting the size of the primer (not including the RE site)
			primer_len <- readline(prompt="enter the length of the primer of the experiment (not including the RE site):\n")
			#adding the new experiment to the list
			Experiments_4C[nrow(Experiments_4C)+1,] <- c(bait,tissue,lane,cis,bait_pos,primer,primer_len)
			#sorting the list of experiments by bait alphabetically (and sorting the row indices)
			Experiments_4C <- Experiments_4C[order(Experiments_4C$Bait),]
			rownames(Experiments_4C) <- seq(length=nrow(Experiments_4C))
			#adding the new data to the file (by erasing the old one and creating a new one)
			system("rm Experiments_4C.txt")
			write.table(Experiments_4C,"Experiments_4C.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
			
			out <- list(Experiments_4C,cis,bait_pos,primer,primer_len)
			return(out)
		}
		else
		{
			cis <- NA
			bait_pos <- NA
			primer <- NA
			primer_len <- NA
			out <- list(Experiments_4C,cis,bait_pos,primer,primer_len)
			return(out)
		}
	}
	else
	{
		#getting the primer from the experiment list if exists already in list
		cis <- Experiments_4C$Cis[exp_ind]
		bait_pos <- Experiments_4C$Bait_position[exp_ind]
		primer <- Experiments_4C$Primer[exp_ind]
		primer_len <- Experiments_4C$Primer_Length[exp_ind]
		
		out <- list(Experiments_4C,cis,bait_pos,primer,primer_len)
		return(out)	
	}
}

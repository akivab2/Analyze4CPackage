#' @export

#this function tests if the bed file data might be translocated
#if it is, it asks if the user would want to remove the the added areas - they can do it manually or automatically (just removing the first section that was added)
#the function does this by chromosome
#the function returns the data after being treated with removing (or not if user choses not to) the translocated areas
#the input to the function has to be in bed file format

#note:
#this function will only be effective if the added translocated areas have indices that put the whole choromosomes indices out of order,
#but if they are all less than the smallest index of the real chromosomes indices, then the function won't be able to identify it

remove_translocated <- function(bed_dat,remove="y",ans1=0) #in the function "tissue_Expression_comparison" the value "y" should be entered, meaning that only the automatic removal should be done
{
	final <- c()
	chroms <- unique(bed_dat[,1])
	num_of_chr <- length(chroms)
	for(i in 1:num_of_chr)
	{
		da <- bed_dat[bed_dat[,1]==i,]
		da_first <- da[,2]
		len <-nrow(da)
		if(len==1)
		{
			final <- rbind(final,da)
		}
		else if(len!=0)
		{	
			difs <- diff(da_first)
			find_dif <- which(difs<=0)
			if(identical(find_dif,integer(0))=="FALSE")
			{
				cat("\nchromosome",i,"is unsorted, or with the same indices in a row\nthis could be an indication that the chromosome has undergone a translocation\n\n")
				if(remove != "y")
				{
					remove <- readline(prompt=cat("\nwould you like to remove the translocated areas\ny/n\n\n"))
				}
				
				if(remove == "y")
				{
					if(ans1==0)
					{
						ans1 <- as.integer(readline(prompt=cat("\nhow would you like to remove the section?\n1) manually\n2) automatically (the first section of the chromosome that seems to have been translocated will be removed)\n\n")))
					}
					
					if(ans1 == 1) #removing manually
					{
						cat("\nremoving the translocated area of chromosome",i,"automatically:\n\n")
						flag <- "y"
						write.table(da,"~/Analyze4C/temp/remove_indices_temp.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
						cat("\nin order to remove the wanted sections, view the file 'remove_indices_temp.bed' which is in the 'temp' folder\n")
						cat("enter the indices of the sections that you would like to KEEP\nthe sections will be kept in the order that you enter the indices\n")
						while(flag=="y")
						{
							first <- as.integer(readline(prompt=cat("\nenter the first index of the section would like to KEEP:\n")))
							last <- as.integer(readline(prompt=cat("\nenter the last index of the section would like to KEEP:\n")))
							if(first < last & first > 0 & last <= nrow(da))
							{
								final <- rbind(final,da[first:last,])
							}
							else
							{
								cat("\nthe indices entered were illegal\n")
							}
							flag <- readline(prompt=cat("\nwould you like to keep another section?\ny/n\n"))
						}
						system("rm ~/Analyze4C/temp/remove_indices_temp.bed")
					}
					else if(ans1 == 2)#removing automatically
					{
						cat("\nremoving the translocated area of chromosome",i,"automatically...\n\n")
						final <- rbind(final,da[(find_dif[1]+1):len,])
					}
				}
			}
			else
			{
				final <- rbind(final,da)
			}
		}
	}
	rownames(final) <- seq(length=nrow(final))
	return(final)
}

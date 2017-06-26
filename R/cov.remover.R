#' @export

cov.remover <- function(data,remove)
{
	data_pos <- data[which(data[,3]>0),]; #getting the positive REs
	probs <- 1-(data_pos[,3]/sum(data_pos[,3])); #getting the relative part of reads for each RE site from sum of reads. we do 1 minus this number since we want those with less reads to get a higher score.
	n <- round(nrow(data_pos)*remove); #getting the number of REs that need to be removed (meaning that they were positive and now will get a 0)
	samps <- sample(data_pos[,2],n,replace=FALSE,prob=probs) #sampling from the positive REs, those with a lower number of reads have a higher chance of being chosen here
	inds <- match(samps,data_pos[,2]) #getting the indices (indices from the data frame "data_pos") of the sampled REs
	data_pos[inds,3] <- 0 #giving the sampled REs a 0
	data[which(data[,3]>0),3] <- data_pos[,3] #replacing the previous number of reads with the updated one
	return(data);
}

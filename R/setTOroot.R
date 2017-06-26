#' @export

#explanation:
#this function checks the current directory to see if it is the root "Analyze4C"
#if not it will change the directory to it

setTOroot <- function()
{
	cur_dir <- system("pwd | rev | cut -d '/' -f1 | rev",intern=TRUE)
	while(cur_dir != "Analyze4C")
	{
		setwd("..")
		cur_dir <- system("pwd | rev | cut -d '/' -f1 | rev",intern=TRUE)
	}
}

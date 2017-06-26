#' @export

#this function is supposed to be done only the first time the program is activated, when the 'Analyze4C' folder is created
#the function copies all the perl files into the proxy folder for later use
create_perl_files <- function()
{
	add_inds.path <- system.file("extdata", "add_inds.pl", package = "Analyze4CPackage")
	system(paste("cp",add_inds.path,"~/Analyze4C/proxy/"))

	edit_notAligned.path <- system.file("extdata", "edit_notAligned.pl", package = "Analyze4CPackage")
	system(paste("cp",edit_notAligned.path,"~/Analyze4C/proxy/"))

	EOL_tab_adder.path <- system.file("extdata", "EOL_tab_adder.pl", package = "Analyze4CPackage")
	system(paste("cp",EOL_tab_adder.path,"~/Analyze4C/proxy/"))
	
	precisionRecallFmeasure.path <- system.file("extdata", "precisionRecallFmeasure.pl", package = "Analyze4CPackage")
	system(paste("cp",precisionRecallFmeasure.path,"~/Analyze4C/proxy/"))

	seq_assignRE.path <- system.file("extdata", "seq_assignRE.pl", package = "Analyze4CPackage")
	system(paste("cp",seq_assignRE.path,"~/Analyze4C/proxy/"))

	sgr_to_bed.path <- system.file("extdata", "sgr_to_bed.pl", package = "Analyze4CPackage")
	system(paste("cp",sgr_to_bed.path,"~/Analyze4C/proxy/"))	
}




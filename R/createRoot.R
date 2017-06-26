#' @export

#if "Analyze4C" is not created then create the folder from the root and create all the internal folder as well
createRoot <- function()
{
	ls_files1 <- system("ls",intern=TRUE)
	ind1 <- pmatch("Analyze4C",ls_files1,nomatch=0)
	if(ind1 == 0)
	{
		system("mkdir Analyze4C")
		#enter the folder "Analyze4C"
		setwd("Analyze4C")
		system("mkdir Reads")
		system("mkdir Seq")
		system("mkdir genomes")
		system("mkdir rawData")
		system("mkdir rawData/original")
		system("mkdir rawData/coverage_removed")
		system("mkdir rawData/rearranged")
		system("mkdir sgr_orgs")
		system("mkdir sgr_exps")
		system("mkdir pScores")
		system("mkdir pScores/with_bp")
		system("mkdir pScores/no_bp")
		system("mkdir contact_bands")
		system("mkdir plots")
		system("mkdir temp")
		#system("mkdir lib") #i might not need this when it is a package
		system("mkdir proxy")
		system("mkdir RNAseq")
		system("mkdir RNAseq/sra")
		system("mkdir RNAseq/FPKM")
		system("mkdir RNAseq/temp")
		system("mkdir FDR")
		system("mkdir ChIPseq")
		system("mkdir ChIPseq/SRA_or_FASTQ")
		system("mkdir ChIPseq/peaks")
		system("mkdir ChIPseq/temp")

		#moving the perl files into proxy
		create_perl_files()	
		
		#creating coverageChanged_Experiments.txt for coverage changing records
		coverageChanged_Experiments <- data.frame(Experiment=character(),Date_and_Time=character(),Removal_type=character(),Removal_location=character(),Percentage_Removed=character(),Minimum_reads=character(),stringsAsFactors=FALSE)
		write.table(coverageChanged_Experiments, "coverageChanged_Experiments.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating rearranged_rawData.txt for rearranged data records
		rearranged_rawData <- data.frame(Experiment=character(),Date_and_Time=character(),Description=character(),Added_Sections_Lines=character(),Note=character(),stringsAsFactors=FALSE)
		write.table(rearranged_rawData, "rearranged_rawData.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating pScores_of_rearranged.txt for rearranged p-score records
		pScores_of_rearranged <- data.frame(Experiment=character(),Date_and_Time=character(),Description=character(),Note=character(),stringsAsFactors=FALSE)
		write.table(pScores_of_rearranged, "pScores_of_rearranged.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)	
		
		#creating coverageVSfmeasure_plots.txt for coverage vs f-measure plots
		coverageVSfmeasure_plots <- data.frame(Coverage_Removal_Experiment=character(),Coverage_Removal_pScore_RE_Window=character(),Coverage_Removal_pScore_bp_Window=character(),Coverage_Removal_RE_Gap=character(),Coverage_Removal_Bp_gap=character(),Coverage_Removal_CO_method=character(),Coverage_Removal_CO=character(),Static_Experiment=character(),Static_pScore_RE_Window=character(),Static_pScore_bp_Window=character(),Static_RE_Gap=character(),Static_Bp_gap=character(),Static_CO_method=character(),Static_CO=character(),Limit_Method=character(),Limit=character(),Cycles=character(),Removal_location=character(),Coverage_Removal_Step=character(),Fmeasure_chromosomes=character(),Beta=character(),Date_and_Time=character(),stringsAsFactors=FALSE)
		write.table(coverageVSfmeasure_plots, "coverageVSfmeasure_plots.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)

		#creating combined_files.txt for combining raw data
		combined_files <- data.frame(Raw1=character(),Raw2=character(),Date=character(),stringsAsFactors=FALSE)
		write.table(combined_files, "combined_files.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating REwindowIntersections_plots.txt for RE windows vs precision, recall, f-measure, and intersection plots
		REwindowIntersections_plots <- data.frame(filename_1=character(),min_RE1=character(),max_RE1=character(),step_RE1=character(),REperBp_1=character(),contactBands_RE_gap_1=character(),contactBands_bp_gap_1=character(),CO_type1=character(),CO1=character(),filename_2=character(),min_RE2=character(),max_RE2=character(),step_RE2=character(),REperBp_2=character(),contactBands_RE_gap_2=character(),contactBands_bp_gap_2=character(),CO_type2=character(),CO2=character(),beta=character(),DandT1=character(),stringsAsFactors=FALSE)
		write.table(REwindowIntersections_plots, "REwindowIntersections_plots.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)

		#creating coIntersections_plots.txt for RE windows vs precision, recall, f-measure, and intersection plots
		coIntersections_plots <- data.frame(filename_1=character(),CO_type1=character(),min_CO1=character(),max_CO1=character(),step_CO1=character(),REperBp_1=character(),contactBands_RE_gap_1=character(),contactBands_bp_gap_1=character(),RE1=character(),filename_2=character(),CO_type2=character(),min_CO2=character(),max_CO2=character(),step_CO2=character(),REperBp_2=character(),contactBands_RE_gap_2=character(),contactBands_bp_gap_2=character(),RE2=character(),beta=character(),DandT1=character(),stringsAsFactors=FALSE)
		write.table(coIntersections_plots, "coIntersections_plots.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
				
		#creating RNAseq_data.txt for creating RNAseq data
		RNAseq_data <- data.frame(Name_of_File=character(),Organism=character(),Tissue=character(),Age=character(),Link=character(),Note=character(),Download_Date=character(),stringsAsFactors=FALSE)
		write.table(RNAseq_data, "RNAseq_data.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating tophat_cufflinks_outputs.txt for tophat and cufflinks commands, this way we can follow what commands and parameters were used for any specific analysis using tophat and cufflinks
		tophat_cufflinks_outputs <- data.frame(Name=character(),Date=character(),Tophat_Command=character(),Cufflinks_Command=character(),stringsAsFactors=FALSE)
		write.table(tophat_cufflinks_outputs, "tophat_cufflinks_outputs.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)

		#creating expressionVScontacts_plots.txt for expression (FPKM) vs. contact bands plots
		expressionVScontacts_plots <- data.frame(Plotfile_name=character(),Test_Structure=character(),Query_file=character(),Intersected_file1=character(),Intersected_file2=character(),BP_or_FPKM=character(),PercentageOf=character(),Date_and_Time=character(),FPKM_CO_type=character(),FPKM_CO=character(),FPKM_CO_appliedTo=character(),Intersection_chromosomes=character(),Notes=character(),stringsAsFactors=FALSE)
		write.table(expressionVScontacts_plots, "expressionVScontacts_plots.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)

		#creating expressionVScontacts_sumOFintersections_plots.txt for expression and contacts comparing by summing intersections
		expressionVScontacts_sumOFintersections_plots <- data.frame(Date_and_Time=character(),FPKM_file=character(),Num_of_experiments=character(),Experiments=character(),intersections_chroms=character(),summed_chroms=character(),RE_gap=character(),bp_gap=character(),stringsAsFactors=FALSE)
		write.table(expressionVScontacts_sumOFintersections_plots,"expressionVScontacts_sumOFintersections_plots.txt", sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)		

		#creating expressionVScontacts_correlation_plots.txt for comparing expression and contacts by correlating them
		expressionVScontacts_correlation_plots <- data.frame(pScore_or_reads=character(),num_of_data_files=character(),data_Experiments=character(),num_of_expression_files=character(),expression_Experiments=character(),data_CO_type=character(),data_CO=character(),data_CO_chroms=character(),expression_CO_type=character(),expression_CO=character(),expression_CO_chroms=character(),intersections_chroms=character(),sections_removed=character(),divided_by_windows=character(),window_sizes=character(),Date_and_Time=character(),stringsAsFactors=FALSE)
		write.table(expressionVScontacts_correlation_plots,"expressionVScontacts_correlation_plots.txt", sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating ChIPseq_data.txt for creating ChIPseq data
		ChIPseq_data <- data.frame(Name_of_File=character(),Organism=character(),Tissue=character(),Age=character(),Link=character(),Note=character(),Download_Date=character(),stringsAsFactors=FALSE)
		write.table(ChIPseq_data, "ChIPseq_data.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)

		#creating MACS_outputs.txt for MACS commands, this way we can follow what commands and parameters were used for any specific analysis using MACS
		MACS_outputs <- data.frame(Name=character(),Date=character(),MACS_Command=character(),stringsAsFactors=FALSE)
		write.table(MACS_outputs, "MACS_outputs.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating ChIPseqVScontacts_plots.txt for ChIPseq (peaks) vs. contact bands plots
		ChIPseqVScontacts_plots <- data.frame(Plotfile_name=character(),Test_Structure=character(),Query_file=character(),Intersected_file1=character(),Intersected_file2=character(),BP_or_peaks=character(),PercentageOf=character(),Date_and_Time=character(),ChIPseq_CO_source=character(),ChIPseq_CO_type=character(),ChIPseq_CO=character(),ChIPseq_CO_appliedTo=character(),ChIPseq_value_source=character(),Intersection_chromosomes=character(),Notes=character(),stringsAsFactors=FALSE)
		write.table(ChIPseqVScontacts_plots, "ChIPseqVScontacts_plots.txt", sep="\t",  row.names = FALSE,col.names = TRUE,quote=FALSE)
		
		#creating ChIPseqVScontacts_sumOFintersections_plots.txt for ChIPseq and contacts comparing by summing intersections
		ChIPseqVScontacts_sumOFintersections_plots <- data.frame(Date_and_Time=character(),peaks_file=character(),Num_of_experiments=character(),Experiments=character(),intersections_chroms=character(),summed_chroms=character(),RE_gap=character(),bp_gap=character(),quantile_source=character(),value_source=character(),stringsAsFactors=FALSE)
		write.table(ChIPseqVScontacts_sumOFintersections_plots,"ChIPseqVScontacts_sumOFintersections_plots.txt", sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)		

		#creating ChIPseqVScontacts_correlation_plots.txt for comparing ChIPseq and contacts by correlating them
		ChIPseqVScontacts_correlation_plots <- data.frame(pScore_or_reads=character(),num_of_data_files=character(),data_Experiments=character(),num_of_ChIPseq_files=character(),ChIPseq_Experiments=character(),data_CO_type=character(),data_CO=character(),data_CO_chroms=character(),ChIPseq_CO_source=character(),ChIPseq_CO_type=character(),ChIPseq_CO=character(),ChIPseq_CO_chroms=character(),intersections_chroms=character(),sections_removed=character(),divided_by_windows=character(),window_sizes=character(),ChIPseq_value_source=character(),Date_and_Time=character(),stringsAsFactors=FALSE)
		write.table(ChIPseqVScontacts_correlation_plots,"ChIPseqVScontacts_correlation_plots.txt", sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)		
	}
	else
	{
		#enter the folder "Analyze4C"
		setwd("Analyze4C")
	}
}

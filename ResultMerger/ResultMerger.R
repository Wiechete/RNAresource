# do merge of EdgeR results and APAlyzer results
## Check if libraries are installed in user's R environment

## Load relevant libraries
library(dplyr)

##### 1. Loading csv files #####
## Function to check and access script arguments
parse_args <- function(required_args = NULL) {
	# Get all user-provided arguments
	args <- commandArgs(TRUE)
	
	# Check if arguments are provided
	if (length(args) == 0) {
		cat("No arguments provided.\n")
		return(NULL)
	}
	
	# Check if the required number of arguments is provided
	if (!is.null(required_args) && length(args) < required_args) {
		stop(paste0("Insufficient arguments provided. Expecting at least ", required_args, 
								" arguments, but got ", length(args), "."))
	}
	
	# Return the arguments
	return(args)
}

## detect whether R script is launched in batch mode or RStudio mode
batch_mode_on <- is.na(Sys.getenv("RSTUDIO", unset = NA))
if(batch_mode_on){
	cat("Running in batch mode\r\n")
	arguments <- parse_args(2)
	
	DE_edger_result_path <- arguments[1]
	UTR_apalyzer_result_path <- arguments[2]
	
}else{
	cat("Running in iteractive mode\r\n")
	
	DE_edger_result_path <- paste(dirname(getwd()),"EdgeR/DE_hs_min1_ult.csv",sep = "/")
	UTR_apalyzer_result_path <- paste(dirname(getwd()),"APAlyzer/3UTR_APAdiff_results.csv",sep = "/")
}

## reading csv files
DE_edger_result <- read.csv(DE_edger_result_path)
UTR_apalyzer_result <- read.csv(UTR_apalyzer_result_path)

#### 2. Merge results by ensembl GeneID ####
## 
# Merge results with gene_mapping for IPA
final_result <- UTR_apalyzer_result %>%
	left_join(DE_edger_result, by = "Geneid") %>%  # Match using gene_symbol
	select(Geneid, gene_symbol, logFC, FDR, RED, pvalue, p_adj,DE, APAreg)

#### 3. Save merged dataframe to csv ####
write.csv(final_result, "final_result.csv", row.names = FALSE, quote = FALSE)

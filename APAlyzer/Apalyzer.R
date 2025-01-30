# 1. Loading relevant packages
# Function to check BiocManager version and update if needed
check_and_update_biocmanager <- function() {
	# Check if BiocManager is installed
	if (!requireNamespace("BiocManager", quietly = TRUE)) {
		message("BiocManager is not installed. Would you like to install it? (yes/no)")
		response <- tolower(readline())
		if (response == "yes") {
			install.packages("BiocManager")
			message("BiocManager has been installed.")
		} else {
			stop("BiocManager installation is required to proceed. Exiting.")
		}
	}
	
	# Check the installed version of BiocManager
	bioc_version <- packageVersion("BiocManager")
	required_version <- "1.30.0"
	
	# Compare versions and update if necessary
	if (bioc_version < required_version) {
		message("Your BiocManager version is ", bioc_version, ". A version of ", required_version, " or higher is required.")
		message("Would you like to update BiocManager? (yes/no)")
		response <- tolower(readline())
		if (response == "yes") {
			install.packages("BiocManager")
			message("BiocManager has been successfully updated to version ", packageVersion("BiocManager"), ".")
		} else {
			stop("BiocManager update is required to proceed. Exiting.")
		}
	} else {
		message("Your BiocManager version (", bioc_version, ") is up-to-date.")
	}
	
	# Optionally, you can check the Bioconductor version
	message("Current Bioconductor version: ", BiocManager::version())
}

# Run the function
check_and_update_biocmanager()

using<-function(...) {
	libs<-unlist(list(...))
	req<-unlist(lapply(libs,require,character.only=TRUE))
	need<-libs[req==FALSE]
	n<-length(need)
	if(n>0){
		
	}
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("APAlyzer")

using("APAlyzer","Rsamtools","txdbmaker","GenomicFeatures",
			"ggplot2","parallel","zeallot", "biomaRt", "dplyr")
#####################################################
library(APAlyzer)
suppressMessages(library("Rsamtools"))
suppressMessages(library(txdbmaker))
suppressPackageStartupMessages(library(GenomicFeatures))
library(ggplot2)
library(parallel)
library(zeallot)
library(biomaRt)
library(dplyr)

# Function to determine the script's directory (works for Rscript in terminal)
get_script_dir <- function() {
	# For Rscript execution, sys.frame(1)$ofile provides the script's path
	script_path <- tryCatch(normalizePath(sys.frame(1)$ofile),
													error = function(e) NULL)
	if (!is.null(script_path)) {
		return(dirname(script_path))
	} else {
		stop("Cannot determine the script directory. Ensure the script is executed with Rscript.")
	}
}

#### 1. Sample data loading ####
parse_args <- function(){
	args=(commandArgs(TRUE))
	controlCount <- args[1]
	treatmentCount <- args[2]
	data_folder <- args[3]
	
	# if(is.null(data_folder) | !dir.exists(data_folder)){
	# 	script_dir <- getwd()
	# 	data_folder <- file.path(dirname(script_dir), "data_folder")
	# }
	# Check if data_folder exists
	if (!dir.exists(data_folder)) {
		script_dir <- getwd() 
		data_folder <- file.path(dirname(script_dir), "data_folder")
		cat("Data folder not found, using default:", data_folder, "\n")
	}
	return(c(controlCount,treatmentCount, data_folder))
}

# detect whether R script is launched in batch mode or RStudio mode
batch_mode_on <- is.na(Sys.getenv("RSTUDIO", unset = NA))
if(batch_mode_on){
	cat("Running in batch mode\r\n")

	parsed_args <- parse_args()
	treatmentCount <- as.numeric(parsed_args[1])
	controlCount <- as.numeric(parsed_args[2])
	data_folder <- parsed_args[3]
										 
}else{
	cat("Running in iteractive mode\r\n")
	data_folder <- file.path(dirname(getwd()), "data_folder")
	treatmentCount <- 6
	controlCount <- 3
}

# extract vector of new sample names from a file
sample_names_path <- paste(getwd(),"sample_names_metadata.txt", sep = "/")
new_sample_names <- scan(sample_names_path,
										 what = "character",
										 sep = ",")
# checking if number of groups is equal
if(length(new_sample_names) != (controlCount + treatmentCount)){
	stop("Error: number of groups in metadata file does not equal number of
			 groups provided as arguments")
}

# checking if .bam files are located in a folder with the script

cat("Data folder resolved to:", data_folder, "\n")

# Get all .bam files with their absolute paths
file_path <- list.files(path = data_folder, pattern = "\\.bam$", full.names = TRUE)
flsall <- BamFileList(file_path)

# Replace the names in the flsall object
names(flsall) <- new_sample_names

# Check if the renaming is successful
# print(names(flsall))
# flsall


#### 2. Constructing PAS reference from GTF file ####

# read data from single GTF in the current folder
GTFfile_path <- list.files(path = data_folder, pattern = "\\.gtf.gz$", full.names = TRUE)
GTFfile_name <- GTFfile_path[1]

#parse the GTF file
if(batch_mode_on){
	GTFfile <- GTFfile_name
	PASREFraw=PAS2GEF(GTFfile)
}else{
	PASREFraw <- readRDS("PASREFraw.rds")
}

## Converting the parsed GTFREF file to ensmbl format by manually changing chromosome names
PASREFraw <- lapply(PASREFraw, function(df) {
    if ("Chrom" %in% names(df)) {
        # Modify the "Chrom" column in each data frame
        df$Chrom <- gsub("^chr", "", df$Chrom)  # Remove "chr" prefix
    }
    return(df)
})

## Passing parsed data to PASREF list of dataframes
refUTRraw=PASREFraw$refUTRraw
dfIPAraw=PASREFraw$dfIPA
dfLEraw=PASREFraw$dfLE

PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)


## Correcting relevant PASREF dataframe column values to be of numeric type

dfIPA = PASREF$dfIPA
dfLE = PASREF$dfLE
dfIPA$Pos = as.numeric(as.character(dfIPA$Pos))
dfIPA$upstreamSS = as.numeric(as.character(dfIPA$upstreamSS))
dfIPA$downstreamSS = as.numeric(as.character(dfIPA$downstreamSS))
dfLE$LEstart = as.numeric(as.character(dfLE$LEstart))
dfLE$TES = as.numeric(as.character(dfLE$TES))


#### 4. Analysis of APA in 3’UTRs ####
## reprocessing paired-end bam files
# Define the folder name
folder_name <- "ThreemostBamFolder"

# Get the full path
ThreemostBamFolderPath <- file.path(data_folder, folder_name)

# Check if the folder exists
if (!dir.exists(ThreemostBamFolderPath)) {
	dir.create(ThreemostBamFolderPath)  # Create the folder if it doesn't exist
	message("Folder created: ", ThreemostBamFolderPath)
} else {
	message("Folder already exists: ", ThreemostBamFolderPath)
}

# Loop through each BAM file and
# Extract 3 prime most alignment of a paired-end 
# bam file and saved into a new bam file
for (bam_file in file_path) {
	message("Extracting three prime most alignment from
					the paired-end bam file: ", bam_file)
	
	ThreeMostPairBam (BamfilePath=flsall, 
									OutDirPath=ThreemostBamFolderPath, 
									StrandType='forward-reverse')
	}



# Load newly reprocessed bam files
bam_file_path <- list.files(path = ThreemostBamFolderPath,
												pattern = "\\.bam$", full.names = TRUE)
bamProcessed <- BamFileList(file_path)


if(batch_mode_on){
## Building aUTR and cUTR references	
	UTRdbraw=REF3UTR(refUTRraw)
	
## Calculation of relative expression
	DFUTRraw = PASEXP_3UTR(UTRdbraw, bamProcessed, Strandtype="NONE")
	
}else{
	DFUTRraw <- readRDS("DFUTRraw.rds")
}

#### 5. Analysis of APA in introns ####

## Calculating the number of threads to use for computing

#nts means Number of ThreadS used for computing used by featureCounts
# just use max number of cores
numCores <- detectCores()
numThreads <- numCores - 1
paste("Used threads for computing equals:", numThreads)

## Calculation of relative expression of 3’UTR APA and IPA

# workaround so v argument is recognized properly as sample_name: path_to_sample
bamToIntrons <- file_path
names(bamToIntrons) <- new_sample_names
bamToIntrons

if(batch_mode_on){
  IPA_OUTraw = PASEXP_IPA(dfIPA, dfLE, bamProcessed, Strandtype="NONE", nts=numThreads, SeqType='ThreeMostPairEnd')
}else{
	IPA_OUTraw <- readRDS("IPA_OUTraw.rds")
}
#### 6. Significance analysis of APA events ####

treatmentGroup <- "Treatment"
controlGroup <- "Control"

sampleTable = data.frame(samplename = c(new_sample_names),
                    condition = c(rep(treatmentGroup,treatmentCount),
                    							rep(controlGroup,controlCount)))

## Significantly regulated APA in 3’UTRs

# Analysis 3'UTR APA between KD and NT group using non-repilicate design
test_3UTRsing=APAdiff(sampleTable,DFUTRraw,  
                        conKET=controlGroup,
                        trtKEY=treatmentGroup,
                        PAS='3UTR',
                        CUTreads=0,
                        p_adjust_methods="fdr")

## Significantly regulated APA in introns

test_IPAsing=APAdiff(sampleTable,
                        IPA_OUTraw, 
                        conKET=controlGroup,
                        trtKEY=treatmentGroup,
                        PAS='IPA',
                        CUTreads=0,
                        p_adjust_methods="fdr") 



#### 7. Visualization of analysis results ####
## 3'UTR APA plotting
png("3UTRVolcano_plot-Daniel.png", width = 800, height = 600)
APAVolcano(test_3UTRsing, PAS='3UTR', Pcol = "pvalue", top=5, plot_title='3UTR APA')

dev.off()


png("3UTRBox_plot-Daniel.png", width = 800, height = 600)
APABox(test_3UTRsing, xlab = "APAreg", ylab = "RED", plot_title = NULL)

dev.off()


## IPA plotting
png("IPAVolcano_plot-Daniel.png", width = 800, height = 600)
APAVolcano(test_IPAsing, PAS='IPA', Pcol = "pvalue", top=5, plot_title='IPA')

dev.off()

png("IPABox_plot-Daniel.png", width = 800, height = 600)
APABox(test_IPAsing, xlab = "APAreg", ylab = "RED", plot_title = NULL)

dev.off()


## 3'UTR APA vs IPA plotting
### building the plotting data frame for violin plots and CDF curves
test_3UTRsing$APA="3'UTR"
test_IPAsing$APA="IPA"
dfplot=rbind(test_3UTRsing[,c('RED','APA')],test_IPAsing[,c('RED','APA')])

###violin
png("ViolinDiff_plot-Daniel.png", width = 800, height = 600)
ggplot(dfplot, aes(x = APA, y = RED)) + 
    geom_violin(trim = FALSE) + 
    geom_boxplot(width = 0.2)+ theme_bw() + 
    geom_hline(yintercept=0, linetype="dashed", color = "red")


dev.off()


###CDF

png("CDFDiff_plot-Daniel.png", width = 800, height = 600)
ggplot(dfplot, aes( x = RED, color = APA)) + 
    stat_ecdf(geom = "step") +
    ylab("cumulative fraction")+ 
    geom_vline(xintercept=0, linetype="dashed", color = "gray")+ theme_bw() + 
    geom_hline(yintercept=0.5, linetype="dashed", color = "gray")


dev.off()

###########
# Retrieve gene_id and gene_symbol mapping from Ensembl
get_gene_mapping <- function() {
	mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Replace dataset if non-human
	gene_mapping <- getBM(
		attributes = c("ensembl_gene_id", "hgnc_symbol"),
		mart = mart
	)
	colnames(gene_mapping) <- c("Geneid", "gene_symbol")
	return(gene_mapping)
}


#### script below works as intended for adding ensembl gene_id for every row###
# Generate gene mapping table
gene_mapping <- get_gene_mapping()

# Merge results with gene_mapping for 3utr
enriched_result_3UTR <- test_3UTRsing %>%
	left_join(gene_mapping, by = "gene_symbol") %>%  # Match using gene_symbol
	select(Geneid, gene_symbol, RED, pvalue, p_adj, APAreg)

# Merge results with gene_mapping for IPA
enriched_result_IPA <- test_IPAsing %>%
	left_join(gene_mapping, by = "gene_symbol") %>%  # Match using gene_symbol
	select(Geneid, gene_symbol,PASid, RED, pvalue, p_adj, APAreg)

#### Saving results ####

write.csv(enriched_result_3UTR, "3UTR_APAdiff_results.csv",row.names = FALSE, 
					quote = FALSE)
write.csv(enriched_result_IPA, "IPA_APAdiff_results.csv",row.names = FALSE, 
					quote = FALSE)



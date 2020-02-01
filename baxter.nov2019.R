# This script contains useful functions to takes care of formatting
# issues with the metadatafile and mothur outputs from J. nigra microbiome 2017,
# allowing the Excell metadata spreadsheet to be read in correctly
# every time one of its functions are called, and deals with
# issues stemming from inconsistant naming conventions between
# the three individual Illumina runs that were used to produce
# the data.

library(tidyverse)
library(readxl)

get_metadata <- function(){
	
	wd <- getwd()
	
	# Put your directory here
	setwd(".")

	# Recode metadata variables

    metadata <- read_excel(path="baxter.metadata.names.xlsx",
    col_types=c(Illumina_Cell_Code = "numeric", Sample_Well = "text", Group = "text", Sample_ID2 = "text", F= "text", R = "text", Date = "text", State = "text", County = "text", Type = "text", Clone = "text", TreeTCDLevel = "text", TreeSize = "text", TreeAge = "text", Soil_Series = "text", Soil_Texture = "text", SoilType = "text", soilpH_rep1 = "numeric", soilEC_rep1 = "numeric", soilpH_rep2 = "numeric", soilEC_rep2 = "numeric", HerbicideUse = "text", Groundcover = "text", New.names= "text" ) )
    metadata[["Soil_Series"]] <- na_if(metadata[["Soil_Series"]], "?")
    metadata[["Soil_Texture"]] <- na_if(metadata[["Soil_Texture"]], "?")

    # Reformat rownames so they match whats in the mothur output files and put them in "sample"
    
    metadata$sample <- gsub(".fastq", "", metadata$F) %>% (function (x) gsub("-", "_", x)) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))

	setwd(wd)

    return(metadata)

}

get_metadata_bacteria <- function() {
	
	wd <- getwd()
	
	# Put your directory here
	setwd(".")
	
	# Recode metadata variables

    metadata <- read_excel(path="baxter.metadata.merged.Aaron.xlsx",
    col_types=c(Illumina_Cell_Code = "numeric", Sample_Well = "text", Group = "text", Sample_ID2 = "text", F= "text", R = "text", Date = "text", State = "text", County = "text", Type = "text", Clone = "text", TreeTCDLevel = "text", TreeSize = "text", TreeAge = "text", Soil_Series = "text", Soil_Texture = "text", SoilType = "text", soilpH_rep1 = "numeric", soilEC_rep1 = "numeric", soilpH_rep2 = "numeric", soilEC_rep2 = "numeric", HerbicideUse = "text", Groundcover = "text", New.names= "text", Bact.names = "text", Bact.names.2 = "text") )
    metadata[["Soil_Series"]] <- na_if(metadata[["Soil_Series"]], "?")
    metadata[["Soil_Texture"]] <- na_if(metadata[["Soil_Texture"]], "?")
    metadata[["SoilType"]] <- recode(.x=metadata[["SoilType"]], "Ultisol - Paleudult"="ultisol")
    metadata[["SoilType"]] <- recode(.x=metadata[["SoilType"]], "Ultisol - Hapludult"="ultisol")
    metadata[["SoilType"]] <- recode(.x=metadata[["SoilType"]], "Mollisol - Palouse SiltLoam"="mollisol")
    metadata[["County"]] <- recode(.x=metadata[["County"]], "Tippacanoe"="Tippecanoe")

    # Reformat rownames so they match group names in the fungal mothur output files
    
    metadata$sample <- gsub(".fastq", "", metadata$F) %>% (function (x) gsub("-", "_", x)) %>% (function (x) gsub("_L001", "", x)) %>% (function (x) gsub("_001", "", x)) %>% (function (x) gsub("_R1", "", x))

	setwd(wd)

    return(metadata)

}


get_metadata_soilmb <- function() {
	
	wd <- getwd()
	
	# Put your directory here
	setwd(".")
	
	# Recode metadata variables
    
    metadata <- read_excel(path="baxter.metadata.names.xlsx",
    col_types=c(Illumina_Cell_Code = "numeric", Sample_Well = "text", Group = "text", Sample_ID2 = "text", F = "text", R = "text", Date = "text", State = "text", County = "text", Type = "text", Clone = "text", TreeTCDLevel = "text", TreeSize = "text", TreeAge = "text", Soil_Series = "text", Soil_Texture = "text", SoilType = "text", soilpH_rep1 = "numeric", soilEC_rep1 = "numeric", soilpH_rep2 = "numeric", soilEC_rep2 = "numeric", HerbicideUse = "text", Groundcover = "text" ,New.names="text") )
    metadata[["Soil_Series"]] <- na_if(metadata[["Soil_Series"]], "?")
    metadata[["Soil_Texture"]] <- na_if(metadata[["Soil_Texture"]], "?")
    
    # Correct sample names were hardcoded into the Excell file in column New.names
    metadata$sample <- metadata$New.names 
  
    columns_to_remove <- which(colnames(metadata) %in% c('Illumina_Cell_Code', 'Sample_Well', 'Sample_Well', 'F', 'R', 'Date', 'Type'))
    
    # Get rid of unwnated duplicate rows
    
    rows_to_remove <- grep('s$', perl=T, metadata$Group)
    
    setwd(wd)
    
    return(metadata[-rows_to_remove, -columns_to_remove])
    
}

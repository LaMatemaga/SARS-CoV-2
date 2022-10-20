# Main.R
# Core code of the research
# Authors: Cynthia Castillo, Joaquin Lopez, Carolina Sanmiguel,
#          Javier Almaguer & Francisco Cabrera


#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#
#BiocManager::install("Pbase")
#

# Initialize libraries, workspace, and functions
packages <- c("seqinr",               # Reading/writing fasta files
              "expm",                 # Matrix to the nth power
              "dplyr",                # Manipulate dataframes
              "plyr",                 # Combine dataframes
              "readr",                # read_csv is way faster than read.csv
              "lubridate",            # Work with dates
              "tibble",               # Auxiliary library for dplyr
              "ttutils")              # Merge is used to merge lists
              #"DECIPHER")            # To obtain consensus
              #"Pbase")               # Proteomic things, I think we are not using this
              #"crayon")              # Concatenate char vectors for mutations
              #"reticulate")          # Run Python on R
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)
rm(packages)
setwd("C:/Users/Elaia/Documents/Academia/Investigación/SARS-CoV-2")
source("Toolkit.R")                   # Call custom functions


# File merging
# WARNING: Prior manual cleaning to metadata need to be done before running
#          lines from 34 to 42.
#          Variable 'metadataToMerge' is expected to contain only clean data.
sequencesToMerge <- c("Sequences2020Jul09.fasta","Sequences2020Jul18.fasta",
                      "Sequences2021Mar22.fasta","Sequences2021May14.fasta",
                      "Sequences2021Oct27.fasta","Sequences2022Feb16.fasta",
                      "Sequences2022Feb18.fasta","Sequences2022Feb21.fasta",
			                "Sequences2022Mar23.fasta","Sequences2022Abr11.Fasta")
metadataToMerge  <- c("Metadata2020Jul09.csv","Metadata2020Nov29.csv",
                      "Metadata2021Mar22.csv","Metadata2021May18.csv",
                      "Metadata2021Oct27.csv","Metadata2022Feb16.csv",
                      "Metadata2022Feb18.csv","Metadata2022Feb21.csv",
			                "Metadata2022Mar23.csv","Metadata2022Abr11.csv")
mergeSequences(sequencesToMerge,metadataToMerge)


# Set-up input files
file <- "input/Sequences.fasta"
sequences <- read.fasta(file, seqtype = "AA", forceDNAtolower = FALSE,
                        set.attributes = FALSE)
metadata  <- read.csv(file="input/Sequences.csv")
metadata  <- column_to_rownames(metadata,var="X.1")
metadata  <- metadata %>% select(Release_Date,Geo_Location,Collection_Date)
rm(file)


# Modify Collection_Date and Release_Date so they match dates
metadata$Release_Date    <- mdy(metadata$Release_Date)     #Corregir porque convierte a NA XD formato incorrecto
metadata$Collection_Date <- mdy(metadata$Collection_Date)


# Clean ambiguous data (excess of 'X')
lengthSeq <- 1273
output    <- cleanData(sequences,metadata,lengthSeq)
sequences <- output[[1]]
metadata  <- output[[2]]
rm(output)


# Save cleaned data
write.fasta(sequences,names(sequences),file.out="input/CleanedSequences.fasta")
write.csv(metadata,file="input/CleanedSequences.csv")


# Obtain consensus of our data and calculate the metric we are going to use
aminoacids        <- c("A","R","N","D","C","Q","E","G","H","I",
                       "L","K","M","F","P","S","T","W","Y","V")
consensusSequence <- getConsensus(sequences,aminoacids)
metricMatrix      <- getMetricMatrix(sequences,consensusSequence,aminoacids)


# Get mutated positions
mutationPositions <- checkMutations(consensusSequence)


#####################################################
# Check different countries in the data base
write(sort(unique(metadata$Geo_Location)),file="auxiliarFiles/Countries.txt")


# Should take care of this during first cleaning, meanwhile use these lines
metadata$Geo_Location[metadata$Geo_Location=="USA: Puerto Rico"] <- "Puerto Rico"
metadata$Geo_Location[metadata$Geo_Location=="USA: Guam"] <- "Guam"


# Set filters per region
fileFilters <- c("AsiaEast","AsiaWest","Europe","Oceania","SouthAmerica",
                 "USA","USA-California","USA-Center","USA-Florida",
                 "USA-Massachusetts","USA-Michigan","USA-Minnesota","USA-North",
                 "USA-NorthEast","USA-Oregon","USA-Others","USA-South","USA-SouthEast",
                 "USA-Virginia","USA-Washington","USA-Wisconsin")
fileFilters <- paste("input/",fileFilters,".txt",sep="")
infoRegions <- filterByRegion(metadata,fileFilters)
metadata    <- infoRegions$metadata
infoRegions$metadata <- NULL

#####################################################


# Transform sequences using the metric matrix and filtrate by regions
# uniqueSeqFilter <- filterUnique(sequences,metadata)
transformWithMetric(sequences,metadata,metricMatrix,consensusSequence,aminoacids,
                    mutationPositions)
filtrateByRegions(fileFilters)
getMutations(sequences,metadata,consensusSequence)

# !! Bug: Se guardan los encabezados como filas xd
# !! Bug: Se generan duplicados
# !! Bug: Se guardan las columnas como .1 y .2, buscar manera de que no se guarden asi
# globalMetricTransformation <- unique(globalMetricTransformation)


save.image(file="temporal/Workspace 11Abr2022.RData")
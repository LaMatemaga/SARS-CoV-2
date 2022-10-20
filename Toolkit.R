# Toolkit.R
# Auxiliary functions made for this research
# Authors: Cynthia Castillo, Joaquin Lopez, Carolina Sanmiguel,
#          Javier Almaguer & Francisco Cabrera

mergeMetadata <- function(metadataToMerge){
  # Auxiliar function for mergeSequences which only merges the metadata.
  
  #Set variables
  filesMet <- length(metadataToMerge)
  
  # Initialize sequences' metadata
  mergedMetadata  <- read.csv(file=paste("merge/",metadataToMerge[1],sep=""))
  names(mergedMetadata)[1] <- "Accession"
  mergedMetadata  <- mergedMetadata %>% select(Accession,Release_Date,Geo_Location,Collection_Date)

  
  # Merge next sequences' metadata with the initialized one
  for(file in 2:filesMet){
    # Read next sequences' metadata.
    metadata       <- read.csv(file=paste("merge/",metadataToMerge[file],sep=""))
    names(metadata)[1] <- "Accession"
    metadata       <- metadata %>% select(Accession,Release_Date,Geo_Location,Collection_Date)

    
    # Join new dataframe with previous one
    # Then avoid duplicates using Accession variable
    mergedMetadata <- union(metadata,mergedMetadata)
    mergedMetadata <- distinct(mergedMetadata,Accession,.keep_all=TRUE)
    rm(metadata)
  }
  
  # Write mergedMetadata for further use
  write.csv(mergedMetadata,file="merge/MergedMetadata.csv")
}

mergeFasta <- function(sequencesToMerge){
  # Auxiliar function for mergeSequences which only merges the fasta files.
  
  #Set variables
  filesSeq <- length(sequencesToMerge)
  
  # Initialize sequences
  mergedSequences <- read.fasta(paste("merge/",sequencesToMerge[1],sep=""),
                                seqtype = "AA", forceDNAtolower = FALSE,
                                as.string=TRUE, set.attributes = FALSE)
  
  # Rename sequences so they match the metadata Accession ID
  dataNames <- names(mergedSequences)
  filterDBN <- rep(c(T,F),length(mergedSequences))
  dataNames <- unlist(strsplit(dataNames,split = "[.]"))[filterDBN]
  names(mergedSequences) <- dataNames
  
  # Merge next sequences with the initialized one
  for(file in 2:filesSeq){
    # Read next sequences
    sequences       <- read.fasta(paste("merge/",sequencesToMerge[file],sep=""),
                                  seqtype = "AA", forceDNAtolower = FALSE,
                                  as.string=TRUE, set.attributes = FALSE)
    # Rename sequences so they match the metadata Accession ID
    dataNames <- names(sequences)
    filterDBN <- rep(c(T,F),length(sequences))
    dataNames <- unlist(strsplit(dataNames,split = "[.]"))[filterDBN]
    names(sequences) <- dataNames
    
    # Combine lists
    mergedSequences <- merge(sequences,mergedSequences,mergeUnnamed=FALSE)
    rm(sequences)
  }
  
  # Write mergedSequences for further use
  write.fasta(mergedSequences,names(mergedSequences),file.out="merge/MergedSequences.fasta")
}

mergeSequences  <- function(sequencesToMerge,metadataToMerge){
  # Merge all the sequences and their metadata in one single file, omitting
  # duplicates and missing information (i.e. missing sequence or missing metadata)
  
  # Set variables
  mergeFasta(sequencesToMerge)
  mergeMetadata(metadataToMerge)
  
  # Initialize the sequences
  mergedSequences <- read.fasta(file="merge/MergedSequences.fasta",
                                seqtype = "AA", forceDNAtolower = FALSE,
                                as.string=TRUE, set.attributes = FALSE)
  mergedMetadata  <- read.csv(file="merge/MergedMetadata.csv")
  mergedMetadata  <- column_to_rownames(mergedMetadata,var="Accession")
  
  # Check the intersection between metadata and sequences
  namesSequences <- names(mergedSequences)
  namesMetadata  <- row.names(mergedMetadata)
  intersection   <- intersect(namesSequences,namesMetadata)
  
  # Save only data on the intersection for both files
  mergedSequences <- mergedSequences[intersection]
  mergedMetadata  <- mergedMetadata[intersection,]
  
  # Save files in /input folder
  write.fasta(mergedSequences,names(mergedSequences),file.out="input/Sequences.fasta")
  write.csv(mergedMetadata,file="input/Sequences.csv")
}

cleanData       <- function(sequences,metadata,lengthSeq,ambiguityTolerance=0.02){
  # Sorts all data by Accession ID and cleans ambiguous data
  # over ambiguityTolerance.
  
  # Set variables
  maxAmbiguity <- floor(ambiguityTolerance*lengthSeq)
  sizeDB <- as.integer(length(sequences))
  
  # Rename sequences so they match the metadata Accession ID
  dataNames <- row.names(metadata)
  metadata <- rownames_to_column(metadata,var="Accession")
  
  # Reorder sequences by dataNames so the next loop will be faster
  sequences <- sequences[order(dataNames)]
  metadata  <- metadata[order(dataNames),]
  
  # Check if there's any data missing from the metadata
  # No longer necessary because it was done on preprocessing and cleaning
  #if(nrow(metadata)>sizeDB | any(!(dataNames %in% row.names(metadata)))){
  #  errorCondition("Entries on both sequences and metadata do not match.")
  #  return()
  #}
  
  # Loop for cleaning ambiguous data over ambiguityTolerance
  seq <- 1                            # Number of samples
  while(seq<=sizeDB){
    aux <- 0
    if(length(sequences[[seq]])!=1273){
      # En caso de que la longitud no sea de 1273
      metadata     <- metadata[-seq,]
      sequences[seq] <- NULL
      sizeDB       <- sizeDB-1
    }else{
      for(AA in 1:lengthSeq){
        if(sequences[[seq]][AA]=='X'){
          aux <- aux+1
        }
        if(aux>maxAmbiguity){
          # Aqui me quita el dato "ambiguo"
          metadata     <- metadata[-seq,]
          sequences[seq] <- NULL
          sizeDB       <- sizeDB-1
          break
        }
      }
    }
    seq <- seq+1
  }
  output <- list(sequences, metadata)
  return(output)
}

filterByRegion  <- function(metadata,fileFilters){
  # Gives a filter of each one of the regions defined by the user.
  
  # Set variables
  numFilters <- length(fileFilters)
  regions    <- gsub("input/|.txt", "", fileFilters)
  regionID    <- c()
  
  # Loop to get the filters for each region
  filters    <- list()
  for(fil in 1:numFilters){
    filters[[regions[fil]]] <- list(filter=c())
    countries      <- scan(file = fileFilters[fil], what = character(),
                           flush = TRUE, sep = "\n")
    filters[[fil]]$filter <- metadata$Geo_Location %in% countries
    regionID[filters[[fil]]$filter] <- fil
  }
  
  # Add new column to classify the metadata
  metadata$RegionID <- regionID
  filters$metadata <- metadata
  
  # Return the filters and the modified metadata. Metadata should be deleted
  # from the list later.
  return(filters)
}

getConsensus    <- function(sequences,aminoacids){
  # Obtains the consensus of all the sequences we have in our database.
  
  # Set variables
  sizeDB         <- as.integer(length(sequences))
  lengthSeq      <- as.integer(length(sequences[[1]]))
  numAminoacids  <- length(aminoacids)
  auxiliarMatrix <- matrix(0, lengthSeq, numAminoacids)
  colnames(auxiliarMatrix) <- aminoacids
  
  # Loop to count which aminoacids are the most common for each aminoacidic
  # position
  for(sample in 1:sizeDB){
    for(location in 1:lengthSeq){
      posAmino <- which(aminoacids==sequences[[sample]][location])
      auxiliarMatrix[location,posAmino] <- auxiliarMatrix[location,posAmino]+1
    }
  }
  
  # Put the 'X' on the overall consensus
  for(amino in 1:lengthSeq){       # Para poner las X en el consenso manualmente
    aux <- sizeDB-sum(auxiliarMatrix[amino,])
    index <- which(auxiliarMatrix[amino,]==max(auxiliarMatrix[amino,]))
    auxiliarMatrix[amino,index] <- aux + auxiliarMatrix[amino,index]
  }
  
  # Get the consensus sequence using the auxiliar matrix
  consensusSequence <- c()
  for(amino in 1:lengthSeq){
    index <- which(auxiliarMatrix[amino,]==max(auxiliarMatrix[amino,]))
    consensusSequence <- c(consensusSequence,aminoacids[index])
  }
  
  # Save the auxiliar matrix and consensus for further use, just in case
  write.csv(auxiliarMatrix, file = "temporal/auxiliarMatrixForConsensus.csv")
  write.fasta(consensusSequence, names="consensus",file = "temporal/consensus.fa")
  
  # Return consensus sequence
  return(consensusSequence)
}

getMetricMatrix <- function(sequences,consensusSequence,aminoacids){
  # Gets the metric matrix we need to perform both TDA and Random Forest
  
  # Set variables
  sizeDB         <- as.integer(length(sequences))
  lengthSeq      <- as.integer(length(sequences[[1]]))
  numAminoacids  <- length(aminoacids)
  freqMatrix     <- matrix(0, nrow=numAminoacids, ncol=numAminoacids)
  colnames(freqMatrix) <- aminoacids
  rownames(freqMatrix) <- aminoacids
  
  # Get frequencies on what AA do we have by sample vs what AA should be there
  for (sample in 1:sizeDB){
    for (amino in 1:lengthSeq){
      if(sequences[[sample]][amino] %in% aminoacids){
        Is <- which(aminoacids == sequences[[sample]][amino])
      }else{
        Is <- which(aminoacids == consensusSequence[amino])
      }
      ShouldBe <- which(aminoacids == consensusSequence[amino])
      freqMatrix[Is,ShouldBe] <- freqMatrix[Is,ShouldBe] + 1
    }
  }
  
  # Calculate the metric matrix we will use for the analysis
  sumFrequencies <- rowSums(freqMatrix)
  relativeMatrix <- t(t(freqMatrix+1)/sumFrequencies)
  metricMatrix   <- -log(relativeMatrix %^% 2, 2)    # Based on PAM2
  
  # Adjust the values so they go from the range 0 - 50
  metricMatrix   <- metricMatrix/max(metricMatrix) * 50
  
  # Adjust the diagonal to 0
  for (i in 1:numAminoacids){
    metricMatrix[i,i] <- 0
  }
  
  # Write temporal files, in case we need them
  write.csv(metricMatrix, file="temporal/metricMatrix.csv")
  write.csv(freqMatrix, file="temporal/freqMatrix.csv")
  
  # Return the metric matrix, as we need it for the next steps
  return(metricMatrix)
}

filterUnique    <- function(sequences,metadata){
  # Obtain the first unique mutations per country and state
  
  # Set variables
  sizeDB     <- as.integer(length(sequences))
  lengthSeq  <- as.integer(length(sequences[[1]]))
  auxDates   <- as.integer(metadata$Collection_Date)
  
  # Get auxiliary string to compare
  auxSequences <- c()
  for(sample in 1:sizeDB){
    auxSequences[sample] <- paste(sequences[[sample]],  collapse= '')
  }
  auxSequences <- paste(auxSequences,metadata$Geo_Location,metadata$USA,sep="")
  
  # Get the first ocurrence of a sample per country and state
  checkedSamples <- as.vector(matrix(FALSE, nrow = sizeDB))
  uniqueSeqID    <- c()
  for (sample in 1:sizeDB){
    if(!checkedSamples[sample]){
      filterSimilar   <- auxSequences==auxSequences[sample]
      firstOccurence  <- min(auxDates[filterSimilar])
      similarAndFirst <- filterSimilar & auxDates==firstOccurence
      getID           <- which(similarAndFirst==TRUE)[1]
      uniqueSeqID    <- c(uniqueSeqID, getID)
      checkedSamples <- checkedSamples==TRUE | filterSimilar
    }
  }
  
  # Organize the data
  uniqueSeqID <- uniqueSeqID[order(uniqueSeqID)]
  write(uniqueSeqID, file = "temporal/UniqueSequencesID.txt")
  
  # Get filter
  uniqueSeqFilter <- as.vector(matrix(FALSE, nrow = sizeDB))
  for(sample in 1:sizeDB){
    if(sample %in% uniqueSeqID){
      uniqueSeqFilter[sample] <- TRUE
    }
  }
  
  # Return unique sequence filter
  return(uniqueSeqFilter)
}

checkMutations  <- function(consensusSequence){
  # Gets the AA positions that are changed in respect with the consensusSequence.
  # Uses the auxiliar matrix we generated before.
  
  # Set variables
  auxiliarMatrix <- read.csv(file = "temporal/auxiliarMatrixForConsensus.csv")
  auxiliarMatrix <- auxiliarMatrix %>% select(-X)
  lengthSeq      <- length(consensusSequence)
  mutationsPerPosition <- c()
  
  # Check which AA positions have a mutation and how many
  for (amino in 1:lengthSeq) {
    mutationsPerPosition <- c(mutationsPerPosition,
                              sum(auxiliarMatrix[amino,])-max(auxiliarMatrix[amino,]))
  }
  write(which(mutationsPerPosition!=0), file = "temporal/changedAminoacids.txt")
  write(mutationsPerPosition, file = "temporal/mutationQuantityPerPosition.txt")
  return(which(mutationsPerPosition!=0))
}

transformToVector <- function(startPoint,endPoint,mutationPositions,aminoacids,
                              consensusSequence,metricMatrix,sequences){
  # Function that helps to fragment all the process in different segments
  # so we don't run out of memory that easily
  
  # Set variables
  nRows <- endPoint-startPoint+1
  metricTransformation <- matrix(0, nrow=nRows, ncol=length(mutationPositions))
  colnames(metricTransformation) <- paste0(consensusSequence[mutationPositions],mutationPositions)
  
  # Vectorize
  for (seq in 1:nRows){
    for (AA in 1:length(mutationPositions)){
      ShouldBe <- consensusSequence[mutationPositions[AA]]
      if(sequences[[seq]][mutationPositions[AA]]=="B"){
        # D o N
        if(metricMatrix["D",ShouldBe] < metricMatrix["N",ShouldBe]){
          Is <- "D"
        }else{
          Is <- "N"
        }
      }else if(sequences[[seq]][mutationPositions[AA]]=="Z"){
        # E o Q
        if(metricMatrix["E",ShouldBe] < metricMatrix["Q",ShouldBe]){
          Is <- "E"
        }else{
          Is <- "Q"
        }    
      }else if(sequences[[seq]][mutationPositions[AA]]=="J"){
        # I o L
        if(metricMatrix["I",ShouldBe] < metricMatrix["L",ShouldBe]){
          Is <- "I"
        }else{
          Is <- "L"
        }
        #}else if(sequences[[seq]][mutationPositions[AA]]=="U"){
        # Selenocysteine
        #Is <- sequences[[seq]][mutationPositions[AA]]
        #}else if(sequences[[seq]][mutationPositions[AA]]=="O"){
        # Pyrrolysine
        #Is <- sequences[[seq]][mutationPositions[AA]]
      }else if(sequences[[seq]][mutationPositions[AA]]=="X"){
        # Unknown
        Is <- consensusSequence[mutationPositions[AA]]
      }else{
        Is <- sequences[[seq]][mutationPositions[AA]]
      }
      metricTransformation[seq,AA] <- metricMatrix[Is,ShouldBe]
    }
  }
  
  # Return vectorized variable
  return(metricTransformation)
}

transformWithMetric <- function(sequences,metadata,metricMatrix,consensusSequence,
                                aminoacids,mutationPositions){
  # Transform the sequences to numeric values using the metricMatrix.
  
  # Set variables
  sizeDB         <- as.integer(length(sequences))
  lengthSeq      <- as.integer(length(sequences[[1]]))
  metricTransformation <- matrix(0, nrow=0, ncol=length(mutationPositions))
  colnames(metricTransformation) <- as.character(mutationPositions)
  
  # Vectorize using the metric matrix
  aux <- 0
  while(aux<sizeDB){
    startPoint <- aux+1
    if(startPoint+2500<=sizeDB){
      endPoint <- aux+2500
      aux      <- aux+2500
    }else{
      endPoint <- sizeDB
      aux      <- sizeDB
    }
    
    # Transform data by blocks of 2500 sequences
    metricTransformation <- transformToVector(startPoint,endPoint,mutationPositions,
                                   aminoacids,consensusSequence,metricMatrix,
                                   sequences[startPoint:endPoint])
    rownames(metricTransformation) <- metadata$Accession[startPoint:endPoint]
    
    # Transform previous matrix into a dataframe
    metricTransformation <- as.data.frame(metricTransformation)
    
    # Join with metadata
    metricTransformation <- rownames_to_column(metricTransformation,var="Accession")
    auxMatrix <- metadata[startPoint:endPoint,] %>% select(Accession,Release_Date,
                                                           Geo_Location,RegionID,
                                                           Collection_Date)
    output    <- inner_join(auxMatrix,metricTransformation)
    output    <- column_to_rownames(output,var="Accession")
    
    # Write output in a file
    write.csv(output, file=paste("temporal/ToMerge/GlobalMetricTransformation",
                                 sprintf("%03d",startPoint%/%2500),".csv",sep=""))
    
    # Erase variables
    rm(auxMatrix,metricTransformation,output)
  }
  
  # Merge all files into a single file
  mergeFiles(filePath="temporal/ToMerge/")
}

mergeFiles <- function(filePath){
  # Takes all CSVs on a folder and merges them into a single file
  # Adapted from https://statisticsglobe.com/merge-csv-files-in-r
  
  # Identify all csv files in folder
  output <- list.files(path=filePath, pattern="*.csv", full.names=TRUE) %>%
            lapply(read_csv) %>%        # Store all files in list
            bind_rows                   # Combine data sets into one data set 
  
  # Write a single file with all the merged data
  write.csv(output, file="temporal/MetricTransformation-Global.csv")
}

filtrateByRegions <- function(fileFilters){
  # Takes GlobalMetricTransformation.csv and splits it by region in different files.
  
  # Set variables
  numFilters <- length(fileFilters)
  regions    <- gsub("input/|.txt", "", fileFilters)

  # Read file
  globalMetricTransformation <- read_csv(file="temporal/MetricTransformation-Global.csv")
  
  # Filter by region and save the filtered data
  for(region in 1:numFilters){
    if(region==6){
      # In case it's from the USA, consider regions from 6 to 21
      filteredData <- filter(globalMetricTransformation, RegionID>=6)
    }else{
      # Else, consider individual regions
      filteredData <- filter(globalMetricTransformation, RegionID==region)
    }
    
    # Write filtered data
    write.csv(filteredData, file=paste0("temporal/MetricTransformation-",regions[region],".csv"))
  }
}

getMutations <- function(sequences,metadata,consensusSequence,start=319,end=541){
  # Function that gets mutations sets and their frequency according to the
  # consensusSequence. Region of interest is from AA 319 to 541.
  
  # Set variables
  sizeDB     <- as.integer(length(sequences))
  mutations  <- list()
  timestamps <- metadata$Collection_Date
  geographic <- metadata$Geo_Location
  sizeMut    <- 0
  rm(metadata)
  
  # Get mutations sets
  for(seq in 1:sizeDB){
    # Set variables
    mutatedAminoacid <- c()
    mutatedPosition  <- c()
    
    # Check mutations in sequence # seq
    for(AA in start:end){
      if(sequences[[seq]][AA]!=consensusSequence[AA] && sequences[[seq]][AA]!="X"){
        mutatedAminoacid <- c(mutatedAminoacid,sequences[[seq]][AA])
        mutatedPosition  <- c(mutatedPosition,AA)
      }
    }
    
    # Store mutations for sequence # seq
    if(length(mutatedPosition)!=0){
      fingerprint <- paste(paste0(consensusSequence[mutatedPosition],mutatedPosition,mutatedAminoacid),collapse=" ")
    }else{
      fingerprint <- "No mutation"
    }
    
    # Store new mutation in list or increase the frequency in case same
    # mutation set has been identified
    if(sizeMut!=0){
      # Set variables
      flag <- sizeMut+1
      
      # Check if mutation set has been identified before
      for(i in 1:sizeMut){
        if(fingerprint==mutations$fingerprint[i]){
          flag <- i
          break
        }
      }
      
      # Check if mutation set is new or if it's been identified before
      if(flag>sizeMut){
        # Add new mutation set
        mutations$frequency[flag]     <- 1
        mutations$fingerprint[flag]   <- fingerprint
        mutations$firstLocation[flag] <- geographic[seq]
        mutations$firstDate[flag]     <- timestamps[seq]
        sizeMut <- sizeMut+1
      }else{
        # Add 1 to the frequency
        mutations$frequency[flag] <- mutations$frequency[flag]+1
        
        # Check if this mutation appeared before
        if(!is.na(timestamps[seq])){
          if(is.na(mutations$firstDate[flag]) | timestamps[seq]<mutations$firstDate[flag]){
            mutations$firstDate[flag]     <- timestamps[seq]
            mutations$firstLocation[flag] <- geographic[seq]
          }
        }
      }
    }else{
      # Add first mutation set
      mutations$frequency[1]     <- 1
      mutations$fingerprint[1]   <- fingerprint
      mutations$firstLocation[1] <- geographic[seq]
      mutations$firstDate        <- mdy(mutations$firstDate)
      mutations$firstDate[1]     <- timestamps[seq]
      sizeMut <- 1
    }
  }
  
  # Store mutations sets and their frequency in a dataframe and store them in a csv
  output <- as.data.frame(mutations)
  write.csv(output, file="docking/Mutations.csv")
}

similarityScore <- function(A, B){
  
  A <- gsub( "-", "_", tolower(A))
  B <- gsub( "-", "_", tolower(B))
  
  return(stringdist(A,B))
}


findSimilar <- function(gtList, ptList){
  
  gtptSimList <- data.frame(phenotype = ptList)
  gtptSimList$genoType <- ""
  #gtSimList <- array("", dim = length(ptList))
  for (i in 1:length(ptList)){
    minDist <- length(ptList[i])
    for (j in 1:length(gtList)){
      tmpDist <- similarityScore(ptList[i],gtList[j])
      if (tmpDist < minDist){
        minDist <- tmpDist
        gtptSimList[i,"genoType"] <- gtList[j]
      }
    }
  }
  return(gtptSimList)
}
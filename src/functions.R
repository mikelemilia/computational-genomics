ma <- function(x, y){
  m <- log2(x) - log2(y)
  a <- (log2(x) + log2(y))/2
  
  return(list(m,a))
}

tmm_normalization <- function(DATA, refsample){
  refsample <- as.numeric(refsample)
  ref <- DATA[,refsample] + 1
  A <- matrix(0, nrow = length(ref), ncol = 1)
  M <- matrix(0, nrow = length(ref), ncol = 1)
  norm <- matrix(ref, nrow = dim(DATA)[1], ncol = 1)
  
  for (i in (3:ncol(DATA[-1])))
  {
    temp <- DATA[i]+1
    temp <- unlist(temp,use.names=FALSE)
    ref <- unlist(ref,use.names=FALSE)
    result <- ma(ref,temp)
    
    SF <- mean(result[[1]],trim=0.1)
    modtemp <- temp+2^(SF)
    
    result2 <- ma(ref,modtemp)
    
    norm <- cbind(norm,modtemp)
    A <- cbind(A,result2[[2]])
    M <- cbind(M,result2[[1]])
  }
  
  return(list(M,A,norm))
}

quantile_normalization <- function(DATA){
  cols <- ncol(DATA[]) # get number of columns
  rows <- nrow(DATA[]) # get number of rows
  dataSort = matrix(0, rows, cols)
  dataIdx  = matrix(0, rows, cols)
  dataNorm = matrix(0, rows, cols)
  
  # setting the header of dataNorm
  colnames(dataNorm) <- names(DATA[])
  
  for(i in 1:cols){ 
    data = DATA[,i]
    dataSorted = sort(data)
    dataSortedIdxs = rank(data, ties.method="average")
    dataSort[,i] = dataSorted
    dataIdx[,i] = dataSortedIdxs
  }
  
  dataMean = apply(dataSort, 1, mean) 
  
  for(i in 1:cols ) { 
    for(k in 1:rows ) { 
      dataNorm[k,i]= dataMean[dataIdx[k,i]]
    }
  }
  
  return(dataNorm)
}

estimateG0<-function(c_pvalue) {
  lambda<-seq(0, 1, 0.01)
  i<-1
  G0_prev<-0
  G0<-vector("integer",length(length(lambda)))
  r<-vector("integer",length(length(lambda)))
  convergenza<-0
  while (convergenza==0 & i<length(lambda)) {
    
    l<-lambda[i]
    
    minori<-(c_pvalue<l )
    selected<-which(minori==TRUE)
    num_sel<-length(selected)

    G0[i]<-(G-num_sel)/(1-l)
    
    r[i]<-(G0[i]-G0_prev)^2
    
    G0_prev<-G0[i]
    
    i<-i+1
  }
  
  G0_estimate<-G0[which.min(r)]
  return(list(G0_estimate,lambda[which.min(r)]))
}

expected_values <- function(G,G0,alpha,selected_genes){
  expected_FP <- min(G0*alpha, selected_genes)
  expected_TP <- max(0, (selected_genes - expected_FP))
  expected_TN <- G0 - expected_FP
  expected_FN <- max(0,G-selected_genes-expected_TN)

  return (c(expected_TP,expected_FP,expected_TN,expected_FN))
}


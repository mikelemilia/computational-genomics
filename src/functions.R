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
  
  for (i in 3:ncol(DATA))
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

#function that removes duplicated in the 'individual' column
remove_duplicates <- function (data, group){
  #TAKING CARE OF DUPLICATES FOR GROUP
  duplicates <- duplicated(group$individual) # get the position of all the duplicated individual 
  #get only the duplicated
  duplicate_group <- group$individual[duplicates == TRUE]
  #find indexes in group of duplicated subjects
  d <- as.numeric(group$individual)
  
  dim <- dim(group)
  group_nodup <- matrix(0,nrow(data),dim[1])
  #initialize a count because matrix at the end WON'T have same columns as dim[1]=41, the group SRR
  count<-0
  for (i in 1:length(duplicate_group))
  {
    indexes <- which(d %in% duplicate_group[i])
    #find samples of duplicated subjects, get the SRR code
    seq_sample <- group[indexes,]$seq.sample
    #mean of duplicated samples
    group_nodup[,i] <- apply(data[, seq_sample], 1, mean)
    count <- count+1;
  }
  
  #finding subjects NOT DUPLICATED
  '%notin%' <- Negate(`%in%`)
  ind<-which(d %notin% duplicate_group)
  
  for (i in 1:length(ind))
  {
    
    #find samples of duplicated subjects, get the SRR code
    seq_sample<-group[ind[i],]$seq.sample
    #mean of duplicated samples
    group_nodup[,length(duplicate_group)+i]<-data[,seq_sample]
    count<-count+1;
  }
  
  group_nodup<-group_nodup[,1:count]
  return(group_nodup)
}

#functions to remove rows where there are only 0 in the two groups (matrixes)
remove_zeros <- function(group1,group2){
  group1_forcomparison<-apply(group1, 1, sum)
  group2_forcomparison<-apply(group2, 1, sum)
  vec_for_comparison<-as.data.frame(group1_forcomparison+group2_forcomparison)
  #usare i quantili PER TROVARE SOGLIA 
  q<-quantile(unlist(vec_for_comparison))
  #median<-median(unlist(vec_for_comparison)) USANDO MEDIANA? TROPPO POCHI NE RESTANO
  YN<-as.data.frame(vec_for_comparison>q[1])
  index<-which(YN==FALSE) #indici di quelli da togliere
  group1_nozero<-(control_nodup[-index,])
  group2_nozero<-disease_nodup[-index, ]
  
  return(list(group1_nozero,group2_nozero))
}

estimateG0<-function(c_pvalue, G, filename) {
  lambda<-seq(0, 1, 0.01)
  i<-1
  G0_prev<-0
  G0<-vector("integer",length(length(lambda)))
  r<-vector("integer",length(length(lambda)))
  
  while (i<length(lambda)) {
    
    l<-lambda[i]
    
    minori<-(c_pvalue<l)
    selected<-which(minori==TRUE)
    num_sel<-length(selected)

    G0[i]<-(G-num_sel)/(1-l)
    
    r[i]<-(G0[i]-G0_prev)^2
    
    G0_prev<-G0[i]
    
    i<-i+1
    
  }
  
  G0_estimate<-G0[which.min(r)]
  lambda_estimate <- lambda[which.min(r)]
  
  # plot the estimates
  png(file = paste(getPlotPath(filename, "G0 estimate"), ".png", sep = ""))
  plot(lambda, G0/G, xlab="lambda", ylab="G0", main=filename)
  points(lambda_estimate, G0_estimate/G, col= "red")
  dev.off()
  
  
  return(list(G0_estimate,lambda_estimate))
}

expected_values <- function(G,G0,alpha,selected_genes){
  expected_FP <- min(G0*alpha, selected_genes)
  expected_TP <- max(0, (selected_genes - expected_FP))
  expected_TN <- G0 - expected_FP
  expected_FN <- max(0,G-selected_genes-expected_TN)

  return (c(expected_TP,expected_FP,expected_TN,expected_FN))
}

fisher_test_matrixes <- function(matrixes){
  pval<-NULL
  for (i in (1:dim(matrixes)[1]))
  {
    a<-as.integer(matrixes[i,2])
    b<-as.integer(matrixes[i,3])
    c<-as.integer(matrixes[i,4])
    d<-as.integer(matrixes[i,5])
    m<-matrix(c(a,c,b,d), 2, 2)
    res<-fisher.test(m, alternative="greater")
    pval<-rbind(pval,res$p.value)
  } 
  return (pval)
}
  

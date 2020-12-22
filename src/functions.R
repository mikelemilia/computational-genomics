MA <- function(x, y){
  
  M <- log2(x) - log2(y)
  A <- (log2(x) + log2(y))/2
  
  return(
    list('M' = M, 'A' = A)
  )
}

produceMvA <- function(x, index, interval, folder, graph_title){
  
  extracted <- x[,index] # extract the sample i
  
  for(i in interval){
      
    png(file = paste(getPlotPath(filename = paste(index, "vs.", i, sep = " "), folder = folder), ".png", sep = ""))
      
    # select the i-th element
    selected <- x[,i]
      
    computed <- MA(extracted, selected)
    
    M <- cbind(M, computed$M) # equivalent to do computed[[1]]
    A <- cbind(A, computed$A) # equivalent to do computed[[2]]
    
    
    plot(computed$A, computed$M, xlab="A", ylab="M", main = graph_title, sub = paste("Sample", index, "vs.", i, sep = " "))
    abline(0,0)
    
    dev.off()
    
      
  }
    
}

tmm_normalization <- function(x, index, interval){
  
  extracted <- x[,index] # extract the sample i
  
  normed_samples <- matrix(extracted, nrow = genes_number, ncol = 1)
  
  for (i in interval) {
    
    selected <- x[,i]
    
    computed <- MA(extracted, selected)
    
    SF <- mean(computed$M, trim = 0.1)
    
    normed <- selected + 2^SF
    
    normed_samples <- cbind(normed_samples, normed)

  }
  
  return(
    list('samples' = normed_samples)
  )
  
}

# codice di Giulia
tmm_normalization2 <- function(DATA, refsample){

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

quantile_normalization <- function(x){
  
  cols <- ncol(x) # get number of columns
  rows <- nrow(x) # get number of rows
  
  # setup matrixes
  dataSort = matrix(0, rows, cols)
  dataIdx  = matrix(0, rows, cols)
  dataNorm = matrix(0, rows, cols)
  
  # setting the header of dataNorm
  colnames(dataNorm) <- names(x)
  
  for(i in 1:cols){ 
    data = x[,i]
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
  
  return(
    list('samples' = dataNorm)
  )
}

# codice di Giulia
quantile_normalization2 <- function(DATA){
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
  G0_estimate<-NULL
  lambda_estimate <- NULL
  
  while (i<=length(lambda)) {
    
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
  
  
  return(list(G0,lambda))
}

G0_value_estimation<-function(lambda_est, eps, res) {
  
  lambda_min<-lambda_est-eps
  lambda_max<-lambda_est+eps
  index_min<-which(res[[2]]==lambda_min)
  index_max<-which(res[[2]]==lambda_max)
  G0_min <- res[[1]][index_min]
  G0_max <- res[[1]][index_max]
  G0_est<-round(mean(G0_min, G0_max)) 
  
  return(G0_est)
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

#Silhouette Statistic
silhouette <- function(points, cluster, k){
  #caso particolare, un cluster
  if(k == 1){
    s <- -1
    return(s)
  }
  
  cat ("Analisi per k=",k)
  
  d <- as.matrix(dist(points))
  
  s <- NULL
  clusters <- vector(mode = "list", length = k)
  
  #divisione degli indici degli elementi in liste, raggruppati per il cluster a cui appartengono
  i <- 1
  while(i <= k) {
    elements <- which(cluster == i)
    clusters[[i]] <- elements 
    cat("Cluster ", i, " has length ", length(clusters[[i]]), "\n")
    cat("Cluster ", i, ": ", clusters[[i]], "\n")
    i <- i + 1
  }
  
  #scandisco ogni singolo punto
  for (i in (1:nrow(points))){
    point <- points[i]
    distances <- d[i,]
    cl <- cluster[i]
    
    if(length(distances[clusters[[cl]]]) == 1) { # ho un singleton
      
      s <- c(s,1)
      
    } else { #non ho un singleton
      
      indmin <- 0
      minb <- Inf
      #scandisco ogni possibile cluster
      for (c in (1:k)){
        cat("Analizzo punto ", i, " in cluster ",cl, " -- Confronto con cluster ", c, "\n")
        #se ho scelto il cluster a cui il punto appartiene, calcolo a
        if (c == cl){

          distances_topoint <- distances[clusters[[c]]]
          a <- (sum(distances_topoint))/(length(distances_topoint) - 1)
        }
        #se ho scelto un cluster a cui il punto non appartiene, calcolo b
        else{
          
          distances_topoint <- distances[clusters[[c]]]
          if (minb > sum(distances_topoint)){
            minb<-sum(distances_topoint)
            indmin<-c
          }
        }
      }
      
      minb <- minb/length(clusters[[indmin]])
      s <- c(s,(minb-a)/max(minb,a))

    }
  }
  return(sum(s))
}


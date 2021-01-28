MA <- function(x, y){
  
  A <- matrix(0, nrow = length(x), ncol = 0)
  M <- matrix(0, nrow = length(x), ncol = 0)
  
  M <- log2(x) - log2(y)
  A <- (log2(x) + log2(y))/2
  
  return(list('M' = M, 'A' = A))
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

producePlots <- function(x, index, interval, folder, graph_title = ""){
  
  extracted <- x[,index] # extract the sample i
  
  for(i in interval){
    png(file = paste(getPlotPath(filename = paste(index, "vs.", i, sep = " "), folder = folder), ".png", sep = ""))
    
    par(mfrow =c(2,2))
    # select the i-th element
    selected <- x[,i]
    
    computed <- MA(extracted, selected)
    M <- cbind(M, computed$M) # equivalent to do computed[[1]]
    A <- cbind(A, computed$A) # equivalent to do computed[[2]]
    
    hist(computed$A, main = graph_title)
    boxplot(computed$A, main = graph_title)
    
    hist(computed$M, main = graph_title)
    boxplot(computed$M, main = graph_title)
    
    
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


remove_duplicates <- function (data, group){
  # function that removes duplicated in the 'individual' column for 'group'
  
  # get all the duplicated individuals
  duplicates <- duplicated(group$individual)  
  duplicate_group <- group$individual[duplicates == TRUE]
  
  d <- as.numeric(group$individual)
  dim <- dim(group)
  group_nodup <- matrix(0, nrow(data), dim[1]) # ncol(group) == dim[1]
  names <- rep(0, dim[1])
  count <- 0
  
  for (i in 1:length(duplicate_group)) {
    #find samples of duplicated subjects and get the SRR code
    indexes <- which(d %in% duplicate_group[i])
    seq_sample <- group[indexes,]$seq.sample
    
    #mean of duplicated samples
    group_nodup[,i] <- apply(data[, seq_sample], 1, mean)
    names[i] <- unique(group[indexes, ]$individual)
    count <- count+1;
  }
  
  #finding subjects not duplicated
  '%notin%' <- Negate(`%in%`)
  ind<-which(d %notin% duplicate_group)
  
  for (i in 1:length(ind))  {
    #find samples of duplicated subjects, get the SRR code
    seq_sample <- group[ind[i],]$seq.sample
    #mean of duplicated samples
    group_nodup[,length(duplicate_group)+i] <- data[,seq_sample]
    names[length(duplicate_group)+i] <- unique(group[ind[i], ]$individual)
    count <- count+1;
  }
  
  group_nodup<-group_nodup[,1:count]
  
  names <- names[1:count]
  
  colnames(group_nodup) <- names
  
  return(group_nodup)
}

remove_zeros <- function(group1,group2){
  # functions to remove rows where there are only 0 in the two groups (matrixes)
  
  # sum for each column of the data
  group1_forcomparison <- apply(group1, 1, sum)
  group2_forcomparison <- apply(group2, 1, sum)
  vec_for_comparison <- as.data.frame(group1_forcomparison+group2_forcomparison)
  
  # use of the quantile to find the threshold
  q <- quantile(unlist(vec_for_comparison))
  
  # extract indexes of data to remove and remove the correspondent rows
  YN <- as.data.frame(vec_for_comparison>q[1])
  index <- which(YN==FALSE)
  group1_nozero <- group1[-index,] 
  group2_nozero <- group2[-index,]
  
  return (list('control' = group1_nozero, 'disease' = group2_nozero, 'removedindexes' = index))
}

remove_zeros_onegroup <- function(group){
  group_forcomparison <- as.data.frame(apply(group, 1, sum))
  
  # use of the quantile to find the threshold
  q <- quantile(unlist(group_forcomparison))
  
  # extract indexes of data to remove and remove the correspondent rows
  YN <- as.data.frame(group_forcomparison>q[1])
  index <- which(YN==FALSE)
  group_nozero <- group[-index,] 
  
  return (list('group' = group_nozero, 'removedindexes' = index))
}

G0values<-function(lambda, c_pvalue, G, filename) {
  G0<-NULL
  for(l in lambda) {
    minori<-(c_pvalue<l)
    selected<-which(minori==TRUE)
    num_sel<-length(selected)
    G0<-c(G0,(G-num_sel)/(1-l))
  }
  
  plot(lambda, G0/G, xlab="lambda", ylab="G0", main=filename)
  
  # plot the estimates
  png(file = paste(getPlotPath(filename, "G0 estimate"), ".png", sep = ""))
  plot(lambda, G0/G, xlab="lambda", ylab="G0", main=filename)
  dev.off()
  
  return(list('G0'=G0,'lambda'=lambda))
}

G0_value_estimation<-function(lambda_est, eps, res) {
  
  lambda_min<-lambda_est-eps
  lambda_max<-lambda_est+eps
  index_min<-which(res$lambda==lambda_min)
  index_max<-which(res$lambda==lambda_max)
  G0_min <- res$G0[index_min]
  G0_max <- res$G0[index_max]
  G0_est<-round(mean(G0_min, G0_max)) 
  
  return(G0_est)
}

expected_values <- function(G,G0,alpha,selected_genes){
  expected_FP <- min(G0*alpha, selected_genes)
  expected_TP <- max(0, (selected_genes - expected_FP))
  expected_TN <- G0 - expected_FP
  expected_FN <- max(0,G-selected_genes-expected_TN)
  
  return (list('TP'=expected_TP,'FP'=expected_FP,'TN'=expected_TN,'FN'=expected_FN))
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

FDR_fisher<-function(pval_fisher, terms){
  FDR <- 0.05
  # values in the range observed as p-value in fisher exact test
  lambda <- seq(min(out[,4]), max(out[,4]), (max(out[,4])-min(out[,4]))/nrow(out))
  FDR_values_fisher <- NULL
  
  # compute FDR for every lambda
  for (i in (1:length(lambda))) {
    less <- (pval_fisher<lambda[i])
    num_sel <- length(which(less==TRUE))
    
    #length(pval_fisher_BP) è il numero dei test fatti
    expected_FP <- min(length(pval_fisher)*lambda[i], num_sel)
    
    if (num_sel==0){
      FDR_values_fisher <- c(FDR_values_fisher,0)} 
    else {
      FDR_values_fisher <- c(FDR_values_fisher,(expected_FP/num_sel))}
  }
  
  #plot(lambda, FDR_values_fisher, main = 'FDR values for different alpha')
  
  # choose the values that are in [0.05-epsilon;0.05+espilon]
  epsilon <- 0.0001
  alpha_index <- which(FDR_values_fisher>=0.05-epsilon)
  alpha_index2 <- which(FDR_values_fisher<=0.05+epsilon)
  alpha_est <- mean(lambda[intersect(alpha_index,alpha_index2)])
  
  # selection of the correspondent indexes and ID genes in the table
  indexes_annotated <- which(pval_fisher<alpha_est) 
  terms_annotated<-terms[indexes_annotated]
  
  return (list(indexes_annotated,terms_annotated, alpha_est))
  
}

annotation_terms<-function(vals, terms){
  p<-NULL
  for(i in(1:length(terms))){
    t<-terms[i]
    p<-c(p,which(vals$GOID==t))
  }
  ann<-vals$TERM[p]
  
  return(ann)
}


extract_tumor_terms<-function(list_terms)
{sent<-list_terms$terms
ind<-grepl(" tumor ", sent)
ind<-which(ind==TRUE)

to_view<-list_terms[ind,]
return(to_view)
}


silhouette <- function(points, cluster, k){
  
  # particular case of only one cluster
  if(k == 1){return(-1)}
  
  # calculate the distances between points and variables initialization
  d <- as.matrix(dist(points))
  s <- NULL
  clusters <- vector(mode = "list", length = k)
  
  # division of the indexes in lists according to their cluster number
  i <- 1
  while(i <= k) {
    elements <- which(cluster == i)
    clusters[[i]] <- elements 
    i <- i + 1
  }
  
  # scan each point in datamatrix
  for (i in (1:nrow(points))){
    # extract the point, the distances of this point from the others in the dataset and its cluster number 
    point <- points[i]
    distances <- d[i,]
    cl <- cluster[i]
    
    if(length(distances[clusters[[cl]]]) == 1) {
      # if it is a singleton cluster, the silhouette is 1
      s <- c(s,1)
    } 
    
    else {
      # the cluster contains more points
      minb <- Inf
      
      # scan each possible cluster
      for (c in (1:k)){
        if (c == cl){
          # if we are consider the same cluster of the point, we calculate the value 'a'
          distances_topoint <- distances[clusters[[c]]]
          a <- (sum(distances_topoint))/(length(distances_topoint) - 1)
        }
        else{
          # with a different cluster, we calculate the value 'b' and we maintain always the max(b)
          distances_topoint <- distances[clusters[[c]]]
          x <- sum(distances_topoint)/length(distances_topoint)
          minb <- min(minb, x)
        }
      }
      
      # calculate the final silhouette for the point
      s <- c(s,(minb-a)/max(minb,a))
    }
  }
  
  # the output is the sum of the silhouettes for all the points in the dataset
  return(sum(s)/nrow(points))
  
}

generateDendogram <- function(x, k_opt, title){
  
  dend <-  as.dendrogram(x)
  
  fviz_dend(x,
            k = k_opt, 
            cex = 0.5, 
            k_colors = c( "#E7B800", "#FC4E07"), #"#2E9FDF", "#00AFBB",
            color_labels_by_k = T, 
            rect = T, 
            rect_fill = T, 
            show_labels = T)  
  

}

recursiveFeatureExtraction <- function(X_train, Y_train, X_test, Y_test, K = 500){
  
  model.svm <- svm(Group ~ ., data = X_train, kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
  
  prediction <- predict(model.svm, newdata = X_test, decision.values = FALSE)
  result <- confusionMatrix(prediction, factor(Y_test))
  
  # initialization
  
  R <- ncol(X_train)-1  # number of features
  accuracy <- result[["overall"]][["Accuracy"]]
  features_retained <- R
  
  cat("Number of features:", R)
  cat("\taccuracy obtained:", accuracy, "\n")
  
  per_rem <- 0.15
  eps <- 0.05
  
  bestSVM <- model.svm
  bestNames <- as.vector(colnames(X_train[,-1]))
  bestAcc <- accuracy
  
  cat("\nRemoving lot of unuseful features\n\n")
  
  while(R > K && last(accuracy) > first(accuracy)*0.9){
    
    w <- t(model.svm$coefs) %*% model.svm$SV
    names <- as.vector(colnames(X_train[,-1]))

    w <- abs(w)
    w_sort <- sort(w, decreasing = TRUE)
    names_sorted <- names[order(w, decreasing=TRUE)]
    
    m <- ceiling(per_rem * R)
    if(m <= 1) break
    w <- w_sort[1:(length(w_sort) - m)] 
    names_sorted <- names_sorted[1:(length(names_sorted) - m)]
    
    tmp_train <- X_train
    X_train <- X_train[,names_sorted]
    X_train <- cbind('Group'=factor(Y_train), X_train)
    
    R <- ncol(X_train)-1  # number of current features
    
    # train the model
    tmp_svm <- model.svm
    model.svm <- svm(Group ~ ., data = X_train, kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
    
    # make prediction on the above model
    prediction <- predict(model.svm, newdata = X_test, decision.values = FALSE)
    result <- confusionMatrix(prediction, factor(Y_test))
    
    current_accuracy <- result[["overall"]][["Accuracy"]]
    
    if (last(accuracy) > (current_accuracy + eps)){
      
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\t[x]\n")
      
      model.svm <- tmp_svm
      X_train   <- tmp_train
      per_rem <- per_rem - 0.002
      
    } else {
      
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\n")
      
      accuracy <- c(current_accuracy, accuracy)
      features_retained <- c(R, features_retained)
      
      if (current_accuracy >= bestAcc && R < length(bestNames)){
        bestSVM <- model.svm
        bestNames <- as.vector(colnames(X_train[,-1]))
        bestAcc <- current_accuracy
      }
      
    }
  }
  
  cat("\nRemoving one feature in each iteration\n\n")
  
  R <- first(features_retained)  # number of features
  
  bests <- list()
  bestsnames <- list()
  bests <- append(list(model.svm), bests)
  bestsnames <- c(list(colnames(X_train[,-1])), bestsnames)
  
  while(R >= 2){
    deltas <- NULL
    
    # find the worst feature
    for(i in 2:ncol(X_train)) {
      tmp_model.svm <- svm(Group ~ ., data = X_train[,-i], kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
      tmp_prediction <- predict(tmp_model.svm, newdata = X_test, decision.values = FALSE)
      tmp_result <- confusionMatrix(tmp_prediction, factor(Y_test))
      
      tmp_accuracy <- tmp_result[["overall"]][["Accuracy"]]
      deltas <- c(deltas, (first(accuracy) - tmp_accuracy))
    }
    
    # select the worst feature
    idx <- which(deltas == min(deltas))
    idx <- sample(idx, 1)
    
    # remove features
    X_train <- X_train[, -(idx+1)]
    
    # train the model
    model.svm <- svm(Group ~ ., data = X_train, kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
    
    # predict
    prediction <- predict(model.svm, newdata = X_test, decision.values = FALSE)
    result <- confusionMatrix(prediction, factor(Y_test))
    
    R <- ncol(X_train)-1
    current_accuracy <- result[["overall"]][["Accuracy"]]
    
    cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
    cat("\tfeatures removed:", 1, "\n")
    
    accuracy <- c(current_accuracy, accuracy)
    features_retained <- c(R, features_retained)
    
    names <- colnames(X_train[,-1])
    
    if (current_accuracy >= bestAcc && R < length(bestNames)){
      bestSVM <- model.svm
      bestNames <- as.vector(colnames(X_train[,-1]))
      bestAcc <- current_accuracy
    }
    
    bests <- append(list(model.svm), bests)
    bestsnames <- c(list(names), bestsnames)
  }
  
  plot (x = features_retained, y = accuracy)
  
  return(list("number_features" = features_retained, 
              "accuracy" = accuracy, 
              "bests" = bests, 
              "bestsnames" = bestsnames,
              "bestmodel" = list("svm" = bestSVM,
                                 "names" = bestNames,
                                 "accuracy" = bestAcc
              )
  )
  )
}


recursiveFeatureExtractionCV <- function(X, Y, K = 500, c = 5){
  
  model.svm <- svm(Group ~ ., data = X, kernel = "linear", 
                   type = 'C-classification', scale = FALSE, 
                   na.action = na.omit, cross = c)
  
  # initialization
  
  R <- ncol(X)-1  # number of features
  accuracy <- model.svm$tot.accuracy
  features_retained <- R
  
  cat("Number of features:", R, "\n")
  cat("Accuracy obtained:", accuracy, "\n")
  
  per_rem <- 0.15
  eps <- 0.05
  
  bestSVM <- model.svm
  bestNames <- as.vector(colnames(X[,-1]))
  bestAcc <- accuracy
  
  cat("Removing lot of unuseful features\n")
  
  while(R > K && last(accuracy) > first(accuracy)*0.9){
    
    w <- t(model.svm$coefs) %*% model.svm$SV
    names <- as.vector(colnames(X[,-1]))
    
    w <- abs(w)
    w_sort <- sort(w, decreasing = TRUE)
    names_sorted <- names[order(w, decreasing=TRUE)]
    
    m <- ceiling(per_rem * R)
    if(m <= 1) break
    w <- w_sort[1:(length(w_sort) - m)] 
    names_sorted <- names_sorted[1:(length(names_sorted) - m)]
    
    tmp_train <- X
    X <- X[,names_sorted]
    X <- cbind('Group' =factor(Y), X)
    
    R <- ncol(X)-1  # number of current features
    
    # train the model
    tmp_svm <- model.svm
    model.svm <- svm(Group ~ ., data = X, kernel = "linear", 
                     type = 'C-classification', scale = FALSE, 
                     na.action = na.omit, cross = c)
    
    current_accuracy <- model.svm$tot.accuracy
    
    if (last(accuracy) > (current_accuracy + eps)){
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\t[x]\n")
      
      model.svm <- tmp_svm
      X <- tmp_train
      per_rem <- per_rem - 0.002
      
    } else {
      
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\n")
      
      accuracy <- c(current_accuracy, accuracy)
      features_retained <- c(R, features_retained)
      
      if (current_accuracy >= bestAcc && R < length(bestNames)){
        bestSVM <- model.svm
        bestNames <- as.vector(colnames(X[,-1]))
        bestAcc <- current_accuracy
      }
      
    }
  }
  
  cat("\nRemoving one feature in each iteration\n")
  
  R <- first(features_retained)  # number of features
  
  bests <- list()
  bestsnames <- list()
  bests <- append(list(model.svm), bests)
  bestsnames <- c(list(colnames(X[,-1])), bestsnames)
  
  while(R >= 2){
    deltas <- NULL
    
    # find the worst feature
    for(i in 2:ncol(X)) {
      
      tmp_model.svm <- svm(Group ~ ., data = X[,-i], kernel = "linear", 
                           type = 'C-classification', scale = FALSE, 
                           na.action = na.omit, cross = c)
      
      tmp_accuracy <- model.svm$tot.accuracy
      deltas <- c(deltas, (first(accuracy) - tmp_accuracy))
      
    }
    
    # select the worst feature
    idx <- which(deltas == min(deltas))
    if (length(idx)!=1){
      w <- model.svm$finalModel@coef[[1]] %*% model.svm$finalModel@xmatrix[[1]]
      w <- w[idx] 
      idx <- which(w == min(w))
    }
    
    # remove features
    X_train <- X_train[, -(idx+1)]
    
    # train the model
    model.svm <- svm(Group ~ ., data = X_train, kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
    
    # predict
    prediction <- predict(model.svm, newdata = X_test, decision.values = FALSE)
    result <- confusionMatrix(prediction, factor(Y_test))
    
    R <- ncol(X_train)-1
    current_accuracy <- result[["overall"]][["Accuracy"]]
    
    cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
    cat("\tfeatures removed:", 1, "\n")
    
    accuracy <- c(current_accuracy, accuracy)
    features_retained <- c(R, features_retained)
    
    names <- colnames(X_train[,-1])
    
    if (current_accuracy >= bestAcc && R < length(bestNames)){
      bestSVM <- model.svm
      bestNames <- as.vector(colnames(X_train[,-1]))
      bestAcc <- current_accuracy
    }
    
    bests <- append(list(model.svm), bests)
    bestsnames <- c(list(names), bestsnames)
  }
  
  plot (x = features_retained, y = accuracy)
  
  return(list("number_features" = features_retained, 
              "accuracy" = accuracy, 
              "bests" = bests, 
              "bestsnames" = bestsnames,
              "bestmodel" = list("svm" = bestSVM,
                                 "names" = bestNames,
                                 "accuracy" = bestAcc
              )
  )
  )
}


recursiveFeatureExtractionCV <- function(X, Y, K = 500, c = 5){
  
  model.svm <- svm(Group ~ ., data = X, kernel = "linear", 
                   type = 'C-classification', scale = FALSE, 
                   na.action = na.omit, cross = c)
  
  # initialization
  
  R <- ncol(X)-1  # number of features
  accuracy <- model.svm$tot.accuracy
  features_retained <- R
  
  cat("Number of features:", R, "\n")
  cat("Accuracy obtained:", accuracy, "\n")
  
  per_rem <- 0.15
  eps <- 0.05
  
  bestSVM <- model.svm
  bestNames <- as.vector(colnames(X[,-1]))
  bestAcc <- accuracy
  
  cat("Removing lot of unuseful features\n")
  
  while(R > K && last(accuracy) > first(accuracy)*0.9){
    
    w <- t(model.svm$coefs) %*% model.svm$SV
    names <- as.vector(colnames(X[,-1]))
    
    w <- abs(w)
    w_sort <- sort(w, decreasing = TRUE)
    names_sorted <- names[order(w, decreasing=TRUE)]
    
    m <- ceiling(per_rem * R)
    if(m <= 1) break
    w <- w_sort[1:(length(w_sort) - m)] 
    names_sorted <- names_sorted[1:(length(names_sorted) - m)]
    
    tmp_train <- X
    X <- X[,names_sorted]
    X <- cbind('Group' =factor(Y), X)
    
    R <- ncol(X)-1  # number of current features
    
    # train the model
    tmp_svm <- model.svm
    model.svm <- svm(Group ~ ., data = X, kernel = "linear", 
                     type = 'C-classification', scale = FALSE, 
                     na.action = na.omit, cross = c)
    
    current_accuracy <- model.svm$tot.accuracy
    
    if (last(accuracy) > (current_accuracy + eps)){
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\t[x]\n")
      
      model.svm <- tmp_svm
      X <- tmp_train
      per_rem <- per_rem - 0.002
      
    } else {
      
      cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
      cat("\tfeatures removed:", m, "\n")
      
      accuracy <- c(current_accuracy, accuracy)
      features_retained <- c(R, features_retained)
      
      if (current_accuracy >= bestAcc && R < length(bestNames)){
        bestSVM <- model.svm
        bestNames <- as.vector(colnames(X[,-1]))
        bestAcc <- current_accuracy
      }
      
    }
  }
  
  cat("\nRemoving one feature in each iteration\n")
  
  R <- first(features_retained)  # number of features
  
  bests <- list()
  bestsnames <- list()
  bests <- append(list(model.svm), bests)
  bestsnames <- c(list(colnames(X[,-1])), bestsnames)
  
  while(R >= 2){
    deltas <- NULL
    
    # find the worst feature
    for(i in 2:ncol(X)) {
      
      tmp_model.svm <- svm(Group ~ ., data = X[,-i], kernel = "linear", 
                           type = 'C-classification', scale = FALSE, 
                           na.action = na.omit, cross = c)
      
      tmp_accuracy <- model.svm$tot.accuracy
      deltas <- c(deltas, (first(accuracy) - tmp_accuracy))
      
    }
    
    # select the worst feature
    idx <- which(deltas == min(deltas))
    if (length(idx)!=1){
      w <- t(model.svm$coefs) %*% model.svm$SV
      w <- w[idx] 
      idx <- which(w == min(w))
    }
    
    # remove features
    X <- X[, -(idx+1)]
    
    # train the model
    model.svm <- svm(Group ~ ., data = X, kernel = "linear", 
                     type = 'C-classification', scale = FALSE, 
                     na.action = na.omit, cross = c)
    
    R <- ncol(X)-1
    current_accuracy <- model.svm$tot.accuracy
    
    cat("Number of features:", R, "\taccuracy obtained:", current_accuracy)
    cat("\tfeatures removed:", 1, "\n")
    
    accuracy <- c(current_accuracy, accuracy)
    features_retained <- c(R, features_retained)
    
    names <- colnames(X[,-1])
    
    if (current_accuracy >= bestAcc && R < length(bestNames)){
      bestSVM <- model.svm
      bestNames <- as.vector(colnames(X[,-1]))
      bestAcc <- current_accuracy
    }
    
    bests <- append(list(model.svm), bests)
    bestsnames <- c(list(names), bestsnames)
  }
  
  plot (x = features_retained, y = accuracy, xlab = "Features retained", ylab = "Accuracy")
  
  return(list("number_features" = features_retained, 
              "accuracy" = accuracy, 
              "bests" = bests, 
              "bestsnames" = bestsnames,
              "bestmodel" = list("svm" = bestSVM,
                                 "names" = bestNames,
                                 "accuracy" = bestAcc
              )
  )
  )
}

recursiveFeatureExtractionCVCARET <- function(X, Y, K = 500, c = 5){
  
  cat("Repeated ", c, "-fold cross validation applied\n\n", sep = "")
  
  model.svm <- train(Group~., data= X,
                     method = "svmLinear",
                     scale = FALSE,  
                     metric = "Accuracy",
                     maximize = TRUE, 
                     trControl = trainControl(method = "repeatedcv", number = c, repeats = 10),
                     na.action=na.omit)
  
  # initialization
  
  R <- ncol(X)-1  # number of features
  accuracy <- model.svm$results$Accuracy
  features_retained <- R
  
  cat("Number of features:", R, "\taccuracy obtained:", accuracy, "\n\n")
  
  per_rem <- 0.15
  eps <- 0.05
  
  bestSVM <- model.svm
  bestNames <- as.vector(colnames(X[,-1]))
  bestAcc <- accuracy
  
  cat("Removing lots of unuseful features\n")
  
  while(R > K && last(accuracy) > first(accuracy) * 0.9){
    
    w <- model.svm$finalModel@coef[[1]] %*% model.svm$finalModel@xmatrix[[1]]
    names <- as.vector(colnames(X[,-1]))
    
    w <- abs(w)
    w_sort <- sort(w, decreasing = TRUE)
    names_sorted <- names[order(w, decreasing = TRUE)]
    
    m <- ceiling(per_rem * R)
    if(m <= 1) break
    w <- w_sort[1:(length(w_sort) - m)] 
    names_sorted <- names_sorted[1:(length(names_sorted) - m)]
    
    tmp_train <- X
    X <- X[, names_sorted]
    X <- cbind('Group' = factor(Y), X)
    
    R <- ncol(X)-1  # number of current features
    
    # train the model
    tmp_svm <- model.svm
    model.svm <- train(Group~., data= X,
                       method = "svmLinear",
                       scale = FALSE,  
                       metric = "Accuracy",
                       maximize = TRUE, 
                       trControl = trainControl(method = "repeatedcv", number = c, repeats = 10),
                       na.action=na.omit)
    
    current_accuracy <- model.svm$results$Accuracy
    
    if (last(accuracy) > (current_accuracy + eps)){
      
      cat("Number of features:", R, "\tfeatures removed:", m)
      cat("\taccuracy obtained:", current_accuracy, "\t[x]\n")
      
      model.svm <- tmp_svm
      X <- tmp_train
      per_rem <- per_rem - 0.002
      
    } else {
      
      cat("Number of features:", R, "\tfeatures removed:", m)
      cat("\taccuracy obtained:", current_accuracy, "\n")
      
      accuracy <- c(current_accuracy, accuracy)
      features_retained <- c(R, features_retained)
      
      if (current_accuracy >= bestAcc && R < length(bestNames)){
        bestSVM <- model.svm
        bestNames <- as.vector(colnames(X[,-1]))
        bestAcc <- current_accuracy
      }
      
    }
  }
  
  cat("\nRemoving one feature in each iteration\n")
  
  R <- first(features_retained)  # number of features
  
  bests <- list()
  bestsnames <- list()
  bests <- append(list(model.svm), bests)
  bestsnames <- c(list(colnames(X[,-1])), bestsnames)
  
  while(R >= 2){
    deltas <- NULL
    
    # find the worst feature
    for(i in 2:ncol(X)) {
      cat(i, " - ")
      tmp_model.svm <- train(Group~., data= X,
                             method = "svmLinear",
                             scale = FALSE,  
                             metric = "Accuracy",
                             maximize = TRUE, 
                             trControl = trainControl(method = "cv", number = c),
                             na.action=na.omit)
      
      tmp_accuracy <- model.svm$results$Accuracy
      deltas <- c(deltas, (first(accuracy) - tmp_accuracy))
      
    }
    
    # select the worst feature
    idx <- which(deltas == min(deltas))
    if (length(idx)!=1){
      w <- model.svm$finalModel@coef[[1]] %*% model.svm$finalModel@xmatrix[[1]]
      w <- w[idx] 
      idx <- which(w == min(w))
    }
    
    # remove features
    X <- X[, -(idx+1)]
    
    # train the model
    model.svm <- train(Group~., data= X,
                       method = "svmLinear",
                       scale = FALSE,  
                       metric = "Accuracy",
                       maximize = TRUE, 
                       trControl = trainControl(method = "repeatedcv", number = c, repeats = 10),
                       na.action=na.omit)
    
    R <- ncol(X)-1
    current_accuracy <- model.svm$results$Accuracy
    
    cat("Number of features:", R, "\tfeatures removed:", 1)
    cat("\taccuracy obtained:", current_accuracy, "\n")
    
    accuracy <- c(current_accuracy, accuracy)
    features_retained <- c(R, features_retained)
    
    names <- colnames(X[,-1])
    
    if (current_accuracy >= bestAcc && R < length(bestNames)){
      bestSVM <- model.svm
      bestNames <- as.vector(colnames(X[,-1]))
      bestAcc <- current_accuracy
    }
    
    bests <- append(list(model.svm), bests)
    bestsnames <- c(list(names), bestsnames)
  }
  
  plot (x = features_retained, y = accuracy, xlab = "Features retained", ylab = "Accuracy")
  
  return(list("number_features" = features_retained, 
              "accuracy" = accuracy, 
              "bests" = bests, 
              "bestsnames" = bestsnames,
              "bestmodel" = list("svm" = bestSVM,
                                 "names" = bestNames,
                                 "accuracy" = bestAcc
              )
  )
  )
}

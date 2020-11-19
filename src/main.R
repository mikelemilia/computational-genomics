# 1 - Loading data

DATA<-read.table(getRawPath("E-GEOD-76987-raw-counts.tsv"), sep = "\t", row.names = 1, header = TRUE)
LABELS <- read.delim(getRawPath("labels.txt"), sep = "\t", header = TRUE)


# 2 - Calculate the sequencing depth of each sample 

depth <- apply(DATA[-1], 2, sum)

# 3 - Produce the MvA plots of each sample vs. sample 1

indsample <- 2;
sample <- DATA[,indsample] + 1 # extract the sample 1

cols <- ncol(DATA[-1]) # get number of columns
rows <- length(sample) # get number of rows

A <- matrix(0, nrow = rows, ncol = 0)
M <- matrix(0, nrow = rows, ncol = 0)

for(i in 3:cols){
  png(file = paste(getPlotPath(i-2, "MvA"), ".png", sep = ""))
  
  # select the i-th element
  selected <- DATA[,i] + 1
  
  result <- ma(sample, selected)
  
  M <- cbind(M, result[[1]])
  A <- cbind(A, result[[2]])
  
  plot(result[[2]], result[[1]], ylab="M", xlab="A", main=paste("MvA Plot", i-2, sep =" "))
  abline(0,0)
  
  dev.off()
}

# 4 - Produce the normalization of data with TMM and print MvA plots

normed <- tmm_normalization(DATA[-1],indsample)
dataNorm <- normed[[3]]

for(i in 3:cols) {
  png(file = paste(getPlotPath(i-2, "MvA - TMM Normalized"), ".png", sep = ""))
  
  plot(normed[[2]],normed[[1]],ylab="M",xlab="A",main=paste("MvA Plot", i-2, sep =" "))
  abline(0,0)
  
  dev.off()
}

# 5 - Produce the normalization of data with Quantiles, calculate new M and A and print MvA plots 

dataNorm <- quantile_normalization(DATA[-1])

for(i in 3:cols){
  png(file = paste(getPlotPath(i-2, "MvA - Quantile Normalized"), ".png", sep = ""))
  
  # select the i-th element
  selected <- dataNorm[,i] + 1
  
  result <- ma(sample, selected)
  
  M <- cbind(M, result[[1]])
  A <- cbind(A, result[[2]])
  
  plot(result[[2]], result[[1]], ylab="M", xlab="A", main=paste("MvA Plot", i-2, sep =" "))
  abline(0,0)
  
  dev.off()
}

# 6 - Remove the duplicated individuals

normal <- LABELS[LABELS$sample_type == c("normal"),] # get all the normal samples
unimuc <- LABELS[LABELS$sample_type == c("uninvolved mucosa"),] # get all the uninvolved mucosa samples
control <- rbind(normal, unimuc) # concatenate them 

# TODO sort 

disease <- LABELS[LABELS$sample_type == c("colon sessile serrated adenoma/polyp"),] # get all the disease samples

#TAKING CARE OF DUPLICATES FOR CONTROL
duplicates <- duplicated(control$individual) # get the position of all the duplicated individual 
#get only the duplicated
duplicate_control <- control$individual[duplicates == TRUE]
#find indexes in control of duplicated subjects
d <- as.numeric(control$individual)
#give names to dataNorm
colnames(dataNorm)<-names((DATA[-1])[1,])

dim <- dim(control)
control_nodup <- matrix(0,nrow(DATA),dim[1])
#initialize a count because matrix at the end WON'T have same columns as dim[1]=41, the control SRR
count<-0
for (i in 1:length(duplicate_control))
{
  indexes <- which(d %in% duplicate_control[i])
  #find samples of duplicated subjects, get the SRR code
  seq_sample <- control[indexes,]$seq.sample
  #mean of duplicated samples
  control_nodup[,i] <- apply(dataNorm[, seq_sample], 1, mean)
  count <- count+1;
}

#finding subjects NOT DUPLICATED
'%notin%' <- Negate(`%in%`)
ind<-which(d %notin% duplicate_control)

for (i in 1:length(ind))
{
  
  #find samples of duplicated subjects, get the SRR code
  seq_sample<-control[ind[i],]$seq.sample
  #mean of duplicated samples
  control_nodup[,length(duplicate_control)+i]<-dataNorm[,seq_sample]
  count<-count+1;
}
control_nodup<-control_nodup[,1:count]

#TAKING CARE OF DUPLICATES FOR DISEASE

duplicates <- duplicated(disease$individual) # get the position of all the duplicated individual 
#get only the duplicated
duplicate_disease <- disease$individual[duplicates == TRUE]
#find indexes in disease of duplicated subjects
d<-as.numeric(disease$individual)

dim<-dim(disease)
disease_nodup<-matrix(0,nrow(DATA),dim[1])
#initialize a count beacuse matrix at the end WON'T have same columns as dim[1]=41, the control SRR
count<-0
for (i in 1:length(duplicate_disease))
{
  indexes<-which(d %in% duplicate_disease[i])
  #find samples of duplicated subjects, get the SRR code
  seq_sample<-disease[indexes,]$seq.sample
  #mean of duplicated samples
  disease_nodup[,i]<-apply(dataNorm[, seq_sample], 1, mean)
  count<-count+1;
}

#finding subjects NOT DUPLICATED
'%notin%' <- Negate(`%in%`)
ind<-which(d %notin% duplicate_disease)
for (i in 1:length(ind))
{
  #find samples of duplicated subjects, get the SRR code
  seq_sample<-disease[ind[i],]$seq.sample
  #mean of duplicated samples
  disease_nodup[,length(duplicate_disease)+i]<-dataNorm[,seq_sample]
  count<-count+1
  
}
disease_nodup<-disease_nodup[,1:count]


# --------------- TAKING CARE OF ZEROS in controls and diseased ---------------

control_forcomparison<-apply(control_nodup, 1, sum)
disease_forcomparison<-apply(disease_nodup, 1, sum)
vec_for_comparison<-as.data.frame(control_forcomparison+disease_forcomparison)
#usare i quantili PER TROVARE SOGLIA 
q<-quantile(unlist(vec_for_comparison))
#median<-median(unlist(vec_for_comparison)) USANDO MEDIANA? TROPPO POCHI NE RESTANO
YN<-as.data.frame(vec_for_comparison>q[1])
index<-which(YN==FALSE) #indici di quelli da togliere
control_nodup_nozero<-(control_nodup[-index,])
disease_nodup_nozero<-disease_nodup[-index, ]
# ------------- tests  ------------------------
Nc<-nrow(control_nodup_nozero)
c_ttest_pvalue <- NULL
c_wilcoxon_pvalue <- NULL
alpha<-0.05
selected_ttest<-0
selected_wilcox<-0
for(i in (1:Nc)){ 
  c_ttest_pvalue <- c(c_ttest_pvalue,t.test(control_nodup_nozero[i,], disease_nodup_nozero[i,], var.equal = FALSE)[[3]])
  
  c_wilcoxon_pvalue <- c(c_wilcoxon_pvalue,wilcox.test(control_nodup_nozero[i,],disease_nodup_nozero[i,], exact=FALSE)[[3]])
  
  #INSERIRE EDGER
}
minori<-(c_ttest_pvalue<0.05 )
selected_ttest<-which(minori==TRUE)
num_sel_ttest<-length(selected_ttest)
minori<-(c_wilcoxon_pvalue<0.05 )
selected_wilcox<-which(minori==TRUE)
num_sel_wilcox<-length(selected_wilcox)

# --------- FP and FN with G0=G --------------------------------
G<-nrow(control_nodup_nozero)
G0<-G
alpha<-0.05
expected_FP_ttest<-min(G0*alpha, num_sel_ttest)
expected_TP_ttest<-max(0, (num_sel_ttest - expected_FP_ttest))
expected_TN_ttest<-G0 - expected_FP_ttest
expected_FN_ttest<- max(0,G -num_sel_ttest- expected_TN_ttest)




expected_FP_wilcox<-min(G0*alpha, num_sel_wilcox)
expected_TP_wilcox<-max(0, (num_sel_wilcox - expected_FP_wilcox))
expected_FN_wilcox<-max(0,G-num_sel_wilcox - expected_FN_wilcox)
expected_TN_wilcox<-G0 - expected_FP_wilcox




# ---------------------estimate G0 and re-estimate FP and FN-----------------------------

G0_est <- estimateG0(c_ttest_pvalue)











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
  png(file = paste(getPlotPath(i-2, "MvA - TMM Normalization"), ".png", sep = ""))
  
  plot(normed[[2]][,i], normed[[1]][,i], ylab="M",xlab="A",main=paste("MvA - TMM Plot", i-2, sep =" "))
  abline(0,0)
  
  dev.off()
}

# 5 - Produce the normalization of data with Quantiles, calculate new M and A and print MvA plots 

dataNorm <- quantile_normalization(DATA[-1])

for(i in 3:cols){
  png(file = paste(getPlotPath(i-2, "MvA - Quantile Normalization"), ".png", sep = ""))
  
  # select the i-th element
  selected <- dataNorm[,i] + 1
  
  result <- ma(sample, selected)
  
  M <- cbind(M, result[[1]])
  A <- cbind(A, result[[2]])
  
  plot(result[[2]], result[[1]], ylab="M", xlab="A", main=paste("MvA - Quantile Plot", i-2, sep =" "))
  abline(0,0)
  
  dev.off()
}

#give names to dataNorm
colnames(dataNorm)<-names((DATA[-1])[1,])

# 6 - Remove the duplicated individuals

normal <- LABELS[LABELS$sample_type == c("normal"),] # get all the normal samples
unimuc <- LABELS[LABELS$sample_type == c("uninvolved mucosa"),] # get all the uninvolved mucosa samples
control <- rbind(normal, unimuc) # concatenate them 
control<-control[order(as.numeric(control$individual)),] #sort in function of 'individual' value

disease <- LABELS[LABELS$sample_type == c("colon sessile serrated adenoma/polyp"),] # get all the disease samples

control_nodup<-remove_duplicates(dataNorm,control) #Remove duplicates in control
disease_nodup<-remove_duplicates(dataNorm,disease) #Remove duplicates in disease

# 7 - Taking care of zeros in controls and diseased

groups_nozero <- remove_zeros(control_nodup,disease_nodup)
control_nodup_nozero<-groups_nozero[[1]]
disease_nodup_nozero<-groups_nozero[[2]]

# 8 - t-test and Wilcoxon test

Nc<-nrow(control_nodup_nozero)
c_ttest_pvalue <- NULL
c_wilcoxon_pvalue <- NULL

selected_ttest<-0
selected_wilcox<-0
for(i in (1:Nc)){ 
  c_ttest_pvalue <- c(c_ttest_pvalue,t.test(control_nodup_nozero[i,], disease_nodup_nozero[i,], var.equal = FALSE)[[3]])
  c_wilcoxon_pvalue <- c(c_wilcoxon_pvalue,wilcox.test(control_nodup_nozero[i,],disease_nodup_nozero[i,], exact=FALSE)[[3]])
}

# 9 - Pre-processing and EdgeR

#Rebuilt control matrix with initial data (not normalized)
control_nodup<-remove_duplicates(DATA,control)

count<-ncol(control_nodup)
nomi<-rep("", count)
for (i in (1:count))
{
  x<-paste("control",as.character(i),sep='_')
  nomi[i]<-x
}
colnames(control_nodup)<-nomi

#Rebuilt disease matrix with initial data (not normalized)

disease_nodup<-remove_duplicates(DATA,disease)

count<-ncol(disease_nodup)
nomi<-rep("", count)
for (i in (1:count))
{
  x<-paste("disease",as.character(i),sep='_')
  nomi[i]<-x
}
colnames(disease_nodup)<-nomi

#------------------------- edgeR-----------------------------------

#unisco le due tabelle
mat <- cbind(disease_nodup, control_nodup)
library(edgeR)

# creo i gruppi che poi mi serviranno per definire "group"; di base mi serve solo che mi separino
#ciò che c'è in tabella in control e disease
#sono due vettori
gruppo_controllo <- rep("control",dim(control_nodup)[2])
gruppo_malato <- rep("disease",dim(disease_nodup)[2])

#li unisco, mi serve vettore unico per usare factor
gruppo <- cbind(t(as.data.frame(gruppo_malato)),t(as.data.frame(gruppo_controllo)))

#uso factor(), ottengo un oggetto diviso in due livelli (contorl e disease), come ci serve
group <- factor(gruppo)

#matrice di design
design <- model.matrix(~0+group) 
rownames(design) <- colnames(mat)  
print(design)

# fit values of phi (we need this step to fit our GLM model)
y <- DGEList(counts=mat, remove.zeros = TRUE)    # y is an object of type DGE
y <- calcNormFactors(y)   # This calculates the SF using the TMM normalization !!!
SF<-y$samples

y <- estimateGLMCommonDisp(y,design, verbose=TRUE) #phi common to the entire dataset
y <- estimateGLMTrendedDisp(y,design) #phi depends on mu
y <- estimateGLMTagwiseDisp(y,design) #phi is gene specific
fit <- glmFit(y,design) #finally the model fit (that accounts for raw NB data and scaling factors and seq. depth) 
summary(fit)

#il test
Confronti<-makeContrasts(Treatment=groupdisease-groupcontrol,levels=design)
RES<-glmLRT(fit,contrast=Confronti[,"Treatment"])

#alcuni output

# The first column of RES reports the log_Fold_Change, i.e.: 
# log2(Normalized_data_average_groupProvadisease / Normalized_data_average_groupProvacontrol)
RES$table[1:5,]

out <- topTags(RES, n = "Inf")$table
out[1:5,]

# 10 - Selection of the genes for alpha = 0.05

alpha<-0.05

minori<-(c_ttest_pvalue<alpha)
selected_ttest<-which(minori==TRUE)
num_sel_ttest<-length(selected_ttest)

minori<-(c_wilcoxon_pvalue<alpha)
selected_wilcox<-which(minori==TRUE)
num_sel_wilcox<-length(selected_wilcox)

#selected usando edgeR 
indSELedgeR<-length(which(out$PValue<alpha)) #i selected

# 11 - E[FP] and E[FN] for t-test and Wilcoxon test with G0 = G

G<-nrow(control_nodup_nozero)
G0<-G

#function that returns a vector with, in order, TP FP TN FN
expected_ttest <- expected_values(G, G0, alpha, num_sel_ttest)
expected_wilcoxontest <- expected_values(G, G0, alpha, num_sel_wilcox)
expected_edger <- expected_values(G, G0, alpha, indSELedgeR)

# 12 - Estimate G0 and re-estimate FP and FN

res <- estimateG0(c_ttest_pvalue, G0, "T test")
G0_est_ttest <- res[[1]]
lambda_est <- res[[2]]

expected_ttest_est <- expected_values(G, G0_est_ttest, alpha, num_sel_ttest)

res <- estimateG0(c_wilcoxon_pvalue, G0, "Wilcoxon test")
G0_est_wilcoxontest <- res[[1]]
lambda_est <- res[[2]]

expected_wilcoxontest_est <- expected_values(G, G0_est_wilcoxontest, alpha, num_sel_wilcox)

res <- estimateG0(out[,4], G0, "edgeR test")
G0_est_edger <- res[[1]]
lambda_est <- res[[2]]

expected_edgeR_est <- expected_values(G, G0_est_edger, alpha, indSELedgeR)

# 13 - Select the final list of DE genes and FDR = 5%

FDR <- 0.05
#values observed as p-value in edgeR (test scelto)
lambda<-seq(min(out[,4]), max(out[,4]), (max(out[,4])-min(out[,4]))/nrow(out))
FDR_values<-NULL

for (i in (1:length(lambda))) {
  #compute FDR for every lambda
  minori<-(out[,4]<lambda[i])
  num_sel<-length(which(minori==TRUE))
  
  expected_val <- expected_values(G, G0_est_edger, lambda[i], num_sel)
  
  if (num_sel==0){
    FDR_values<-c(FDR_values,0)} 
  else {
    FDR_values<-c(FDR_values,(expected_val[2]/num_sel))}
}

plot(lambda,FDR_values)

#choose the values that are in [0.05-epsilon;0.05+espilon]
#non possiamo usare = 0.05 perchè nessuno lo ritorna esattamente
epsilon <- 0.0001
alpha_index <- which(FDR_values>=0.05-epsilon)
alpha_index2 <- which(FDR_values<=0.05+epsilon)
alpha_est <- mean(lambda[intersect(alpha_index,alpha_index2)])

indexes<-which(out$PValue<alpha_est) #i selected
index_genes_selected<-sort(as.numeric(rownames(out[indexes,])))
names_genes_selected<-DATA[index_genes_selected,1]
number_genes_selected<-length(names_genes_selected)
names_genes_notselected <- DATA[-index_genes_selected,1]

# 14 - estrazione di tutti i GOterm associati ai geni selezionati e creazione delle tabelle associate  

library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)

alldata <- select(org.Hs.eg.db, names_genes_selected, columns = c("SYMBOL","ENTREZID","GOALL"), keytype="SYMBOL")
terms <- unique(alldata[,3])
terms <- terms[!is.na(terms)]
matrixes <- NULL

#TODO: eliminare da tutte le liste i geni non annotati

for (i in (1:length(terms))){
  GOterm <- terms[i]
  GOterm_indexes <- which(alldata$GOALL==GOterm)
  a <- length(intersect(alldata[GOterm_indexes,1],names_genes_selected))
  b <- number_genes_selected - a
  c <- length(GOterm_indexes) - a
  d <- length(names_genes_notselected)- c
  type<-alldata[(which(alldata[,3]==GOterm))[1],5]
  #matrice che ha nelle righe i GOterms associati e nelle colonne il tipo di GOterm e i valori di a,b,c,d per il fisher test 
  matrixes <- rbind(matrixes,c(type,a,b,c,d))
}

colnames(matrixes)<-c("type","a","b","c","d")
rownames(matrixes)<-terms
matrixes<-matrixes[order(rownames(matrixes)),]

# 15 - divisione delle tabelle per type e computazione del fisher test

matrixesCC <- matrixes[which(matrixes[,1]=="CC"),]
colnames(matrixesCC)<-c("type","a","b","c","d")
matrixesBP <- matrixes[which(matrixes[,1]=="BP"),]
colnames(matrixesBP)<-c("type","a","b","c","d")
matrixesMF <- matrixes[which(matrixes[,1]=="MF"),]
colnames(matrixesMF)<-c("type","a","b","c","d")

pval_fisherCC<-fisher_test_matrixes(matrixesCC)
pval_fisherBP<-fisher_test_matrixes(matrixesBP)
pval_fisherMF<-fisher_test_matrixes(matrixesMF)

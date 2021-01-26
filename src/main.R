# Loading all the useful libraries

library(FactoMineR)
library(edgeR)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(EDASeq)
library(cluster)
library(caret)
library(e1071)

# Load all the function used

source(paste(getwd(), 'utilities.R', sep = ""))
source(paste(getwd(), 'functions.R', sep = ""))

# Part 1

DATA   <- read.table("E-GEOD-76987-raw-counts.tsv", sep = "\t", row.names = 1, header = TRUE)
LABELS <- read.delim("labels.txt", sep = "\t", header = TRUE)

genes <- DATA[,1]                
samples <- DATA[,2:ncol(DATA)]
genes_number <- length(genes)         
samples_number <- length(samples)  

# Part 2

# get the 'normal' and 'uninvolved mucosa' samples (samples of the control group) 
normal <- LABELS[LABELS$sample_type == c("normal"),]            
unimuc <- LABELS[LABELS$sample_type == c("uninvolved mucosa"),] 

# concatenate and sort samples
control <- rbind(normal, unimuc)                               
control <- control[order(as.numeric(control$individual)),]    

# get the 'colon sessile serrated adenoma/polyp' samples (samples of the disease group)
disease <- LABELS[LABELS$sample_type == c("colon sessile serrated adenoma/polyp"),]

# remove duplicates in the groups
control_nodup <- remove_duplicates(quantile_normed$samples, control)
disease_nodup <- remove_duplicates(quantile_normed$samples, disease) 

groups_nozero <- remove_zeros(control_nodup, disease_nodup)
control_nodup_nozero<-groups_nozero$control
disease_nodup_nozero<-groups_nozero$disease

Nc <- nrow(control_nodup_nozero)
c_ttest_pvalue <- NULL
c_wilcoxon_pvalue <- NULL

# use of the t.test and wilcoxon.test functions for each gene
for(i in (1:Nc)){ 
  c_ttest_pvalue <- c(c_ttest_pvalue,t.test(control_nodup_nozero[i,], disease_nodup_nozero[i,], var.equal = FALSE)$p.value)
  c_wilcoxon_pvalue <- c(c_wilcoxon_pvalue,wilcox.test(control_nodup_nozero[i,],disease_nodup_nozero[i,], exact=FALSE)$p.value)
}

# rebuilt control matrix and label with the proper name
control_nodup <- remove_duplicates(DATA, control)
control_nodup <- renameColumns(control_nodup, "control")

# rebuilt disease matrix and label with the proper name 
disease_nodup <- remove_duplicates(DATA, disease)
disease_nodup <- renameColumns(disease_nodup, "disease")

# merge the two tables of disease and control 
mat <- cbind(disease_nodup, control_nodup)

# create groups and merge them
control_group <- rep("control",dim(control_nodup)[2])
disease_group <- rep("disease",dim(disease_nodup)[2])
group_all <- cbind(t(as.data.frame(disease_group)),t(as.data.frame(control_group)))

# with factor() we get an object divided in two levels (control and disease)
group <- factor(group_all)

# design matrix
design <- model.matrix(~0+group) 
rownames(design) <- colnames(mat)

# fit values of phi (step to fit our GLM model)
y <- DGEList(counts=mat, remove.zeros = TRUE)    
y <- calcNormFactors(y)   # scaling factors with TMM approach
SF <- y$samples

y <- estimateGLMCommonDisp(y,design, verbose=TRUE) # phi common to the entire dataset
y <- estimateGLMTrendedDisp(y,design) # phi depends on mu
y <- estimateGLMTagwiseDisp(y,design) # phi is gene specific
fit <- glmFit(y,design) # the model fit 

# the test
Confr <- makeContrasts(Treatment=groupdisease-groupcontrol,levels=design)
RES <- glmLRT(fit,contrast=Confr[,"Treatment"])

# some outputs to give an idea of the test's results
RES$table[1:5,]

# final values
out <- topTags(RES, n = "Inf")$table

alpha<-0.05

# selected values for t-test
lower <- (c_ttest_pvalue<alpha)
selected_ttest <- which(lower==TRUE)
num_sel_ttest <- length(selected_ttest)

# selected values for Wilcoxon test
lower <- (c_wilcoxon_pvalue<alpha)
selected_wilcox <- which(lower==TRUE)
num_sel_wilcox <- length(selected_wilcox)

#selected values for edgeR
num_sel_edgeR <- length(which(out$PValue<alpha)) #i selected

G <- nrow(control_nodup_nozero)
G0 <- G

#function that returns a vector with, in order, TP FP TN FN
expected_ttest <- expected_values(G, G0, alpha, num_sel_ttest)
expected_wilcoxontest <- expected_values(G, G0, alpha, num_sel_wilcox)
expected_edger <- expected_values(G, G0, alpha, num_sel_edgeR)

lambda<-seq(0, 0.99, 0.01)

# t-test analysis
res <- G0values(lambda,c_ttest_pvalue, G0, "T test")
lambda_est_ttest <- 0.8
eps <- 0.03
G0_est_ttest <- G0_value_estimation(lambda_est_ttest, eps, res)

expected_ttest_est <- expected_values(G, G0_est_ttest, alpha, num_sel_ttest)

# Wilcoxon analysis
res <- G0values(lambda,c_wilcoxon_pvalue, G0, "Wilcoxon test")
lambda_est_wilcoxon <- 0.8
eps <- 0.03
G0_est_wilcoxontest <- G0_value_estimation(lambda_est_wilcoxon, eps, res)

expected_wilcoxontest_est <- expected_values(G, G0_est_wilcoxontest, alpha, num_sel_wilcox)

# edgeR analysis
res <- G0values(lambda,out[,4], G0, "edgeR test")
lambda_est_edger <- 0.65
eps <- 0.03
G0_est_edger <- G0_value_estimation(lambda_est_edger, eps, res)

expected_edgeR_est <- expected_values(G, G0_est_edger, alpha, num_sel_edgeR)

FDR <- 0.05
# values in the range observed as p-value in edgeR
lambda <- seq(min(out[,4]), max(out[,4]), (max(out[,4])-min(out[,4]))/nrow(out))
FDR_values <- NULL

# compute FDR for every lambda
for (i in (1:length(lambda))) {
  less <- (out[,4]<lambda[i])
  num_sel <- length(which(less==TRUE))
  
  expected_val <- expected_values(G, G0_est_edger, lambda[i], num_sel)
  
  if (num_sel==0) FDR_values <- c(FDR_values, 0)
  else FDR_values <- c(FDR_values, (expected_val$FP/num_sel))
  
}

# choose the values that are in [0.05-epsilon;0.05 + espilon]
epsilon <- 0.0001
alpha_idx_lower <- which(FDR_values >= FDR - epsilon)
alpha_idx_upper <- which(FDR_values <= FDR + epsilon)
alpha_est <- mean(lambda[intersect(alpha_idx_lower, alpha_idx_upper)])

indexes <- which(out$PValue<alpha_est) 
index_genes_selected <- sort(as.numeric(rownames(out[indexes,])))

ID_genes_selected <- rownames(DATA[index_genes_selected,])
number_genes_selected <- length(ID_genes_selected)
ID_genes_notselected <- setdiff(rownames(DATA),ID_genes_selected)
number_genes_notselected <- length(ID_genes_notselected)

# Part 3

# extraction of the associated terms in function of the ENSEMBL ID 
alldata <- select(org.Hs.eg.db, ID_genes_selected, columns = c("SYMBOL","ENTREZID", "ENSEMBL","GOALL"), keytype="ENSEMBL")

# remove the genes for which we have the NA term associated 
GOALL_NA <- which(is.na(alldata$GOALL))
ID_goall_na <- alldata$ENSEMBL[GOALL_NA]

ID_genes_selected_notna <- setdiff(ID_genes_selected,ID_goall_na)
number_genes_selected_notna <- length(ID_genes_selected_notna)
ID_genes_notselected_notna <- setdiff(ID_genes_notselected,ID_goall_na)
number_genes_notselected_notna <- length(ID_genes_notselected_notna)

# remove the duplicates in the extracted terms 
terms <- unique(alldata[,4])
terms <- terms[!is.na(terms)]

#creation of the matrix for the Fisher test, one row for each GOterm
matrixes <- NULL
for (i in (1:length(terms))){
  GOterm <- terms[i]
  GOterm_indexes <- which(alldata$GOALL==GOterm)
  a <- length(intersect(alldata[GOterm_indexes,1],ID_genes_selected_notna))
  b <- number_genes_selected - a
  c <- length(GOterm_indexes) - a
  d <- length(ID_genes_notselected_notna)- c
  type<-alldata[(which(alldata[,4]==GOterm))[1],6]
  
  matrixes <- rbind(matrixes,c(type,a,b,c,d))
}

colnames(matrixes)<-c("type","a","b","c","d")
rownames(matrixes)<-terms
matrixes<-matrixes[order(rownames(matrixes)),]

# creation of three sub-matrices in function of the GOterm's type  
matrixesCC <- matrixes[which(matrixes[,1]=="CC"),]
indexCC<-which(matrixes[,1]=="CC")
terms_CC<-terms[indexCC]
colnames(matrixesCC)<-c("type","a","b","c","d")
matrixesBP <- matrixes[which(matrixes[,1]=="BP"),]
indexBP<-which(matrixes[,1]=="BP")
terms_BP<-terms[indexBP]
colnames(matrixesBP)<-c("type","a","b","c","d")
matrixesMF <- matrixes[which(matrixes[,1]=="MF"),]
indexMF<-which(matrixes[,1]=="MF")
terms_MF<-terms[indexMF]
colnames(matrixesMF)<-c("type","a","b","c","d")

pval_fisherCC<-fisher_test_matrixes(matrixesCC)
pval_fisherBP<-fisher_test_matrixes(matrixesBP)
pval_fisherMF<-fisher_test_matrixes(matrixesMF)

# correction for multiple testing in fisher

fisher_analysis_BP<-FDR_fisher(pval_fisherBP,terms_BP)
number_terms_annotatedBP<-length(fisher_analysis_BP[[1]])

fisher_analysis_CC<-FDR_fisher(pval_fisherCC,terms_CC)
number_terms_annotatedCC<-length(fisher_analysis_CC[[1]])

fisher_analysis_MF<-FDR_fisher(pval_fisherMF,terms_MF)
number_terms_annotatedMF<-length(fisher_analysis_MF[[1]])

# ATTENIONE: NON BISOGNA AVERE DPLYR IN LIBRERIA! (SOVRASCRIVE SELECT)
vals = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))

annotation_terms_BP<-as.data.frame(annotation_terms(vals, fisher_analysis_BP[[2]]))
colnames(annotation_terms_BP)<-"terms"
annotation_terms_CC<-as.data.frame(annotation_terms(vals, fisher_analysis_CC[[2]]))
colnames(annotation_terms_CC)<-"terms"
annotation_terms_MF<-as.data.frame(annotation_terms(vals, fisher_analysis_MF[[2]]))
colnames(annotation_terms_MF)<-"terms"

#most important == lowest pval

pval_fisherCC_sorted<-sort(fisher_test_matrixes(matrixesCC))
pval_fisherBP_sorted<-sort(fisher_test_matrixes(matrixesBP))
pval_fisherMF_sorted<-sort(fisher_test_matrixes(matrixesMF))

# correction for multiple testing in fisher 

fisher_analysis_BP_sorted<-FDR_fisher(pval_fisherBP_sorted,terms_BP)
fisher_analysis_CC_sorted<-FDR_fisher(pval_fisherCC_sorted,terms_CC)
fisher_analysis_MF_sorted<-FDR_fisher(pval_fisherMF_sorted,terms_MF)


annotation_terms_BP_sorted<-as.data.frame(annotation_terms(vals, fisher_analysis_BP_sorted[[2]]))
colnames(annotation_terms_BP_sorted)<-"terms"
annotation_terms_CC_sorted<-as.data.frame(annotation_terms(vals, fisher_analysis_CC_sorted[[2]]))
colnames(annotation_terms_CC_sorted)<-"terms"
annotation_terms_MF_sorted<-as.data.frame(annotation_terms(vals, fisher_analysis_MF_sorted[[2]]))
colnames(annotation_terms_MF_sorted)<-"terms"

# Fourth part

d <- getGeneLengthAndGCContent(ID_genes_selected, "hsa")

# remove the first column of data (gene names), extract the selected genes and normalize them
data_normalized <- samples
data_normalized <- data_normalized[index_genes_selected,]
data_normalized <- t(t(data_normalized)/d[[1]])

# remove the duplicates and bind them
dataNorm_nodup_control <- remove_duplicates(data_normalized,control) 
dataNorm_nodup_disease <- remove_duplicates(data_normalized,disease) 
dataNorm_clustering <- cbind(dataNorm_nodup_control,dataNorm_nodup_disease)

K <- seq(1,10)

WITHIN_SS_gene_kmeans <- NULL
clus_km <- NULL
s_genes_kmeans <- NULL

for(i in K){
  k <- K[i]
  cl_kmeans_genes <- kmeans(x=dataNorm_clustering,centers=k,iter.max=100,nstart=3)
  clus_km <- c(clus_km,cl_kmeans_genes)
  WITHIN_SS_gene_kmeans <- rbind(WITHIN_SS_gene_kmeans, cl_kmeans_genes$tot.withinss)
  s_genes_kmeans <- rbind(s_genes_kmeans, silhouette(dataNorm_clustering,cl_kmeans_genes[[1]], k))
}

D <- dist(dataNorm_clustering) 
cl_hclust_ward <- hclust(d=D, method="ward.D2")

s_genes_hierar <- NULL

for (i in (1:length(K))){
  k <- K[i]
  clusters_hclust_ward <- cutree(cl_hclust_ward, k=k)
  s_genes_hierar <- c(s_genes_hierar, silhouette(dataNorm_clustering,clusters_hclust_ward,k))
}

WITHIN_SS_sample <- NULL
clus_km_sample <- NULL
s_samples_kmeans <- NULL

for(i in K) {
  k<-K[i]
  cl_kmeans_samples<-kmeans(x=t(dataNorm_clustering),centers=k,iter.max=100,nstart=100)
  clus_km_sample<-c(clus_km_sample,cl_kmeans_samples)
  WITHIN_SS_sample<-rbind(WITHIN_SS_sample, cl_kmeans_samples$tot.withinss)
  s_samples_kmeans <- rbind(s_samples_kmeans, silhouette(t(dataNorm_clustering),cl_kmeans_samples[[1]], k))
}

D<-dist(t(dataNorm_clustering)) 
cl_hclust_ward_S<-hclust(d=D,method="ward.D2")

s_samples_hierar <- NULL

for (i in (1:length(K))){
  k <- K[i]
  clusters_hclust_ward_S<-cutree(cl_hclust_ward_S, k=k)
  s_samples_hierar <- c(s_samples_hierar, silhouette(t(dataNorm_clustering),clusters_hclust_ward_S,k))
}

# results for the gap statistic
#functions for the subsequent analysis
test_hclust <- function(x, k) list(cluster=cutree(hclust(dist(x), method = "average"),k=k))
test_kmeans <- function(x, k) (kmeans(x=x,centers=k,iter.max=100,nstart=100))

# K-Means - genes

gap_res<-clusGap(dataNorm_clustering, test_kmeans, length(K), B=20)

for (i in (2:(nrow(gap_res$Tab)-1))){
  if ((gap_res$Tab[i,3])>(gap_res$Tab[i+1,3]+gap_res$Tab[i+1,4]))
    break;
}
kopt_g_kmeans_gap <- i

# K-Means - samples

gap_res<-clusGap(t(dataNorm_clustering), test_kmeans, length(K), B=20)

for (i in (2:(nrow(gap_res$Tab)-1))){
  if ((gap_res$Tab[i,3])>(gap_res$Tab[i+1,3]+gap_res$Tab[i+1,4]))
    break;
}
kopt_s_kmeans_gap <- i

# Hierarchical - genes

gap_res<-clusGap(dataNorm_clustering, test_hclust, length(K), B=20)

for (i in (2:(nrow(gap_res$Tab)-1))){
  if ((gap_res$Tab[i,3])>(gap_res$Tab[i+1,3]+gap_res$Tab[i+1,4]))
    break;
}
kopt_g_hier_gap <- i

# Hierarchical - samples

gap_res <-clusGap(t(dataNorm_clustering), test_hclust, length(K), B=20)

for (i in (2:(nrow(gap_res$Tab)-1))){
  if ((gap_res$Tab[i,3])>(gap_res$Tab[i+1,3]+gap_res$Tab[i+1,4]))
    break;
}
kopt_s_hier_gap <- i


# Fifth Part

set.seed(3738)

# extract data, remove duplicates and zeros
data_nodup_control <- remove_duplicates(samples[index_genes_selected,], control) 
data_nodup_disease <- remove_duplicates(samples[index_genes_selected,], disease)
groups_nozero <- remove_zeros(data_nodup_control, data_nodup_disease)
data_nodup_control <- groups_nozero$control
data_nodup_disease <- groups_nozero$disease
data_SVM <- cbind(data_nodup_control, data_nodup_disease)
# assign its ID to each gene
rownames(data_SVM) <- ID_genes_selected[-(groups_nozero$removedindexes)]
# construct the vector for control and disease 
namegroup <- c(rep("control", ncol(data_nodup_control)),rep("disease", ncol(data_nodup_disease)))

# separate in train and test, both data and labels
trainIndex <- createDataPartition((1:ncol(data_SVM)), p=0.7, list=FALSE, times=1)
data_train <- data_SVM[,trainIndex]
data_test <- data_SVM[,-trainIndex]

label_train <- namegroup[trainIndex]
label_test <- namegroup[-trainIndex]

# we can still have genes in data_train and data_test for which all values are 0: we remove these genes (if in one group it generates this problem, we have to remove that gene also from the other group).
train_removed <- remove_zeros_onegroup(data_train)
data_train <- train_removed$group
data_test <- data_test[-train_removed$removedindexes,]

# finally, take the transpose to have samples on the rows
data_test<-t(data_test)

dataNorm_train <- scale(data_train)

mean_train <- apply(data_train,2,mean)
sd_train <- apply(data_train,2,sd)
dataNorm_test <- scale(data_test, center = mean_train, scale = sd_train)

# create data.frames in which in the first element we have the factor element
dataNorm_train<-as.data.frame(dataNorm_train)
dataNorm_train<-cbind('Group'=factor(label_train),dataNorm_train)

dataNorm_test<-as.data.frame(dataNorm_test)
dataNorm_test<-cbind('Group'=factor(label_test),dataNorm_test)
data_train<-t(data_train)

rfe <- recursiveFeatureExtractionCV(dataNorm_train, label_train, 500, 10)

best_genes <- (rfe$bestmodel)$names
best_model <- (rfe$bestmodel)$svm
matrix_train <- dataNorm_train[,best_genes]
matrix_train <- cbind('Group'=factor(label_train),matrix_train)

# Performances on the test set 
pred_prova <- predict(best_model, newdata=dataNorm_test, decision.values = FALSE)
res_prova <- confusionMatrix(pred_prova, factor(label_test))  
print(res_prova)

# Check the performances
svmfit_prova <- svm(Group ~ ., data = matrix_train, kernel = "linear", type = 'C-classification', scale = FALSE, na.action = na.omit)
pred_prova <- predict(svmfit_prova, newdata=matrix_train, decision.values = FALSE)
res_prova <- confusionMatrix(pred_prova, factor(label_train))  

accuracy_test <- NULL
par <- NULL
for (i in (1:100)){
  svmfit_prova <- rfe$bests[[i]]
  pred_prova <- predict(svmfit_prova, newdata=dataNorm_test, decision.values = FALSE)
  res_prova <- confusionMatrix(pred_prova, factor(label_test)) 
  accuracy_test <- c(accuracy_test,res_prova[["overall"]][["Accuracy"]])
  par <- c(par,rfe$number_features[i])
}

w <- t(best_model$coefs) %*% best_model$SV
w<-abs(w)
names_sorted <- best_genes[order(w, decreasing=TRUE)]
print(names_sorted)
x <- dataNorm_train[,best_genes]


save(file = "workspace.Rdata")
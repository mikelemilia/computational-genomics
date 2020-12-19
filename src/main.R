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
lambda_est_ttest <- 0.8
eps<-0.1
G0_est_ttest<-G0_value_estimation(lambda_est_ttest, eps, res)

expected_ttest_est <- expected_values(G, G0_est_ttest, alpha, num_sel_ttest)

res <- estimateG0(c_wilcoxon_pvalue, G0, "Wilcoxon test")
lambda_est_wilcoxon <- 0.8
eps<-0.1
G0_est_wilcoxontest<-G0_value_estimation(lambda_est_wilcoxon, eps, res)

expected_wilcoxontest_est <- expected_values(G, G0_est_wilcoxontest, alpha, num_sel_wilcox)

res <- estimateG0(out[,4], G0, "edgeR test")
lambda_est_edger <- 0.8
eps<-0.1
G0_est_edger<-G0_value_estimation(lambda_est_edger, eps, res)

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

#TODO: ho tolto le istruzioni 'unique' nella selezione dei nomi dei geni 
#ma bisogna capire come gestire i nomi doppi con geneID differenti!!!!!

ID_genes_selected<-rownames(DATA[index_genes_selected,])
number_genes_selected<-length(ID_genes_selected)
ID_genes_notselected <- setdiff(rownames(DATA),ID_genes_selected)
number_genes_notselected<-length(ID_genes_notselected)


## ----- GESTIONE DEI NOMI DOPPI CON GENEID DIFFERENTI  --------------
#può andare bene?? NB: I NOMI SONO ANCHE DOPPI, MA GLI ENSEMBL (ENSG...) NON LO SONO
# 14 - estrazione di tutti i GOterm associati ai geni selezionati e creazione delle tabelle associate  

library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)

alldata <- select(org.Hs.eg.db, ID_genes_selected, columns = c("SYMBOL","ENTREZID", "ENSEMBL","GOALL"), keytype="ENSEMBL")
GOALL_NA<-which(is.na(alldata$GOALL))
#which(!(GOALL_NA==alldata$ONTOLOGYALL)) --> si nota che se NA su GOALL allora NA anche su ONTOLOGYALL 
ID_goall_na<-alldata$ENSEMBL[GOALL_NA]
GOALL_NA<-unique(GOALL_NA)

terms <- unique(alldata[,4])
terms <- terms[!is.na(terms)]

#QUESTO SOSTITUISCE I DUE CICLI FOR SUCCESSIVI, DA CAPIRE COSA TOGLIERE
ID_genes_selected_notna <- setdiff(ID_genes_selected,ID_goall_na)
number_genes_selected_notna <- length(ID_genes_selected_notna)
ID_genes_notselected_notna <- setdiff(ID_genes_notselected,ID_goall_na)
number_genes_notselected_notna <- length(ID_genes_notselected_notna)

#for(i in (1:length(GOALL_NA))){
#  term<-names_goall_na[i]
#  j<-0
#  l<-length(names_genes_selected)
#  while(l>0 && j<length(names_genes_selected))
#  {
#    j<-j+1
#    name_sel<-names_genes_selected[j]
#    if (term==name_sel)
#    {
#      names_genes_selected<-names_genes_selected[-j]
#      number_genes_selected<-number_genes_selected-1
#    }
#    l<-l-1
#  }
#}
#for(i in (1:length(GOALL_NA))){
#  term<-names_goall_na[i]
#  j<-0
#  l<-length(names_genes_notselected)
#  while(l>0 && j<length(names_genes_notselected))
#  {
#    j<-j+1
#    name_sel<-names_genes_notselected[j]
#    if (term==name_sel)
#    {
#      names_genes_notselected<-names_genes_notselected[-j]
#      number_genes_notselected<-number_genes_notselected-1
#    }
#    l<-l-1
#  }
#}

matrixes <- NULL

for (i in (1:length(terms))){
  GOterm <- terms[i]
  GOterm_indexes <- which(alldata$GOALL==GOterm)
  a <- length(intersect(alldata[GOterm_indexes,1],ID_genes_selected_notna))
  b <- number_genes_selected - a
  c <- length(GOterm_indexes) - a
  d <- length(ID_genes_notselected_notna)- c
  type<-alldata[(which(alldata[,4]==GOterm))[1],6]
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

# 16 - lunghezza geni e clustering

#library(goseq)
#lengths_genes_selected<-getlength(ID_genes_selected, 'hg19', 'ensGene')
#l<-as.matrix(lengths_genes_selected)

library (EDASeq)
ensembl_list <- ID_genes_selected
d<-getGeneLengthAndGCContent(ensembl_list, "hsa")
#d[[1]] contiene le lunghezz dei geni

data_normalized<-DATA[,-1]
data_normalized<-data_normalized[index_genes_selected,]
data_normalized<-t(t(data_normalized)/d[[1]])

# da qui in avanti non ci sono gli oggetti nel dataset salvato

# 17 - clustering

## IDEA
# FACENDO CLUSTERING DEI GENES MI ASPETTO CHE IL CLUSTER SIA 1 O COMUNQUE POCHI CLUSTERS.... SE SONO TUTTI 
# I SELEZIONATI NON DOVREI VEDERE DIFFERENZE TRA LORO, SONO TUTTI ASSOCIATI ALLA MALATTIA

# FACENDO CLUSTERING DEI SAMPLES, DOVREI VEDERE CHE SI DIVIDONO IN DUE CLASSI, MALATI E SANI
# BISOGNA RICONSIDERARE I SAMPLES CAMPIONATI DUE VOLTE?


dataNorm_nodup_control<-remove_duplicates(data_normalized,control) 
dataNorm_nodup_disease<-remove_duplicates(data_normalized,disease) 

dataNorm_clustering<-cbind(dataNorm_nodup_control,dataNorm_nodup_disease)


#------------------- HIERARCHICAL CLUSTERING ----------------------
#CLUSTERING GENES 
D<-dist(dataNorm_clustering) #D is an object of class "dist". To get a matrix one needs to use "as.matrix(D)"


# USIAMO WARD PERCHè SFRUTTA DISTANZA EUCLIDEA E CPSì LO COMPARIAMO BENE CON KMEANS CHE USA SEMPRE LA EUCLIDEA
cl_hclust_ward<-hclust(d=D,method="ward.D2")
plot(cl_hclust_ward, hang=-1) 

sk <- NULL
K<-seq(1,10,by=1)
for (i in (1:length(K))){
  k <- K[i]
  clusters_hclust_ward<-cutree(cl_hclust_ward, k=k)
  sk <- c(sk,silhouette(dataNorm_clustering,clusters_hclust_ward,k))
  #print("i: ",i," - sk: ",sk)
}
cat("Hierarchial clusters over genes!\n")
print(sk)
cat(max(sk), " - optimal number of clusters is :", K[which(sk == max(sk))])

#CLUSTERING SAMPLES
D<-dist(t(dataNorm_clustering))  #D is an object of class "dist". To get a matrix one needs to use "as.matrix(D)"

cl_hclust_ward_S<-hclust(d=D,method="ward.D2")
plot(cl_hclust_ward_S, hang=-1) 

sk <- NULL
K<-seq(1,10,by=1)
for (i in (1:length(K))){
  k <- K[i]
  clusters_hclust_ward_S<-cutree(cl_hclust_ward_S, k=k)
  sk <- c(sk,silhouette(t(dataNorm_clustering),clusters_hclust_ward_S,k))
}
cat("Hierarchial clusters over samples!\n")
print(sk)
cat(max(sk), " - optimal number of clusters is :", K[which(sk == max(sk))])

#------------------- k MEANS ----------------------
#CLUSTERING GENES 

K<-seq(1,3,by=1)
WITHIN_SS<-NULL
clus_km<-NULL
sk <- NULL
for(i in (1:length(K)))
{
  k_i<-K[i]
  cl_kmeans_genes<-kmeans(x=dataNorm_clustering,centers=k_i,iter.max=100,nstart=1)
  clus_km<-c(clus_km,cl_kmeans_genes)
  WITHIN_SS<-rbind(WITHIN_SS, cl_kmeans_genes$tot.withinss)
  sk <- rbind(sk, silhouette(dataNorm_clustering,cl_kmeans_genes[[1]], k_i))
}
print(sk)
cat("K-Means over samples!\n")
print(sk)
cat(max(sk), " - optimal number of clusters is :", K[which(sk == max(sk))])
plot(K, WITHIN_SS)

#CLUSTERING SAMPLES
K<-seq(1,10,by=1)
WITHIN_SS_sample<-NULL
clus_km_sample<-NULL
sk <- NULL
for(i in (1:length(K))) {
  k_i<-K[i]
  cl_kmeans_samples<-kmeans(x=t(dataNorm_clustering),centers=k_i,iter.max=100,nstart=100)
  clus_km_sample<-c(clus_km_sample,cl_kmeans_samples)
  WITHIN_SS_sample<-rbind(WITHIN_SS_sample, cl_kmeans_samples$tot.withinss)
  sk <- rbind(sk, silhouette(t(dataNorm_clustering),cl_kmeans_samples[[1]], k_i))
}
cat("K-Means over samples!\n")
print(sk)
cat(max(sk), " - optimal number of cluster is :", K[which(sk == max(sk))])
plot(K, WITHIN_SS_sample)


#GAP STATISTIC 
library(cluster)

test_hclust <- function(x, k) list(cluster=cutree(hclust(dist(x), method = "average"),k=k))

#SAMPLES
prova<-clusGap(t(dataNorm_clustering), test_hclust, length(K), B=100)

for (i in (2:(nrow(prova[[1]])-1))){
  if (prova[[1]][i,3]>prova[[1]][i+1,3]+prova[[1]][i+1,4])
    break;
}
cat("Numero ottimo per hierachical clustering dei samples secondo Gap Statistics: ",i)

#GENES
prova<-clusGap(dataNorm_clustering, test_hclust, length(K), B=100)

for (i in (1:(nrow(prova[[2]])-1))){
  if (prova[[1]][i,3]>prova[[1]][i+1,3]+prova[[1]][i+1,4])
    break;
}
cat("Numero ottimo per hierachical clustering dei geni secondo Gap Statistics: ",i)

# ---------------

#SAMPLES

test_kmeans <- function(x, k) (kmeans(x=x,centers=k,iter.max=100,nstart=100))

prova<-clusGap(t(dataNorm_clustering), test_kmeans, length(K), B=100)

for (i in (2:(nrow(prova[[1]])-1))){
  if ((prova[[1]][i,3])>(prova[[1]][i+1,3]+prova[[1]][i+1,4]))
    break;
}
cat("Numero ottimo per kmeans dei samples secondo Gap Statistics: ",i)

#GENES

prova<-clusGap(dataNorm_clustering, test_kmeans, length(K), B=2)

for (i in (2:(nrow(prova[[1]])-1))){
  if ((prova[[1]][i,3])>(prova[[1]][i+1,3]+prova[[1]][i+1,4]))
    break;
}
cat("Numero ottimo per kmeans dei geni secondo Gap Statistics: ",i)









## PLOTS-----------------------------------------------------------------
#kmeans
#genes
library(factoextra)
cl_genes<-kmeans(x=dataNorm_clustering, centers = 2, iter.max = 100, nstart = 100)
fviz_cluster(cl_genes, data = dataNorm_clustering,
             palette = c("#2E9FDF", "#00AFBB"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#samples
cl_samples<-kmeans(x=t(dataNorm_clustering), centers = 3, iter.max = 100, nstart = 100)
fviz_cluster(cl_samples, data = t(dataNorm_clustering),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#hierarchical
#samples
fviz_dend(cl_hclust_ward_S, rect = TRUE)
#----oppure------
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")

hcd_S <- as.dendrogram(cl_hclust_ward_S)
plot(hcd_S, type = "rectangle", horiz = FALSE, xlab = "Height", nodePar = nodePar, leaflab = "none")

#genes
fviz_dend(cl_hclust_ward, rect = TRUE)
#----oppure------
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.1, col = "blue")

hcd_G <- as.dendrogram(cl_hclust_ward)
plot(hcd_G, type = "rectangle", horiz = FALSE, ylab = "Height", nodePar = nodePar, leaflab = "none")





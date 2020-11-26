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

  plot(normed[[2]][,i], normed[[1]][,i], ylab="M",xlab="A",main=paste("MvA Plot", i-2, sep =" "))
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

#Rebuilt control matrix
control_nodup<-remove_duplicates(DATA,control)

nomi<-rep("", count)
for (i in (1:count))
{
  x<-paste("control",as.character(i),sep='_')
  nomi[i]<-x
}
colnames(control_nodup)<-nomi

#Rebuilt disease matrix

disease_nodup<-remove_duplicates(DATA,disease)

nomi<-rep("", count)
for (i in (1:count))
{
  x<-paste("disease",as.character(i),sep='_')
  nomi[i]<-x
}
colnames(disease_nodup)<-nomi

#------------------------- edgeR-----------------------------------
#PROVA DI MARTI: ho provato a implementare edgeR; i gruppi sono due, quindi ho usato il primo
#esempio di implemetazione di edgeR della professoressa (L11-EdgeR.Rmd, righe 96-158)

#tutte le variabili hanno nomi di prova, poi se va bene possiamo cambiarle

#unisco le due tabelle
#sono diverse rispetto a quelle create prima perchè mi serve che abbiano tutti nomi diversi nelle colonne
#del tipo "control_1", "control_2" ecc ecc

prova <- cbind(disease_nodup, control_nodup)
library(edgeR)

# creo i gruppi che poi mi serviranno per definire "group"; di base mi serve solo che mi separino
#ciò che c'è in tabella in control e disease
#sono due vettori
gruppo_controllo <- rep("control",dim(control_nodup)[2])
gruppo_malato <- rep("disease",dim(disease_nodup)[2])

#li unisco, mi serve vettore unico per usare factor
gruppo <- cbind(t(as.data.frame(gruppo_malato)),t(as.data.frame(gruppo_controllo)))

#uso factor(), ottengo un oggetto diviso in due livelli (contorl e disease), come ci serve
groupProva <- factor(gruppo)

#matrice di design
designprova <- model.matrix(~0+groupProva) 
rownames(designprova) <- colnames(prova)  
print(designprova)


# fit values of phi (we need this step to fit our GLM model)
yprova <- DGEList(counts=prova, remove.zeros = TRUE)    # y is an object of type DGE
yprova <- calcNormFactors(yprova)   # This calculates the SF using the TMM normalization !!!
SF<-yprova$samples

y <- estimateGLMCommonDisp(yprova,designprova, verbose=TRUE) #phi common to the entire dataset
y <- estimateGLMTrendedDisp(y,designprova) #phi depends on mu
y <- estimateGLMTagwiseDisp(y,designprova) #phi is gene specific
fit <- glmFit(y,designprova) #finally the model fit (that accounts for raw NB data and scaling factors and seq. depth) 
summary(fit)

#il test
Confronti<-makeContrasts(Treatment=groupProvadisease-groupProvacontrol,levels=designprova)
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

#selected usando edgeR e il suo FDR... credo non possiamo, ma giusto per capire
indSELedgeR<-which(out$FDR<0.05) #i selected
print(length(indSELedgeR))

# 11 - E[FP] and E[FN] for t-test and Wilcoxon test with G0 = G

G<-nrow(control_nodup_nozero)
G0<-G

#function that returns a vector with, in order, TP FP TN FN
expected_ttest <- expected_values(G, G0, alpha, num_sel_ttest)
expected_wilcoxontest <- expected_values(G, G0, alpha, num_sel_wilcox)
expected_edger <- expected_values(G, G0, alpha, indSELedgeR)

# 12 - Estimate G0 and re-estimate FP and FN

res <- estimateG0(c_ttest_pvalue)
G0_est_ttest <- res[[1]]
lambda_est <- res[[2]]

expected_ttest_est <- expected_values(G, G0_est_ttest, alpha, num_sel_ttest)

res <- estimateG0(c_wilcoxon_pvalue)
G0_est_wilcoxontest <- res[[1]]
lambda_est <- res[[2]]

expected_wilcoxontest_est <- expected_values(G, G0_est_wilcoxontest, alpha, num_sel_wilcox)

res <- estimateG0(out[,5])
G0_est_edger <- res[[1]]
lambda_est <- res[[2]]

expected_edgeR_est <- expected_values(G, G0_est_edger, alpha, num_sel_wilcox)

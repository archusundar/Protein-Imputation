library(corrr)
library(caret)
library(ggplot2)
  setwd('/Users/archana/Downloads/R_Code_for_package/')
# Loading the data

mRNA <- readRDS("mRNA.rds")
protein_vector<-readRDS("protein.rds")
model<-readRDS("model.rds")
test<-readRDS("test.rds")

# Functions
# This function is to plot raw mRNA and protein correlations; 
#This first column of mRNA_dataframe must be mRNA names
# The protein name is given for x-axis legent 
mRNA_name = protein_name ='BDNF'
Rho_cutoff=0.75
raw_plot <- function (mRNA, protein_vector, mRNA_name, protein_name) {
  mRNA<-mRNA[,mRNA_name]
  ptn<-protein_vector[,protein_name]
  data=data.frame(mRNA, ptn)
  cor<-cor(mRNA,ptn,method="spearman")
  cor<-round(cor,3)
  ggplot(data, aes(x = mRNA, y = ptn)) + xlab("mRNA(Normalized TPMs)") + ylab("Protein(pg/mg)")+  
    geom_point()  + ggtitle(paste0(mRNA_name,"\ncorrelation =", cor)) +
    theme_bw() + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                                    axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
                                    axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                                    axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
                                    plot.title = element_text(size = 25,  hjust = 0.5))
  #dev.off()

}



# This function is to reduce multicollinearity
collinearityReduce <- function (mRNA, Rho_cutoff) {
  gene2 <- mRNA[apply(mRNA!=0, 2, all),]
  col.names(gene2)<-gene2[,1]
  gene2<-gene2[,-1]
  
  gene3<-data.frame(t(mRNA))
  row.names(gene3)<-colnames(mRNA)
  colnames(gene3)<-rownames(mRNA)
  gene2 <- gene3[apply(gene3!=0, 1, all),]
  gene2<-as.matrix(gene2)
  data_cor <- cor(gene2, method = "spearman")
  data_highcor1 <- findCorrelation(data_cor, cutoff=Rho_cutoff)
  data_out<-gene3[,-data_highcor1]
  
}


# This function is to optimize the imputing model parameters via cross-validations
trainOpt <- function (data_out, Rho_cutoff, name) {
  #ptn.rds <- train(BDNF~., data=data_out, method = "svmRadial")
  data1<-data.frame(data_out,protein_vector)
  ctrl <- trainControl(method = "CV", number = 5, savePred=T) 
  ptn.rds <- train(protein_name~., data=data1, method = "svmRadial", trControl = ctrl)
  saveRDS(ptn.rds,"model.rds")
}



# This function is to
# Return Rho (obs vs imputed)

modelPerformance <- function(model,test) {
  model<-readRDS(model)
  test<-readRDS(test)
  x1<-test[protein_name]
  y1<-predict(model,test)
  cor<-cor(x1,y1,method="spearman")
  cor<-round(cor,3)
  cat("correlation=",cor)
}

# This function is for independent dataset testing
# Return imputed values

modelPred <- function (model, test) {
    model<-readRDS(model)
    test<-readRDS(test)
    x1<-test[protein_name]
    y1<-predict(model,test)
    cor<-cor(x1,y1,method="spearman")
    cor<-round(cor,3)
    cat("correlation=",cor)
    data=data.frame(x1, y1)
    names(data)<-c('mRNA','ptn')
    ggplot(data, aes(x = mRNA, y = ptn)) + xlab("Observed (pg/mg)") + ylab("Predicted (pg/mg)")+  
      geom_point()  + ggtitle(paste0(protein_name,"\ncorrelation =", cor)) +
      theme_bw() + theme(axis.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                         axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
                         axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                         axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
                         plot.title = element_text(size = 25,  hjust = 0.5))
  
}
raw_plot(mRNA,protein_vector,mRNA_name,protein_name)
collinearityReduce(mRNA, Rho_cutoff)
modelPerformance('model.rds', 'test.rds')
modelPred('model.rds', 'test.rds')

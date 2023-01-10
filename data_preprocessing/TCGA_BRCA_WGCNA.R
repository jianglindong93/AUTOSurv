library("WGCNA")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

setwd("D:/DL/TCGA data/TCGA-GDC-BRCA/3rd_pipeline/train_test_split/no_normalization")
# genes
data_gene<-read.csv("data_train_gene_overall.csv", h = T)
datExpr<-data_gene[,-c(1:7)]
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
x11()
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

net = blockwiseModules(datExpr, power = 5, maxBlockSize = 2699,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TCGABRCATOM_gene",
                       verbose = 3)
table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "TCGABRCA_gene-01-networkConstruction-auto.RData")
rownames(MEs)<-c(1:nrow(MEs))

data_gene_test<-read.csv("data_test_gene_overall.csv", h = T)
datExpr_test<-data_gene_test[,-c(1:7)]

MEs_train = moduleEigengenes(datExpr, moduleColors)$eigengenes
data_train_eigengene<-cbind(data_gene[,c(1:7)], MEs_train)
write.csv(data_train_eigengene, file = "data_train_eigengene_overall.csv", row.names = F)
MEs_test = moduleEigengenes(datExpr_test, moduleColors)$eigengenes
data_test_eigengene<-cbind(data_gene_test[, c(1:7)], MEs_test)
write.csv(data_test_eigengene, file = "data_test_eigengene_overall.csv", row.names = F)

# miRNA
data_mirna<-read.csv("data_train_mirna_overall.csv", h = T)
datExpr_mirna<-data_mirna[,-c(1:7)]
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft_mirna = pickSoftThreshold(datExpr_mirna, powerVector = powers, verbose = 5)
x11()
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft_mirna$fitIndices[,1], -sign(sft_mirna$fitIndices[,3])*sft_mirna$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft_mirna$fitIndices[,1], -sign(sft_mirna$fitIndices[,3])*sft_mirna$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft_mirna$fitIndices[,1], sft_mirna$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_mirna$fitIndices[,1], sft_mirna$fitIndices[,5], labels=powers, cex=cex1,col="red")

net_mirna = blockwiseModules(datExpr_mirna, power = 3, maxBlockSize = 516,
                             TOMType = "unsigned", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = "TCGABRCATOM_mirna",
                             verbose = 3)
table(net_mirna$colors)

moduleLabels2 = net_mirna$colors
moduleColors2 = labels2colors(net_mirna$colors)
MEs2 = net_mirna$MEs
geneTree2 = net_mirna$dendrograms[[1]]
save(MEs2, moduleLabels2, moduleColors2, geneTree2,
     file = "TCGABRCA_mirna-01-networkConstruction-auto.RData")
rownames(MEs2)<-c(1:nrow(MEs2))

data_mirna_test<-read.csv("data_test_mirna_overall.csv", h = T)
datExpr_mirna_test<-data_mirna_test[,-c(1:7)]

MEs2_train = moduleEigengenes(datExpr_mirna, moduleColors2)$eigengenes
data_train_eigenmirna<-cbind(data_mirna[,c(1:7)], MEs2_train)
write.csv(data_train_eigenmirna, file = "data_train_eigenmirna_overall.csv", row.names = F)
MEs2_test = moduleEigengenes(datExpr_mirna_test, moduleColors2)$eigengenes
data_test_eigenmirna<-cbind(data_mirna_test[, c(1:7)], MEs2_test)
write.csv(data_test_eigenmirna, file = "data_test_eigenmirna_overall.csv", row.names = F)

# divide and export data
export_matched_data<-function(root_path, file_names, data_train_gene, data_train_mirna){
  for(folder in file_names){
    setwd(paste(paste(root_path, folder, sep = "/"), "no_normalization", sep = "/"))
    data_match_train<-read.csv(paste(paste("data_train", folder, sep="_"),".csv",sep=""), h = T)
    data_match_valid<-read.csv(paste(paste("data_valid", folder, sep="_"),".csv",sep=""), h = T)
    matched_id_train<-match(data_match_train$patient_id, data_train_gene$patient_id)
    matched_id_valid<-match(data_match_valid$patient_id, data_train_gene$patient_id)
    write.csv(data_train_gene[matched_id_train,], file = paste(paste("data_train_eigengene", folder, sep="_"),".csv",sep=""), row.names = F)
    write.csv(data_train_gene[matched_id_valid,], file = paste(paste("data_valid_eigengene", folder, sep="_"),".csv",sep=""), row.names = F)
    write.csv(data_train_mirna[matched_id_train,], file = paste(paste("data_train_eigenmirna", folder, sep="_"),".csv",sep=""), row.names = F)
    write.csv(data_train_mirna[matched_id_valid,], file = paste(paste("data_valid_eigenmirna", folder, sep="_"),".csv",sep=""), row.names = F)
  }
}

file_names_export<-paste("tune", seq(1:10), sep = "_")
root_dir<-"D:/DL/TCGA data/TCGA-GDC-BRCA/3rd_pipeline"

export_matched_data(root_dir, file_names_export, data_train_eigengene, data_train_eigenmirna)

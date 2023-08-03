setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")

# read and check data status
ov_au_seqexp<-read.table("exp_seq.OV-AU.tsv", header = T, sep = "\t")
colnames(ov_au_seqexp)
unique(ov_au_seqexp$gene_model)
length(unique(ov_au_seqexp$icgc_donor_id))
table_ids<-matrix(table(ov_au_seqexp$icgc_donor_id))
unique(table_ids[,1])

# change data display format
divide_donor_data<-function(data_raw, var_names){
  donor_list<-list()
  for(i in unique(data_raw$icgc_donor_id)){
    donor_index<-grep(i, data_raw$icgc_donor_id, fixed = T)
    donor_list[[i]]<-data_raw[donor_index, var_names]
  }
  return(donor_list)
}

var_names_run<-c("icgc_donor_id", "icgc_specimen_id",
                 "icgc_sample_id", "submitted_sample_id",
                 "analysis_id", "gene_id",
                 "normalized_read_count")

donor_data<-divide_donor_data(ov_au_seqexp, var_names_run)

get_rep_sample_names<-function(donor_list){
  rep_specimen_names<-list()
  rep_sample_names<-list()
  rep_subsample_names<-list()
  for(i in names(donor_list)){
    if(length(unique(donor_list[[i]]$icgc_specimen_id))>1){
      rep_specimen_names[[i]]<-unique(donor_list[[i]]$icgc_specimen_id)
    }
    if(length(unique(donor_list[[i]]$icgc_sample_id))>1){
      rep_sample_names[[i]]<-unique(donor_list[[i]]$icgc_sample_id)
    }
    if(length(unique(donor_list[[i]]$submitted_sample_id))>1){
      rep_subsample_names[[i]]<-unique(donor_list[[i]]$submitted_sample_id)
    }
  }
  rep_names<-list()
  rep_names[["specimen"]]<-rep_specimen_names
  rep_names[["sample"]]<-rep_sample_names
  rep_names[["subsample"]]<-rep_subsample_names
  return(rep_names)
}

rep_names<-get_rep_sample_names(donor_data)
all.equal(names(rep_names[["specimen"]]), names(rep_names[["sample"]]))
rep_specimen<-rep_names[["specimen"]]

get_corr_rep_samples<-function(donor_list, rep_name_list){
  var_names<-c("gene_id", "normalized_read_count")
  corr_cache<-list()
  for(i in names(rep_name_list)){
    print(i)
    temp_corr<-c()
    rep_name_1<-rep_name_list[[i]][1]
    temp_data_1<-donor_list[[i]][grep(rep_name_1, donor_list[[i]]$icgc_specimen_id), var_names]
    colnames(temp_data_1)[2]<-"normalized_read_count_1"
    for(j in 2:length(rep_name_list[[i]])){
      rep_name_c<-rep_name_list[[i]][j]
      temp_data_c<-donor_list[[i]][grep(rep_name_c, donor_list[[i]]$icgc_specimen_id), var_names]
      colnames(temp_data_c)[2]<-"normalized_read_count_c"
      temp_data_comb<-merge(temp_data_1, temp_data_c, by = "gene_id")
      if(anyNA(temp_data_comb)){
        print("Warning! Merging did not go well.")
      }
      temp_corr<-c(temp_corr, 
                   cor(temp_data_comb$normalized_read_count_1, 
                       temp_data_comb$normalized_read_count_c,
                       method = "spearman"))
    }
    corr_cache[[i]]<-temp_corr
  }
  return(corr_cache)
}

get_corr_test_rep_samples<-function(donor_list, rep_name_list){
  var_names<-c("gene_id", "normalized_read_count")
  corr_test_cache<-list()
  for(i in names(rep_name_list)){
    print(i)
    temp_corr_test<-c()
    rep_name_1<-rep_name_list[[i]][1]
    temp_data_1<-donor_list[[i]][grep(rep_name_1, donor_list[[i]]$icgc_specimen_id), var_names]
    colnames(temp_data_1)[2]<-"normalized_read_count_1"
    for(j in 2:length(rep_name_list[[i]])){
      rep_name_c<-rep_name_list[[i]][j]
      temp_data_c<-donor_list[[i]][grep(rep_name_c, donor_list[[i]]$icgc_specimen_id), var_names]
      colnames(temp_data_c)[2]<-"normalized_read_count_c"
      temp_data_comb<-merge(temp_data_1, temp_data_c, by = "gene_id")
      if(anyNA(temp_data_comb)){
        print("Warning! Merging did not go well.")
      }
      temp_corr_test<-c(temp_corr_test, 
                        cor.test(temp_data_comb$normalized_read_count_1, 
                                 temp_data_comb$normalized_read_count_c,
                                 method = "kendall")$p.value)
    }
    corr_test_cache[[i]]<-temp_corr_test
  }
  return(corr_test_cache)
}

corr_cache_run<-get_corr_rep_samples(donor_data, rep_specimen)

#corr_test_cache_run<-get_corr_test_rep_samples(donor_data, rep_specimen)

reform_exp_data<-function(donor_list, rep_name_list, path){
  setwd(path)
  specimen_ids<-c()
  combinedData<-c()
  var_names<-c("gene_id", "normalized_read_count")
  for(i in names(donor_list)){
    if(i %in% names(rep_name_list)){
      specimen_ids<-c(specimen_ids, rep_name_list[[i]][1])
      temp_data<-donor_list[[i]][grep(rep_name_list[[i]][1], donor_list[[i]]$icgc_specimen_id), var_names]
    }else{
      specimen_ids<-c(specimen_ids, i)
      temp_data<-donor_list[[i]][, var_names]
    }
    colnames(temp_data)[2]<-i
    if(length(combinedData)==0){
      combinedData<-temp_data
    }else{
      combinedData<-merge(combinedData, temp_data, by = "gene_id")
    }
  }
  specimen_id_cache<-data.frame(donor = names(donor_list), specimen = specimen_ids)
  write.csv(specimen_id_cache, file = "specimens_used.csv", row.names = F)
  return(combinedData)
}

reform_exp_data_alt<-function(donor_list, rep_name_list, corr_test_cache, path){
  setwd(path)
  specimen_ids<-c()
  combinedData<-c()
  var_names<-c("gene_id", "normalized_read_count")
  for(i in names(donor_list)){
    if(i %in% names(rep_name_list)){
      if(max(corr_test_cache[[i]])>0.05){
        specimen_ids<-c(specimen_ids, "mean")
        combinedData_sub<-c()
        for(j in rep_name_list[[i]]){
          temp_data_sub<-donor_list[[i]][grep(j, donor_list[[i]]$icgc_specimen_id), var_names]
          colnames(temp_data_sub)[2]<-paste("normalized_read_count", j, sep = "_")
          if(length(combinedData_sub)==0){
            combinedData_sub<-temp_data_sub
          }else{
            combinedData_sub<-merge(combinedData_sub, temp_data_sub, by = "gene_id")
          }
        }
        temp_data<-data.frame(gene_id = combinedData_sub$gene_id,
                              normalized_read_count = rowMeans(combinedData_sub[, -1]))
      }else{
        specimen_ids<-c(specimen_ids, rep_name_list[[i]][1])
        temp_data<-donor_list[[i]][grep(rep_name_list[[i]][1], donor_list[[i]]$icgc_specimen_id), var_names]
      }
    }else{
      specimen_ids<-c(specimen_ids, i)
      temp_data<-donor_list[[i]][, var_names]
    }
    colnames(temp_data)[2]<-i
    if(length(combinedData)==0){
      combinedData<-temp_data
    }else{
      combinedData<-merge(combinedData, temp_data, by = "gene_id")
    }
  }
  specimen_id_cache<-data.frame(donor = names(donor_list), specimen = specimen_ids)
  write.csv(specimen_id_cache, file = "gene_specimens_used.csv", row.names = F)
  return(combinedData)
}

reformed_exp_data<-reform_exp_data(donor_data, rep_specimen, path = "D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")

#reformed_exp_data<-reform_exp_data_alt(donor_data, rep_specimen, corr_test_cache_run, path = "D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related_cor_test_version")

length(unique(reformed_exp_data$gene_id))
length(unique(colnames(reformed_exp_data)))
getwd()
write.csv(reformed_exp_data, file = "ov_au_reformed_exp_data.csv", row.names = F)

length(grep(".", reformed_exp_data$gene_id, fixed = T))

# get pathway information
setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related/gene_name_batches")
divide_index<-seq(1, nrow(reformed_exp_data), 6000)
gene_id_divided<-list()
ensembl_list<-reformed_exp_data$gene_id
for(i in 1:length(divide_index)){
  batch_name<-paste("gene_id_batch", i, sep = "_")
  if(i<length(divide_index)){
    gene_id_divided[[batch_name]]<-ensembl_list[c(divide_index[i]:(divide_index[i+1]-1))]
  }
  else{
    gene_id_divided[[batch_name]]<-ensembl_list[c(divide_index[i]:nrow(reformed_exp_data))]
  }
}
batch_name_cache<-paste("gene_id_batch", seq(1,length(gene_id_divided),1), sep = "_")
for(i in 1:length(gene_id_divided)){
  write.csv(gene_id_divided[[batch_name_cache[i]]], file = paste(batch_name_cache[i], ".csv", sep = ""), row.names = F)
}

merge_cache<-function(data_cache){
  merged_data<-c()
  for(i in names(data_cache)){
    merged_data<-rbind(merged_data, as.data.frame(data_cache[[i]]))
  }
  return(merged_data)
}

setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related/pathway_batches")
pathway_batch_cache<-list()
pathway_batch_names<-paste("pathway_batch", seq(1,length(gene_id_divided),1), sep = "_")
for(i in pathway_batch_names){
  pathway_batch_cache[[i]]<-read.csv(paste(i, ".csv", sep = ""), h = T)
}

pathway_data_merged<-merge_cache(pathway_batch_cache)
length(pathway_data_merged$ID)
length(unique(pathway_data_merged$ID))

which(table(pathway_data_merged$ID)>1, arr.ind = T)
length(grep("ENSG00000115239", pathway_batch_cache[[1]]$ID))
grep("ENSG00000115239", pathway_batch_cache[[1]]$ID)
pathway_batch_cache[[1]][458,]
pathway_batch_cache[[1]][1221,]
pathway_batch_cache[[1]]<-pathway_batch_cache[[1]][-458,]

length(grep("ENSG00000137843", pathway_batch_cache[[2]]$ID))
grep("ENSG00000137843", pathway_batch_cache[[2]]$ID)
pathway_batch_cache[[2]][167,]
pathway_batch_cache[[2]][2419,]
pathway_batch_cache[[2]]<-pathway_batch_cache[[2]][-167,]

length(grep("ENSG00000145979", pathway_batch_cache[[2]]$ID))
grep("ENSG00000145979", pathway_batch_cache[[2]]$ID)
pathway_batch_cache[[2]][824,]
pathway_batch_cache[[2]][826,]
pathway_batch_cache[[2]]<-pathway_batch_cache[[2]][-826,]

length(grep("ENSG00000203668", pathway_batch_cache[[4]]$ID))
grep("ENSG00000203668", pathway_batch_cache[[4]]$ID)
pathway_batch_cache[[4]][28,]
pathway_batch_cache[[4]][393,]
pathway_batch_cache[[4]]<-pathway_batch_cache[[4]][-393,]

length(grep("ENSG00000213999", pathway_batch_cache[[4]]$ID))
grep("ENSG00000213999", pathway_batch_cache[[4]]$ID)
pathway_batch_cache[[4]][16,]
pathway_batch_cache[[4]][340,]
pathway_batch_cache[[4]]<-pathway_batch_cache[[4]][-16,]

length(grep("ENSG00000261740", pathway_batch_cache[[9]]$ID))
grep("ENSG00000261740", pathway_batch_cache[[9]]$ID)
pathway_batch_cache[[9]][9,]
pathway_batch_cache[[9]][10,]
pathway_batch_cache[[9]]<-pathway_batch_cache[[9]][-10,]

get_path_mask_r<-function(path_data){
  path_list<-strsplit(path_data[, 4], split = ",")
  path_row<-list()
  path_row_temp<-c()
  for(i in 1:length(path_list)){
    for(j in 1:length(path_list[[i]])){
      if(length(grep("~", path_list[[i]][j], fixed = T)) != 0){
        path_unit<-strsplit(path_list[[i]][j], split = "~")[[1]][1]
        path_row_temp<-c(path_row_temp, path_unit)
      }
    }
    path_row[[i]]<-path_row_temp
    path_row_temp<-c()
  }
  path_comp_list<-unique(unlist(path_row))
  pathway_mask<-matrix(rep(0, nrow(path_data)*length(path_comp_list)), nrow = nrow(path_data))
  for(i in 1:nrow(pathway_mask)){
    pathway_mask[i, which(path_comp_list %in% path_row[[i]])]<-1
  }
  pathway_mask<-as.data.frame(pathway_mask, row.names = path_data[, 1])
  colnames(pathway_mask)<-path_comp_list
  return(pathway_mask)
}

filter_pathway_index<-function(mask){
  ind_remove<-c()
  for (i in 1:ncol(mask)) {
    temp_sum<-sum(mask[, i])
    if(temp_sum < 15 | temp_sum > 300){
      ind_remove<-c(ind_remove, i)
    }
  }
  return(ind_remove)
}

pathway_data_merged<-merge_cache(pathway_batch_cache)
length(pathway_data_merged$ID)
length(unique(pathway_data_merged$ID))

reactome_pathway_mask<-get_path_mask_r(pathway_data_merged)
removable_pathway_index<-filter_pathway_index(reactome_pathway_mask)
length(removable_pathway_index)
reactome_pathway_mask<-reactome_pathway_mask[,-removable_pathway_index]
test_col_sum<-colSums(reactome_pathway_mask)
min(test_col_sum)
max(test_col_sum)

row_sum_genes<-rowSums(reactome_pathway_mask)
removable_gene_index<-which(row_sum_genes==0,arr.ind = T)
reactome_pathway_mask<-reactome_pathway_mask[-removable_gene_index,]

test_index_gene<-rowSums(reactome_pathway_mask)
min(test_index_gene)

data_gene_pathway<-reformed_exp_data[match(rownames(reactome_pathway_mask), reformed_exp_data$gene_id),]
all.equal(data_gene_pathway$gene_id, rownames(reactome_pathway_mask))

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_to_chromosomes_2nd<-getBM(filters = "ensembl_gene_id",
                                attributes = c("ensembl_gene_id", "chromosome_name"),
                                values = data_gene_pathway$gene_id,
                                mart = mart)
head(genes_to_chromosomes_2nd)
length(unique(genes_to_chromosomes_2nd$chromosome_name))
unique(genes_to_chromosomes_2nd$chromosome_name)
removable_chro<-c("Y", "MT")
removable_chro_index<-which(genes_to_chromosomes_2nd$chromosome_name %in% removable_chro, arr.ind = T)
genes_to_chro_shrinked<-genes_to_chromosomes_2nd[-removable_chro_index,]
length(unique(genes_to_chro_shrinked$chromosome_name))
unique(genes_to_chro_shrinked$chromosome_name)

reactome_pathway_mask_shrinked<-reactome_pathway_mask[match(genes_to_chro_shrinked$ensembl_gene_id, rownames(reactome_pathway_mask)),]
all.equal(rownames(reactome_pathway_mask_shrinked), genes_to_chro_shrinked$ensembl_gene_id)
removable_pathway_index_2<-which(colSums(reactome_pathway_mask_shrinked)==0, arr.ind = T)
length(removable_pathway_index_2)
min(colSums(reactome_pathway_mask_shrinked))
setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")
#setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related_cor_test_version")
write.csv(reactome_pathway_mask_shrinked, file = "pathway_mask_general.csv")

data_gene_pathway_shrinked<-data_gene_pathway[match(rownames(reactome_pathway_mask_shrinked), data_gene_pathway$gene_id),]
all.equal(data_gene_pathway_shrinked$gene_id, rownames(reactome_pathway_mask_shrinked))

exp_matrix<-as.matrix(data_gene_pathway_shrinked[, -1])
patient_id<-colnames(data_gene_pathway_shrinked)[-1]
head(patient_id)
exp_matrix<-t(exp_matrix)
colnames(exp_matrix)<-data_gene_pathway_shrinked$gene_id
rownames(exp_matrix)<-c(1:nrow(exp_matrix))
data_exp_pathway_t<-data.frame(cbind(patient_id, exp_matrix))
all.equal(as.numeric(data_exp_pathway_t$ENSG00000006210), as.numeric(data_gene_pathway_shrinked[grep("ENSG00000006210", data_gene_pathway_shrinked$gene_id), -1]))
write.csv(data_exp_pathway_t, file = "ov_au_exp_data_run.csv", row.names = F)

# miRNA data preparation
setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")
mirna_data_raw<-read.csv("mirna_seq.OV-AU.csv", h = T)
length(unique(mirna_data_raw$icgc_donor_id))
table_ids_2<-matrix(table(mirna_data_raw$icgc_donor_id))
unique(table_ids_2[,1])

colnames(mirna_data_raw)
var_names_run_2<-c("icgc_donor_id", "icgc_specimen_id",
                 "icgc_sample_id", "submitted_sample_id",
                 "analysis_id", "mirna_id",
                 "normalized_read_count")

donor_data_2<-divide_donor_data(mirna_data_raw, var_names_run_2)

rep_names_2<-get_rep_sample_names(donor_data_2)
all.equal(names(rep_names_2[["specimen"]]), names(rep_names_2[["sample"]]))
rep_specimen_2<-rep_names_2[["specimen"]]

get_corr_rep_samples_2<-function(donor_list, rep_name_list){
  var_names<-c("mirna_id", "normalized_read_count")
  corr_cache<-list()
  for(i in names(rep_name_list)){
    print(i)
    temp_corr<-c()
    rep_name_1<-rep_name_list[[i]][1]
    temp_data_1<-donor_list[[i]][grep(rep_name_1, donor_list[[i]]$icgc_specimen_id), var_names]
    colnames(temp_data_1)[2]<-"normalized_read_count_1"
    for(j in 2:length(rep_name_list[[i]])){
      rep_name_c<-rep_name_list[[i]][j]
      temp_data_c<-donor_list[[i]][grep(rep_name_c, donor_list[[i]]$icgc_specimen_id), var_names]
      colnames(temp_data_c)[2]<-"normalized_read_count_c"
      temp_data_comb<-merge(temp_data_1, temp_data_c, by = "mirna_id")
      if(anyNA(temp_data_comb)){
        print("Warning! Merging did not go well.")
      }
      temp_corr<-c(temp_corr, 
                   cor(temp_data_comb$normalized_read_count_1, 
                       temp_data_comb$normalized_read_count_c,
                       method = "spearman"))
    }
    corr_cache[[i]]<-temp_corr
  }
  return(corr_cache)
}

get_corr_test_rep_samples_2<-function(donor_list, rep_name_list){
  var_names<-c("mirna_id", "normalized_read_count")
  corr_test_cache<-list()
  for(i in names(rep_name_list)){
    print(i)
    temp_corr_test<-c()
    rep_name_1<-rep_name_list[[i]][1]
    temp_data_1<-donor_list[[i]][grep(rep_name_1, donor_list[[i]]$icgc_specimen_id), var_names]
    colnames(temp_data_1)[2]<-"normalized_read_count_1"
    for(j in 2:length(rep_name_list[[i]])){
      rep_name_c<-rep_name_list[[i]][j]
      temp_data_c<-donor_list[[i]][grep(rep_name_c, donor_list[[i]]$icgc_specimen_id), var_names]
      colnames(temp_data_c)[2]<-"normalized_read_count_c"
      temp_data_comb<-merge(temp_data_1, temp_data_c, by = "mirna_id")
      if(anyNA(temp_data_comb)){
        print("Warning! Merging did not go well.")
      }
      temp_corr_test<-c(temp_corr_test, 
                        cor.test(temp_data_comb$normalized_read_count_1, 
                                 temp_data_comb$normalized_read_count_c,
                                 method = "kendall")$p.value)
    }
    corr_test_cache[[i]]<-temp_corr_test
  }
  return(corr_test_cache)
}

corr_cache_run_2<-get_corr_rep_samples_2(donor_data_2, rep_specimen_2)

#corr_test_cache_run_2<-get_corr_test_rep_samples_2(donor_data_2, rep_specimen_2)

reform_exp_data_2<-function(donor_list, rep_name_list, corr_cache, path){
  setwd(path)
  specimen_ids<-c()
  combinedData<-c()
  var_names<-c("mirna_id", "normalized_read_count")
  for(i in names(donor_list)){
    if(i %in% names(rep_name_list)){
      if(min(corr_cache[[i]])<0.8){
        specimen_ids<-c(specimen_ids, "mean")
        combinedData_sub<-c()
        for(j in rep_name_list[[i]]){
          temp_data_sub<-donor_list[[i]][grep(j, donor_list[[i]]$icgc_specimen_id), var_names]
          colnames(temp_data_sub)[2]<-paste("normalized_read_count", j, sep = "_")
          if(length(combinedData_sub)==0){
            combinedData_sub<-temp_data_sub
          }else{
            combinedData_sub<-merge(combinedData_sub, temp_data_sub, by = "mirna_id")
          }
        }
        temp_data<-data.frame(mirna_id = combinedData_sub$mirna_id,
                              normalized_read_count = rowSums(combinedData_sub[, -1])/length(rep_name_list[[i]]))
      }else{
        specimen_ids<-c(specimen_ids, rep_name_list[[i]][1])
        temp_data<-donor_list[[i]][grep(rep_name_list[[i]][1], donor_list[[i]]$icgc_specimen_id), var_names]
      }
    }else{
      specimen_ids<-c(specimen_ids, i)
      temp_data<-donor_list[[i]][, var_names]
    }
    colnames(temp_data)[2]<-i
    if(length(combinedData)==0){
      combinedData<-temp_data
    }else{
      combinedData<-merge(combinedData, temp_data, by = "mirna_id")
    }
  }
  specimen_id_cache<-data.frame(donor = names(donor_list), specimen = specimen_ids)
  write.csv(specimen_id_cache, file = "miRNA_specimens_used.csv", row.names = F)
  return(combinedData)
}

reform_exp_data_alt_2<-function(donor_list, rep_name_list, corr_test_cache, path){
  setwd(path)
  specimen_ids<-c()
  combinedData<-c()
  var_names<-c("mirna_id", "normalized_read_count")
  for(i in names(donor_list)){
    if(i %in% names(rep_name_list)){
      if(max(corr_test_cache[[i]])>0.05){
        specimen_ids<-c(specimen_ids, "mean")
        combinedData_sub<-c()
        for(j in rep_name_list[[i]]){
          temp_data_sub<-donor_list[[i]][grep(j, donor_list[[i]]$icgc_specimen_id), var_names]
          colnames(temp_data_sub)[2]<-paste("normalized_read_count", j, sep = "_")
          if(length(combinedData_sub)==0){
            combinedData_sub<-temp_data_sub
          }else{
            combinedData_sub<-merge(combinedData_sub, temp_data_sub, by = "mirna_id")
          }
        }
        temp_data<-data.frame(mirna_id = combinedData_sub$mirna_id,
                              normalized_read_count = rowMeans(combinedData_sub[, -1]))
      }else{
        specimen_ids<-c(specimen_ids, rep_name_list[[i]][1])
        temp_data<-donor_list[[i]][grep(rep_name_list[[i]][1], donor_list[[i]]$icgc_specimen_id), var_names]
      }
    }else{
      specimen_ids<-c(specimen_ids, i)
      temp_data<-donor_list[[i]][, var_names]
    }
    colnames(temp_data)[2]<-i
    if(length(combinedData)==0){
      combinedData<-temp_data
    }else{
      combinedData<-merge(combinedData, temp_data, by = "mirna_id")
    }
  }
  specimen_id_cache<-data.frame(donor = names(donor_list), specimen = specimen_ids)
  write.csv(specimen_id_cache, file = "mirna_specimens_used.csv", row.names = F)
  return(combinedData)
}

reformed_mirna_data<-reform_exp_data_2(donor_data_2, rep_specimen_2, corr_cache_run_2, path = "D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")

#reformed_mirna_data<-reform_exp_data_alt_2(donor_data_2, rep_specimen_2, corr_test_cache_run_2, path = "D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related_cor_test_version")

length(unique(reformed_mirna_data$mirna_id))
length(unique(colnames(reformed_mirna_data)))
getwd()
write.csv(reformed_mirna_data, file = "ov_au_reformed_mirna_data.csv", row.names = F)

length(grep(".", reformed_mirna_data$mirna_id, fixed = T))

exp_matrix2<-as.matrix(reformed_mirna_data[, -1])
patient_id2<-colnames(reformed_mirna_data)[-1]
head(patient_id2)
exp_matrix2<-t(exp_matrix2)
colnames(exp_matrix2)<-reformed_mirna_data$mirna_id
rownames(exp_matrix2)<-c(1:nrow(exp_matrix2))
data_mirna_t<-data.frame(cbind(patient_id2, exp_matrix2))
colnames(data_mirna_t)[1]<-"patient_id"
all.equal(as.numeric(data_mirna_t$hsa.miR.99b.5p), as.numeric(reformed_mirna_data[grep("hsa-miR-99b-5p", reformed_mirna_data$mirna_id, fixed = T), -1]))
write.csv(data_mirna_t, file = "ov_au_mirna_data_run.csv", row.names = F)

# phenotype data processing
setwd("D:/AUTOSurv/ICGC_data")
donor_info<-read.csv("donor.OV-AU.csv", h = T)
specimen_info<-read.csv("specimen.OV-AU.csv", h = T)
sum(!is.na(specimen_info$tumour_grade))
donor_treatment_info<-specimen_info[, c("icgc_specimen_id", "icgc_donor_id", "specimen_donor_treatment_type", "specimen_donor_treatment_type_other")]
inspect_index<-which(specimen_info$specimen_donor_treatment_type=="other therapy")
unique(specimen_info$specimen_donor_treatment_type_other[inspect_index])

get_treatment_list<-function(specimen_info){
  donor_treatment_list<-list()
  disagree_list<-c()
  for(i in unique(specimen_info$icgc_donor_id)){
    donor_treatment_list[[i]]<-specimen_info[grep(i, specimen_info$icgc_donor_id),]
    if(length(unique(donor_treatment_list[[i]]$specimen_donor_treatment_type))!=1){
      if("no treatment" %in% unique(donor_treatment_list[[i]]$specimen_donor_treatment_type)){
        donor_treatment_list[[i]]<-donor_treatment_list[[i]][-grep("no treatment", donor_treatment_list[[i]]$specimen_donor_treatment_type, fixed = T),]
      }
    }
    if(length(unique(donor_treatment_list[[i]]$specimen_donor_treatment_type))!=1 | length(unique(donor_treatment_list[[i]]$specimen_donor_treatment_type_other))!=1){
     disagree_list<-c(disagree_list, i)
    }
  }
  if(length(disagree_list)==0){
    print("All agreed!")
  }else{
    print("Variation exist between specimens from the same donor.")
  }
  donor_treatment_list[["disagree_list"]]<-disagree_list
  return(donor_treatment_list)
}

donor_treatment_list<-get_treatment_list(donor_treatment_info)

reform_donor_treatment<-function(donor_list){
  donor<-c()
  treatment<-c()
  for(i in names(donor_list)){
    donor<-c(donor, i)
    if("other therapy" %in% unique(donor_list[[i]]$specimen_donor_treatment_type)){
      treatment<-c(treatment, "Surgery and Chemotherapy")
    }else if("surgery" %in% unique(donor_list[[i]]$specimen_donor_treatment_type)){
      treatment<-c(treatment, "surgery")
    }else{
      treatment<-c(treatment, "no treatment")
    }
  }
  donor_treatment<-data.frame(donor = donor,
                              treatment = treatment)
  return(donor_treatment)
}

var_names_run_3<-c("icgc_specimen_id", "icgc_donor_id", "specimen_donor_treatment_type", "specimen_donor_treatment_type_other")
donor_data_3<-divide_donor_data(specimen_info, var_names_run_3)
donor_treatment<-reform_donor_treatment(donor_data_3)
dummy_treatment_info<-data.frame(donor = donor_treatment$donor,
                                 treatment = ifelse(donor_treatment$treatment=="surgery" | donor_treatment$treatment=="Surgery and Chemotherapy", 1, 0))

donor_var_run<-c("icgc_donor_id", "donor_vital_status", "donor_age_at_diagnosis", 
                 "donor_diagnosis_icd10", "donor_tumour_stage_at_diagnosis", "donor_survival_time")
selected_donor_info<-donor_info[, donor_var_run]

survival_info<-data.frame(donor = selected_donor_info$icgc_donor_id,
                          OS = ifelse(selected_donor_info$donor_vital_status=="deceased", 1, 0),
                          OS.time = selected_donor_info$donor_survival_time,
                          age = selected_donor_info$donor_age_at_diagnosis,
                          icd10_c56 = ifelse(selected_donor_info$donor_diagnosis_icd10=="C56", 1, 0),
                          stage_iii = ifelse(selected_donor_info$donor_tumour_stage_at_diagnosis=="III", 1, 0))

phenotype_data<-merge(survival_info, dummy_treatment_info, by = "donor")
phenotype_data[grep("DO46586", phenotype_data$donor),]
selected_donor_info[grep("DO46586", selected_donor_info$icgc_donor_id),]
specimen_info[grep("DO46586", specimen_info$icgc_donor_id), var_names_run_3]
colnames(phenotype_data)[1]<-"patient_id"
write.csv(phenotype_data, file = "ov_au_phenotype_data_run.csv", row.names = F)

# merge data
data_exp_pathway_t
data_mirna_t
phenotype_data

merge_data_1<-merge(phenotype_data, data_exp_pathway_t, by = "patient_id")
merge_data_run<-merge(merge_data_1, data_mirna_t, by = "patient_id")
grep("hsa", colnames(merge_data_run), fixed = T)[1]
colnames(merge_data_run)[9654]
all.equal(as.numeric(merge_data_run[grep("DO46606", merge_data_run$patient_id, fixed = T), 8:9653]), 
          as.numeric(data_exp_pathway_t[grep("DO46606", data_exp_pathway_t$patient_id), -1]))
all.equal(as.numeric(merge_data_run[grep("DO46606", merge_data_run$patient_id, fixed = T), 9654:ncol(merge_data_run)]), 
          as.numeric(data_mirna_t[grep("DO46606", data_mirna_t$patient_id), -1]))
write.csv(merge_data_run, file = "ov_au_merged_data_run.csv", row.names = F)

# dividing data
rescale_zero_one<-function(data, index){
  data_mody<-data
  for(i in index:ncol(data)){
    vec_min<-min(as.numeric(data_mody[, i]))
    vec_max<-max(as.numeric(data_mody[, i]))
    data_mody[, i]<-(as.numeric(data_mody[, i])-vec_min)/(vec_max-vec_min+1e-12)
  }
  return(data_mody)
}

low_variance_index_third_pipeline<-function(data, omics_index){
  low_var_index<-c()
  for(i in omics_index:ncol(data)){
    if(var(data[,i])<0.02){
      low_var_index<-c(low_var_index, i)
    }
  }
  return(low_var_index)
}

dividing_parts_third_pipeline<-function(data, p){
  
  div<-rmultinom(nrow(data), 1, p)
  classify<-c(1,2)
  tvt<-t(div)%*%classify
  tvt<-as.numeric(tvt)
  ind_tr<-which(tvt==1,arr.ind = T)
  ind_val<-which(tvt==2,arr.ind = T)
  data_tr<-data[ind_tr, ]
  data_val<-data[ind_val, ]
  divided_data<-list(data_tr, data_val)
  names(divided_data)<-c("training", "validation")
  
  return(divided_data)
}

data_dividing_third_pipeline<-function(data, p){
  
  carry_on<-TRUE
  survey<-c("OS", "icd10_c56", "stage_iii", "treatment")
  try_count<-0
  while (carry_on) {
    probes<-0
    bad_split<-FALSE
    divided_data<-dividing_parts_third_pipeline(data, p)
    for(section in names(divided_data)){
      if(bad_split){
        break
      }else{
        for(blobs in survey){
          if(sum(divided_data[[section]][, blobs]) > (nrow(divided_data[[section]]) - 2) | sum(divided_data[[section]][, blobs]) < 2){
            bad_split<-TRUE
            break
          }else{
            probes<-probes + 1
          }
        }
      }
    }
    if(probes >= length(names(divided_data)) * length(survey)){
      carry_on<-FALSE
    }
    try_count<-try_count + 1
    if(try_count >= 200 & carry_on == TRUE){
      print("Warning: reached 200 times maximum tryout limit, division still not desirable.")
      break
    }
  }
  
  return(divided_data)
}

std_normalize_tv<-function(data_raw_tr, data_raw_val, index_age, index_omics){
  data_process_tr<-data_raw_tr
  data_process_val<-data_raw_val
  
  #normalize age
  mean_age<-mean(as.numeric(data_raw_tr[, index_age]))
  sd_age<-sd(as.numeric(data_raw_tr[, index_age]))
  
  if(sd_age==0){
    print("Warning: constant age observed!")
  }
  
  data_process_tr[, index_age]<-(as.numeric(data_process_tr[, index_age]) - mean_age)/(sd_age+1e-12)
  data_process_val[, index_age]<-(as.numeric(data_process_val[, index_age]) - mean_age)/(sd_age+1e-12)
  
  #normalize omics
  sd_omics_cache<-c()
  for(i in index_omics:ncol(data_process_tr)){
    mean_omics<-mean(as.numeric(data_raw_tr[, i]))
    sd_omics<-sd(as.numeric(data_raw_tr[, i]))
    sd_omics_cache<-c(sd_omics_cache, sd_omics)
    
    data_process_tr[, i]<-(as.numeric(data_process_tr[, i]) - mean_omics)/(sd_omics+1e-12)
    data_process_val[, i]<-(as.numeric(data_process_val[, i]) - mean_omics)/(sd_omics+1e-12)
  }
  
  if(length(which(sd_omics_cache==0, arr.ind=T))!=0){
    print(paste(length(which(sd_omics_cache==0, arr.ind=T)), " constant omics observed!", sep=""))
  }else{
    print("Relax! No constant omics observed.")
  }
  
  data_cache<-list()
  data_cache[["training"]]<-data_process_tr
  data_cache[["validation"]]<-data_process_val
  if(length(which(sd_omics_cache==0, arr.ind=T))!=0){
    data_cache[["constant_omics_index"]]<-which(sd_omics_cache==0, arr.ind=T)+(index_omics-1)
  }
  return(data_cache)
}

get_std_data_tv<-function(data_raw_tr, data_raw_val, index_age, index_omics, pathway_mask){
  data_std<-std_normalize_tv(data_raw_tr, data_raw_val, index_age, index_omics)
  data_train_std<-data_std[["training"]]
  data_valid_std<-data_std[["validation"]]
  gene_start_std<-min(grep("ENSG", colnames(data_train_std), fixed = T))
  gene_end_std<-max(grep("ENSG", colnames(data_train_std), fixed = T))
  mirna_start_std<-min(grep("hsa", colnames(data_train_std), fixed = T))
  mirna_end_std<-max(grep("hsa", colnames(data_train_std), fixed = T))
  data_cache<-list()
  data_cache[["data_train_std"]]<-data_train_std
  data_cache[["data_train_gene_std"]]<-data_train_std[,1:gene_end_std]
  data_cache[["data_train_mirna_std"]]<-data_train_std[,c(1:(gene_start_std-1),mirna_start_std:mirna_end_std)]
  data_cache[["data_valid_std"]]<-data_valid_std
  data_cache[["data_valid_gene_std"]]<-data_valid_std[,1:gene_end_std]
  data_cache[["data_valid_mirna_std"]]<-data_valid_std[,c(1:(gene_start_std-1),mirna_start_std:mirna_end_std)]
  data_cache[["pathway_mask_std"]]<-pathway_mask
  if("constant_omics_index" %in% names(data_std)){
    data_cache[["exist_constant_omics"]]<-TRUE
  }else{
    data_cache[["exist_constant_omics"]]<-FALSE
  }
  return(data_cache)
}

third_pipeline_data_division<-function(out_dir, data_raw, pathway_mask, age_index = 4, omics_index = 8){
  print("Train/test split.")
  first_div<-data_dividing_third_pipeline(data_raw, c(0.8, 0.2))
  data_train_overall<-first_div[["training"]]
  data_test_overall<-first_div[["validation"]]
  
  data_train_minmax_overall<-rescale_zero_one(data_train_overall, omics_index)
  data_test_minmax_overall<-rescale_zero_one(data_test_overall, omics_index)
  
  low_var_index<-low_variance_index_third_pipeline(data_train_minmax_overall, omics_index)
  data_train_overall<-data_train_overall[,-low_var_index]
  data_test_overall<-data_test_overall[,-low_var_index]
  data_train_minmax_overall<-data_train_minmax_overall[,-low_var_index]
  data_test_minmax_overall<-data_test_minmax_overall[,-low_var_index]
  
  gene_start<-min(grep("ENSG", colnames(data_train_overall), fixed = T))
  gene_end<-max(grep("ENSG", colnames(data_train_overall), fixed = T))
  mirna_start<-min(grep("hsa", colnames(data_train_overall), fixed = T))
  mirna_end<-max(grep("hsa", colnames(data_train_overall), fixed = T))
  
  pathway_mask_run<-pathway_mask[match(colnames(data_train_overall)[gene_start:gene_end], rownames(pathway_mask)),]
  ns_pathway_index<-which(colSums(pathway_mask_run)==0, arr.ind=T)
  if(length(ns_pathway_index)!=0){
    pathway_mask_run<-pathway_mask_run[,-ns_pathway_index]
  }
  
  data_std_overall_cache<-get_std_data_tv(data_train_overall, data_test_overall, age_index, omics_index, pathway_mask_run)
  
  setwd(out_dir)
  if(!("train_test_split" %in% list.files())){
    dir.create("train_test_split")
  }
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("no_normalization" %in% list.files())){
    dir.create("no_normalization")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "no_normalization", sep = "/"))
  write.csv(data_train_overall, file="data_train_overall.csv", row.names=F)
  write.csv(data_train_overall[,1:gene_end], file="data_train_gene_overall.csv", row.names=F)
  write.csv(data_train_overall[,c(1:(gene_start-1),mirna_start:mirna_end)], file="data_train_mirna_overall.csv", row.names=F)
  write.csv(data_test_overall, file="data_test_overall.csv", row.names=F)
  write.csv(data_test_overall[,1:gene_end], file="data_test_gene_overall.csv", row.names=F)
  write.csv(data_test_overall[,c(1:(gene_start-1),mirna_start:mirna_end)], file="data_test_mirna_overall.csv", row.names=F)
  write.csv(pathway_mask_run, file="pathway_mask.csv")
  
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("minmax_normalized" %in% list.files())){
    dir.create("minmax_normalized")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "minmax_normalized", sep = "/"))
  write.csv(data_train_minmax_overall, file="data_train_minmax_overall.csv", row.names=F)
  write.csv(data_train_minmax_overall[,1:gene_end], file="data_train_gene_minmax_overall.csv", row.names=F)
  write.csv(data_train_minmax_overall[,c(1:(gene_start-1),mirna_start:mirna_end)], file="data_train_mirna_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall, file="data_test_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall[,1:gene_end], file="data_test_gene_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall[,c(1:(gene_start-1),mirna_start:mirna_end)], file="data_test_mirna_minmax_overall.csv", row.names=F)
  write.csv(pathway_mask_run, file="pathway_mask.csv")
  
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("std_normalized" %in% list.files())){
    dir.create("std_normalized")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "std_normalized", sep = "/"))
  for(i in names(data_std_overall_cache)){
    if(i != "pathway_mask_std" & i != "exist_constant_omics"){
      write.csv(data_std_overall_cache[[i]], file=paste(paste(i,"overall",sep="_"),".csv",sep=""), row.names=F)
    }else if(i != "exist_constant_omics"){
      write.csv(data_std_overall_cache[[i]], file=paste(paste(i,"overall",sep="_"),".csv",sep=""))
    }
  }
  
  file_names<-paste("tune", seq(1:5), sep = "_")
  for(folder in file_names){
    print(paste(folder, " split.", sep=""))
    setwd(out_dir)
    if(!(folder %in% list.files())){
      dir.create(folder)
    }
    
    carry_on_tv<-TRUE
    try_count_tv<-0
    while(carry_on_tv){
      tune_div<-data_dividing_third_pipeline(data_train_overall, c(0.8, 0.2))
      data_train_tune<-tune_div[["training"]]
      data_valid_tune<-tune_div[["validation"]]
      data_std_tune_cache<-get_std_data_tv(data_train_tune, data_valid_tune, age_index, omics_index, pathway_mask_run)
      if(data_std_tune_cache[["exist_constant_omics"]] == FALSE){
        carry_on_tv<-FALSE
      }
      try_count_tv<-try_count_tv+1
      if(try_count_tv >= 50 & carry_on_tv == TRUE){
        print("Warning: reached 50 times maximum tryout limit, division still not desirable.")
        break
      }
    }
    
    data_train_minmax_tune<-rescale_zero_one(data_train_tune, omics_index)
    data_valid_minmax_tune<-rescale_zero_one(data_valid_tune, omics_index)
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("no_normalization" %in% list.files())){
      dir.create("no_normalization")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "no_normalization", sep = "/"))
    write.csv(data_train_tune, file=paste(paste("data_train", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_tune[,1:gene_end], file=paste(paste("data_train_gene", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_tune[,c(1:(gene_start-1),mirna_start:mirna_end)], file=paste(paste("data_train_mirna", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune, file=paste(paste("data_valid", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune[,1:gene_end], file=paste(paste("data_valid_gene", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune[,c(1:(gene_start-1),mirna_start:mirna_end)], file=paste(paste("data_valid_mirna", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(pathway_mask_run, file="pathway_mask.csv")
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("minmax_normalized" %in% list.files())){
      dir.create("minmax_normalized")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "minmax_normalized", sep = "/"))
    write.csv(data_train_minmax_tune, file=paste(paste("data_train_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_minmax_tune[,1:gene_end], file=paste(paste("data_train_gene_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_minmax_tune[,c(1:(gene_start-1),mirna_start:mirna_end)], file=paste(paste("data_train_mirna_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune, file=paste(paste("data_valid_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune[,1:gene_end], file=paste(paste("data_valid_gene_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune[,c(1:(gene_start-1),mirna_start:mirna_end)], file=paste(paste("data_valid_mirna_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(pathway_mask_run, file="pathway_mask.csv")
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("std_normalized" %in% list.files())){
      dir.create("std_normalized")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "std_normalized", sep = "/"))
    for(i in names(data_std_tune_cache)){
      if(i != "pathway_mask_std" & i != "exist_constant_omics"){
        write.csv(data_std_tune_cache[[i]], file=paste(paste(i,folder,sep="_"),".csv",sep=""), row.names=F)
      }else if(i != "exist_constant_omics"){
        write.csv(data_std_tune_cache[[i]], file=paste(paste(i,folder,sep="_"),".csv",sep=""))
      }
    }
  }
}

#out_dir<-"D:/AUTOSurv/ICGC_data/divided_data"
out_dir<-"D:/AUTOSurv/ICGC_data/OV_AU/divided_data_new/div_5"
all.equal(colnames(merge_data_run)[8:9653],rownames(reactome_pathway_mask_shrinked))
third_pipeline_data_division(out_dir, merge_data_run, reactome_pathway_mask_shrinked, 4, 8)

# Uni-variate Cox-PH analysis data preparation
setwd("D:/AUTOSurv/ICGC_data/OV_AU/uni_coxph_analysis")
age_os_data<-merge_data_run[,c("OS", "OS.time", "age")]
write.csv(age_os_data, file = "age_os_data.csv", row.names = F)

c56_os_data<-merge_data_run[,c("OS", "OS.time", "icd10_c56")]
write.csv(c56_os_data, file = "c56_os_data.csv", row.names = F)

stage_os_data<-merge_data_run[,c("OS", "OS.time", "stage_iii")]
write.csv(stage_os_data, file = "stage_os_data.csv", row.names = F)

treatment_os_data<-merge_data_run[,c("OS", "OS.time", "treatment")]
write.csv(treatment_os_data, file = "treatment_os_data.csv", row.names = F)

# code testing
x<-c(1,2,3,4,5)
y<-c(1,1,5)
z<-c()
which(x %in% y)
which(table(pathway_data_merged$ID)>1, arr.ind = T)
test_data_1<-data.frame(donor = c("D1", "D2", "D3", "D4", "D5"),
                        X1 = c(1, 1, 1, 1, 1),
                        X2 = c(2, 2, 2, 2, 2),
                        X3 = c(3, 3, 3, 3, 3))
X_sum<-rowSums(test_data_1[, -1])
rowMeans(test_data_1[, -1])
get_exp_mean<-function(donor_data, specimens, feature){
  combinedData<-c()
  var_names<-c(feature, "normalized_read_count")
  for(i in specimens){
    temp_data<-donor_data[grep(i, donor_data$icgc_specimen_id, fixed = T), var_names]
    colnames(temp_data)[2]<-i
    if(length(combinedData)==0){
      combinedData<-temp_data
    }else{
      combinedData<-merge(combinedData, temp_data, by = feature)
    }
  }
  combinedData_mean<-data.frame(feature = combinedData[, 1],
                                normalized_read_count = rowMeans(combinedData[, -1]))
  colnames(combinedData_mean)[1]<-feature
  return(combinedData_mean)
}
test_exp_means<-get_exp_mean(donor_data_2[["DO46329"]], rep_specimen_2[["DO46329"]], "mirna_id")
all.equal(test_exp_means$normalized_read_count, reformed_mirna_data$DO46329)
setwd("D:/AUTOSurv/ICGC_data/OV_AU/data_processing_related")
pathway_mask_check<-read.table("pathway_mask_general.csv", header = T, row.names = 1, sep = ",")
for(i in 1:ncol(pathway_mask_check)){
  colnames(pathway_mask_check)[i]<-paste(unlist(strsplit(colnames(pathway_mask_check)[i], split = ".", fixed = T)), collapse = "-")
}
all.equal(colnames(pathway_mask_check), colnames(reactome_pathway_mask))

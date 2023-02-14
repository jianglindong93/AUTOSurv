# load gene expression data that were obtained from UCSC Xena data portal, example shown here is from the TCGA-BRCA dataset
data_gene_exp<-read.csv("TCGA-BRCA.htseq_fpkm.csv", h = T)

ensembl_list<-unlist(strsplit(data_gene_exp$Ensembl_ID, split = ".", fixed = T))[seq(1,2*length(data_gene_exp$Ensembl_ID),2)]
length(grep("ENSG", ensembl_list))
length(unique(ensembl_list))
data_gene_exp$Ensembl_ID<-ensembl_list

divide_index<-seq(1, nrow(data_gene_exp), 6000)
gene_id_divided<-list()
for(i in 1:length(divide_index)){
  batch_name<-paste("gene_id_batch", i, sep = "_")
  if(i<length(divide_index)){
    gene_id_divided[[batch_name]]<-ensembl_list[c(divide_index[i]:(divide_index[i+1]-1))]
  }
  else{
    gene_id_divided[[batch_name]]<-ensembl_list[c(divide_index[i]:nrow(data_gene_exp))]
  }
}
batch_name_cache<-paste("gene_id_batch", seq(1,length(gene_id_divided),1), sep = "_")
for(i in 1:length(gene_id_divided)){
  write.csv(gene_id_divided[[batch_name_cache[i]]], file = paste(batch_name_cache[i], ".csv", sep = ""), row.names = F)
}

anyNA(data_gene_exp)

# load Reactome pathway information that is obtained from DAVID
merge_cache<-function(data_cache){
  merged_data<-c()
  for(i in names(data_cache)){
    merged_data<-rbind(merged_data, as.data.frame(data_cache[[i]]))
  }
  return(merged_data)
}

pathway_batch_cache<-list()
pathway_batch_names<-paste("path_batch", seq(1,length(gene_id_divided),1), sep = "_")
for(i in pathway_batch_names){
  pathway_batch_cache[[i]]<-read.csv(paste(i, ".csv", sep = ""), h = T)
}
pathway_data_merged<-merge_cache(pathway_batch_cache)

ind_split<-grep(",",pathway_data_merged$ID,fixed = T)
ind_split

pathway_data_merged$ID[2617]<-"ENSG00000174483"
pathway_data_merged$ID[2680]<-"ENSG00000100890"
pathway_data_merged$ID[7522]<-"ENSG00000107077"

ind_split_test<-grep(",",pathway_data_merged$ID,fixed = T)
length(ind_split_test)

get_path_mask_r<-function(path_data){
  path_list<-strsplit(path_data[, 4], split = ",")
  path_row<-list()
  path_row_temp<-c()
  for(i in 1:length(path_list)){
    for(j in 1:length(path_list[[i]])){
      path_unit<-strsplit(path_list[[i]][j], split = ":")[[1]][1]
      path_row_temp<-c(path_row_temp, path_unit)
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

reactome_pathway_mask<-get_path_mask_r(pathway_data_merged)
removable_pathway_index<-filter_pathway_index(reactome_pathway_mask)
length(removable_pathway_index)
reactome_pathway_mask<-reactome_pathway_mask[,-removable_pathway_index]
test_col_sum<-colSums(reactome_pathway_mask)
min(test_col_sum)
max(test_col_sum)

row_sum_genes<-rowSums(reactome_pathway_mask)
removable_gene_index<-which(row_sum_genes==0,arr.ind = T)
#sum(reactome_pathway_mask[removable_gene_index[1000],])
reactome_pathway_mask<-reactome_pathway_mask[-removable_gene_index,]

test_index_gene<-rowSums(reactome_pathway_mask)
min(test_index_gene)

data_gene_pathway<-data_gene_exp[match(rownames(reactome_pathway_mask),data_gene_exp$Ensembl_ID),]
all.equal(data_gene_pathway$Ensembl_ID, rownames(reactome_pathway_mask))
#write.csv(data_gene_pathway, file = "temp_data_gene_exp.csv", row.names = F)

#BiocManager::install("biomaRt")
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes_to_chromosomes_2nd<-getBM(filters = "ensembl_gene_id",
                                attributes = c("ensembl_gene_id", "chromosome_name"),
                                values = data_gene_pathway$Ensembl_ID,
                                mart = mart)
head(genes_to_chromosomes_2nd)
length(unique(genes_to_chromosomes_2nd$chromosome_name))
removable_chro<-c("Y", "MT")
removable_chro_index<-which(genes_to_chromosomes_2nd$chromosome_name %in% removable_chro, arr.ind = T)
genes_to_chro_shrinked<-genes_to_chromosomes_2nd[-removable_chro_index,]
length(unique(genes_to_chro_shrinked$chromosome_name))
unique(genes_to_chro_shrinked$chromosome_name)

reactome_pathway_mask_shrinked<-reactome_pathway_mask[match(genes_to_chro_shrinked$ensembl_gene_id,rownames(reactome_pathway_mask)),]
all.equal(rownames(reactome_pathway_mask_shrinked), genes_to_chro_shrinked$ensembl_gene_id)
removable_pathway_index_2<-which(colSums(reactome_pathway_mask_shrinked)==0, arr.ind = T)
length(removable_pathway_index_2)
min(colSums(reactome_pathway_mask_shrinked))
write.csv(reactome_pathway_mask_shrinked, file = "pathway_mask_general.csv", row.names = F)

data_gene_pathway_shrinked<-data_gene_pathway[match(rownames(reactome_pathway_mask_shrinked), data_gene_pathway$Ensembl_ID),]
all.equal(data_gene_pathway_shrinked$Ensembl_ID,rownames(reactome_pathway_mask_shrinked))

exp_matrix<-as.matrix(data_gene_pathway_shrinked[,-1])
patient_id<-colnames(data_gene_pathway_shrinked)[-1]
for(i in 1:length(patient_id)){
  patient_id[i]<-paste(unlist(strsplit(patient_id[i],split = ".",fixed = T)),collapse = "-")
}
head(patient_id)
exp_matrix<-t(exp_matrix)
colnames(exp_matrix)<-data_gene_pathway_shrinked$Ensembl_ID
rownames(exp_matrix)<-c(1:nrow(exp_matrix))
data_exp_pathway_t<-data.frame(cbind(patient_id,exp_matrix))
data_exp_pathway_t$patient_id<-patient_id

#miRNA data
data_mirna_raw<-read.csv("TCGA-BRCA.mirna.csv", h = T)

mirna_matrix<-as.matrix(data_mirna_raw[,-1])
patient_id_mirna<-colnames(data_mirna_raw)[-1]
for(i in 1:length(patient_id_mirna)){
  patient_id_mirna[i]<-paste(unlist(strsplit(patient_id_mirna[i],split = ".",fixed = T)),collapse = "-")
}
head(patient_id_mirna)
mirna_matrix<-t(mirna_matrix)
colnames(mirna_matrix)<-data_mirna_raw$miRNA_ID
rownames(mirna_matrix)<-c(1:nrow(mirna_matrix))
data_mirna_t<-data.frame(cbind(patient_id_mirna, mirna_matrix))
colnames(data_mirna_t)[1]<-"patient_id"
#data_mirna_t$patient_id<-patient_id_mirna

head(as.numeric(data_mirna_t[,96]))
head(as.numeric(data_mirna_raw[95,-1]))

#phenotype data
data_pheno<-read.csv("TCGA-BRCA.GDC_phenotype.csv", h = T)
colnames(data_pheno)[grep("stage",colnames(data_pheno))]
colnames(data_pheno)[grep("age",colnames(data_pheno))]

disease_stage<-data_pheno$tumor_stage.diagnoses
patient_age<-data_pheno$age_at_initial_pathologic_diagnosis
primary_site<-data_pheno$primary_site
gender<-data_pheno$gender.demographic
race<-data_pheno$race.demographic

length(unique(primary_site))
unique(primary_site)

data_pheno_extract<-data.frame(patient_id = data_pheno$submitter_id.samples,
                               age = patient_age,
                               gender = gender,
                               race = race,
                               disease_stage = disease_stage)

anyNA(data_pheno_extract)
anyNA(data_pheno_extract$patient_id)
na_index<-c()
for(i in 2:(ncol(data_pheno_extract))){
  temp_index<-which(is.na(data_pheno_extract[,i]),arr.ind = T)
  na_index<-c(na_index, temp_index)
}
na_index<-unique(na_index)
data_pheno_extract<-data_pheno_extract[-na_index,]
anyNA(data_pheno_extract)

male_index<-which(data_pheno_extract$gender=="male",arr.ind = T)
length(male_index)

data_pheno_extract_f<-data_pheno_extract[-male_index,]
unique(data_pheno_extract_f$gender)

race_removable_index<-which(data_pheno_extract_f$race=="not reported",arr.ind = T)
stage_removable_index<-which(data_pheno_extract_f$disease_stage=="not reported"|data_pheno_extract_f$disease_stage=="stage x",arr.ind = T)
removable_index_pheno<-union(stage_removable_index, race_removable_index)
data_pheno_complete<-data_pheno_extract_f[-removable_index_pheno,]
unique(data_pheno_complete$race)
unique(data_pheno_complete$disease_stage)
unique(data_pheno_complete$gender)

stage_i_cache<-c("stage i", "stage ia", "stage ib")
stage_i<-ifelse(data_pheno_complete$disease_stage%in%stage_i_cache,1,0)
stage_ii_cache<-c("stage ii", "stage iia", "stage iib")
stage_ii<-ifelse(data_pheno_complete$disease_stage%in%stage_ii_cache,1,0)
race_white<-ifelse(data_pheno_complete$race=="white",1,0)

data_pheno_run<-data.frame(patient_id = data_pheno_complete$patient_id,
                           age = data_pheno_complete$age,
                           race_white = race_white,
                           stage_i = stage_i,
                           stage_ii = stage_ii)
table(data_pheno_run$race_white)
table(data_pheno_run$stage_i)
table(data_pheno_run$stage_ii)

#survival data
data_survival_pre<-read.csv("TCGA-BRCA.survival.csv", h = T)
anyNA(data_survival_pre)
removable_os_index<-which(is.na(data_survival_pre$OS),arr.ind = T)
length(removable_os_index)
removable_ostime_index<-which(is.na(data_survival_pre$OS.time),arr.ind = T)
length(removable_ostime_index)
#removable_survival_index<-union(removable_os_index,removable_ostime_index)

data_survival_run<-data_survival_pre[,c("sample", "OS", "OS.time")]
colnames(data_survival_run)[1]<-"patient_id"

#merge data
data_merge_1<-merge(data_survival_run, data_pheno_run, by = "patient_id", all=T)
data_merge_2_expand<-merge(data_merge_1, data_mirna_t, by = "patient_id", all=T)
data_merge_expand<-merge(data_merge_2_expand, data_exp_pathway_t, by = "patient_id", all=T)

na_index_merge_expand<-c()
for(i in 1:nrow(data_merge_expand)){
  if(anyNA(data_merge_expand[i,])){
    na_index_merge_expand<-c(na_index_merge_expand, i)
  }
}
length(na_index_merge_expand)
data_merge_expand<-data_merge_expand[-na_index_merge_expand,]
anyNA(data_merge_expand)

colnames(data_merge_expand)[1889]
head(as.numeric(data_merge_expand[500,1889:ncol(data_merge_expand)]))
head(as.numeric(data_exp_pathway_t[grep(data_merge_expand$patient_id[500], data_exp_pathway_t$patient_id),-1]))
all.equal(as.numeric(data_merge_expand[1,1889:ncol(data_merge_expand)]), as.numeric(data_exp_pathway_t[grep(data_merge_expand$patient_id[1], data_exp_pathway_t$patient_id),-1]))

head(as.numeric(data_exp_pathway_t[,501]))
head(as.numeric(data_gene_pathway_shrinked[500,-1]))
all.equal(as.numeric(data_exp_pathway_t[,2]), as.numeric(data_gene_pathway_shrinked[1,-1]))

head(as.numeric(data_mirna_t[,201]))
head(as.numeric(data_mirna_raw[200,-1]))
all.equal(as.numeric(data_mirna_t[,2]), as.numeric(data_mirna_raw[1,-1]))

head(as.numeric(data_merge_expand[500,8:1888]))
head(as.numeric(data_mirna_t[grep(data_merge_expand$patient_id[500], data_mirna_t$patient_id),-1]))
all.equal(as.numeric(data_merge_expand[500,8:1888]), as.numeric(data_mirna_t[grep(data_merge_expand$patient_id[500], data_mirna_t$patient_id),-1]))

write.csv(data_merge_expand, file = "data_merged_expand.csv", row.names = F)

# data division
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
  survey<-c("OS", "race_white", "stage_i", "stage_ii")
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
          if(sum(divided_data[[section]][, blobs]) > (nrow(divided_data[[section]]) - 10) | sum(divided_data[[section]][, blobs]) < 10){
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
  if("constant_omics_index" %in% names(data_std)){
    removable_omics<-data_std[["constant_omics_index"]]
    data_train_std<-data_train_std[,-removable_omics]
    data_valid_std<-data_valid_std[,-removable_omics]
  }
  gene_start_std<-min(grep("ENSG", colnames(data_train_std)))
  gene_end_std<-max(grep("ENSG", colnames(data_train_std)))
  mirna_start_std<-min(grep("hsa", colnames(data_train_std)))
  mirna_end_std<-max(grep("hsa", colnames(data_train_std)))
  pathway_mask_std<-pathway_mask[match(colnames(data_train_std)[gene_start_std:gene_end_std], rownames(pathway_mask)),]
  ns_pathway_std_index<-which(colSums(pathway_mask_std)==0, arr.ind=T)
  if(length(ns_pathway_std_index)!=0){
    pathway_mask_std<-pathway_mask_std[,-ns_pathway_std_index]
  }
  data_cache<-list()
  data_cache[["data_train_std"]]<-data_train_std
  data_cache[["data_train_mirna_std"]]<-data_train_std[,1:mirna_end_std]
  data_cache[["data_train_gene_std"]]<-data_train_std[,c(1:(mirna_start_std-1),gene_start_std:gene_end_std)]
  data_cache[["data_valid_std"]]<-data_valid_std
  data_cache[["data_valid_mirna_std"]]<-data_valid_std[,1:mirna_end_std]
  data_cache[["data_valid_gene_std"]]<-data_valid_std[,c(1:(mirna_start_std-1),gene_start_std:gene_end_std)]
  data_cache[["pathway_mask_std"]]<-pathway_mask_std
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
  
  gene_start<-min(grep("ENSG", colnames(data_train_overall)))
  gene_end<-max(grep("ENSG", colnames(data_train_overall)))
  mirna_start<-min(grep("hsa", colnames(data_train_overall)))
  mirna_end<-max(grep("hsa", colnames(data_train_overall)))
  
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
  write.csv(data_train_overall[,1:mirna_end], file="data_train_mirna_overall.csv", row.names=F)
  write.csv(data_train_overall[,c(1:(mirna_start-1),gene_start:gene_end)], file="data_train_gene_overall.csv", row.names=F)
  write.csv(data_test_overall, file="data_test_overall.csv", row.names=F)
  write.csv(data_test_overall[,1:mirna_end], file="data_test_mirna_overall.csv", row.names=F)
  write.csv(data_test_overall[,c(1:(mirna_start-1),gene_start:gene_end)], file="data_test_gene_overall.csv", row.names=F)
  write.csv(pathway_mask_run, file="pathway_mask.csv")
  
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("minmax_normalized" %in% list.files())){
    dir.create("minmax_normalized")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "minmax_normalized", sep = "/"))
  write.csv(data_train_minmax_overall, file="data_train_minmax_overall.csv", row.names=F)
  write.csv(data_train_minmax_overall[,1:mirna_end], file="data_train_mirna_minmax_overall.csv", row.names=F)
  write.csv(data_train_minmax_overall[,c(1:(mirna_start-1),gene_start:gene_end)], file="data_train_gene_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall, file="data_test_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall[,1:mirna_end], file="data_test_mirna_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall[,c(1:(mirna_start-1),gene_start:gene_end)], file="data_test_gene_minmax_overall.csv", row.names=F)
  write.csv(pathway_mask_run, file="pathway_mask.csv")
  
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("std_normalized" %in% list.files())){
    dir.create("std_normalized")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "std_normalized", sep = "/"))
  for(i in names(data_std_overall_cache)){
    write.csv(data_std_overall_cache[[i]], file=paste(paste(i,"overall",sep="_"),".csv",sep=""), row.names=F)
  }
  
  file_names<-paste("tune", seq(1:10), sep = "_")
  for(folder in file_names){
    print(paste(folder, " split.", sep=""))
    setwd(out_dir)
    if(!(folder %in% list.files())){
      dir.create(folder)
    }
    
    tune_div<-data_dividing_third_pipeline(data_train_overall, c(0.8, 0.2))
    data_train_tune<-tune_div[["training"]]
    data_valid_tune<-tune_div[["validation"]]
    
    data_train_minmax_tune<-rescale_zero_one(data_train_tune, omics_index)
    data_valid_minmax_tune<-rescale_zero_one(data_valid_tune, omics_index)
    
    data_std_tune_cache<-get_std_data_tv(data_train_tune, data_valid_tune, age_index, omics_index, pathway_mask_run)
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("no_normalization" %in% list.files())){
      dir.create("no_normalization")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "no_normalization", sep = "/"))
    write.csv(data_train_tune, file=paste(paste("data_train", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_tune[,1:mirna_end], file=paste(paste("data_train_mirna", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_tune[,c(1:(mirna_start-1),gene_start:gene_end)], file=paste(paste("data_train_gene", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune, file=paste(paste("data_valid", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune[,1:mirna_end], file=paste(paste("data_valid_mirna", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_tune[,c(1:(mirna_start-1),gene_start:gene_end)], file=paste(paste("data_valid_gene", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(pathway_mask_run, file="pathway_mask.csv")
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("minmax_normalized" %in% list.files())){
      dir.create("minmax_normalized")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "minmax_normalized", sep = "/"))
    write.csv(data_train_minmax_tune, file=paste(paste("data_train_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_minmax_tune[,1:mirna_end], file=paste(paste("data_train_mirna_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_train_minmax_tune[,c(1:(mirna_start-1),gene_start:gene_end)], file=paste(paste("data_train_gene_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune, file=paste(paste("data_valid_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune[,1:mirna_end], file=paste(paste("data_valid_mirna_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune[,c(1:(mirna_start-1),gene_start:gene_end)], file=paste(paste("data_valid_gene_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(pathway_mask_run, file="pathway_mask.csv")
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("std_normalized" %in% list.files())){
      dir.create("std_normalized")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "std_normalized", sep = "/"))
    for(i in names(data_std_tune_cache)){
      write.csv(data_std_tune_cache[[i]], file=paste(paste(i,folder,sep="_"),".csv",sep=""), row.names=F)
    }
  }
}

out_dir<-"D:/DL/TCGA data/TCGA-GDC-BRCA/3rd_pipeline"
all.equal(colnames(data_merge_expand)[1889:ncol(data_merge_expand)],rownames(reactome_pathway_mask_shrinked))
third_pipeline_data_division(out_dir, data_merge_expand, reactome_pathway_mask_shrinked, 4, 8)

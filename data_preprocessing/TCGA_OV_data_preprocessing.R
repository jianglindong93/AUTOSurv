# gene expression & pathway information
data_gene_exp<-read.csv("TCGA-OV.htseq_fpkm.csv", h = T)
anyNA(data_gene_exp)
pathway_mask_general<-read.table("pathway_mask_general.csv", header = T, row.names = 1, sep = ",")

ensembl_list<-unlist(strsplit(data_gene_exp$Ensembl_ID, split = ".", fixed = T))[seq(1,2*length(data_gene_exp$Ensembl_ID),2)]
length(grep("ENSG", ensembl_list))
length(unique(ensembl_list))
data_gene_exp$Ensembl_ID<-ensembl_list

keep_id_index<-match(rownames(pathway_mask_general), data_gene_exp$Ensembl_ID)
data_gene_exp_run<-data_gene_exp[keep_id_index,]
all.equal(data_gene_exp_run$Ensembl_ID, rownames(pathway_mask_general))

exp_matrix<-as.matrix(data_gene_exp_run[,-1])
patient_id<-colnames(data_gene_exp_run)[-1]
for(i in 1:length(patient_id)){
  patient_id[i]<-paste(unlist(strsplit(patient_id[i],split = ".",fixed = T)),collapse = "-")
}
head(patient_id)
exp_matrix<-t(exp_matrix)
colnames(exp_matrix)<-data_gene_exp_run$Ensembl_ID
rownames(exp_matrix)<-c(1:nrow(exp_matrix))
data_gene_exp_run_t<-data.frame(cbind(patient_id,exp_matrix))

length(which(data_gene_exp_run_t=="", arr.ind = T))
all.equal(as.numeric(data_gene_exp_run_t$ENSG00000002016), as.numeric(data_gene_exp_run[5,-1]))
head(as.numeric(data_gene_exp_run_t$ENSG00000002016))
head(as.numeric(data_gene_exp_run[5,-1]))

# miRNA expression
data_mirna_exp<-read.csv("TCGA-OV.mirna.csv", h = T)
anyNA(data_mirna_exp)

mirna_exp_matrix<-as.matrix(data_mirna_exp[,-1])
patient_id_mirna<-colnames(data_mirna_exp)[-1]
for(i in 1:length(patient_id_mirna)){
  patient_id_mirna[i]<-paste(unlist(strsplit(patient_id_mirna[i],split = ".",fixed = T)),collapse = "-")
}
head(patient_id_mirna)
mirna_exp_matrix<-t(mirna_exp_matrix)
colnames(mirna_exp_matrix)<-data_mirna_exp$miRNA_ID
rownames(mirna_exp_matrix)<-c(1:nrow(mirna_exp_matrix))
data_mirna_exp_t<-data.frame(cbind(patient_id_mirna,mirna_exp_matrix))
colnames(data_mirna_exp_t)[1]<-"patient_id"

length(which(data_mirna_exp_t=="", arr.ind = T))
all.equal(as.numeric(data_mirna_exp_t$hsa.let.7c), as.numeric(data_mirna_exp[5,-1]))
head(as.numeric(data_mirna_exp_t$hsa.let.7c))
head(as.numeric(data_mirna_exp[5,-1]))

# phenotype
data_pheno<-read.csv("TCGA-OV.GDC_phenotype.csv", h = T)
colnames(data_pheno)
selected_pheno_var<-c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "clinical_stage", "neoplasm_histologic_grade", 
                      "gender.demographic", "race.demographic")
data_pheno_run<-data_pheno[, selected_pheno_var]
anyNA(data_pheno_run)
gender_empty_index<-which(data_pheno_run$gender.demographic=="")
unique(data_pheno_run$neoplasm_histologic_grade[gender_empty_index])

length(grep("TCGA", data_pheno_run$submitter_id.samples))
length(unique(data_pheno_run$submitter_id.samples))
nrow(data_pheno_run)

get_removable_index<-function(data){
  index_cache<-c()
  for(i in 1:ncol(data)){
    temp_index<-which(data[,i]==""|data[,i]=="GX"|is.na(data[,i])|data[,i]=="not reported")
    index_cache<-union(index_cache, temp_index)
  }
  return(index_cache)
}

pheno_removable_index<-get_removable_index(data_pheno_run)
data_pheno_run_shrinked<-data_pheno_run[-pheno_removable_index,]
unique(data_pheno_run_shrinked$age_at_initial_pathologic_diagnosis)
length(which(data_pheno_run_shrinked=="", arr.ind = T))
anyNA(data_pheno_run_shrinked)
length(grep(".", data_pheno_run_shrinked$submitter_id.samples, fixed = T))

table(data_pheno_run_shrinked$clinical_stage)

stage_h<-ifelse(data_pheno_run_shrinked$clinical_stage=="Stage IIIA"|data_pheno_run_shrinked$clinical_stage=="Stage IIIB"|data_pheno_run_shrinked$clinical_stage=="Stage IIIC"|data_pheno_run_shrinked$clinical_stage=="Stage IV", 1, 0)
grade_h<-ifelse(data_pheno_run_shrinked$neoplasm_histologic_grade=="G3"|data_pheno_run_shrinked$neoplasm_histologic_grade=="G4", 1, 0)
race_white<-ifelse(data_pheno_run_shrinked$race.demographic=="white", 1, 0)
age<-data_pheno_run_shrinked$age_at_initial_pathologic_diagnosis

data_pheno_premerge<-data.frame(cbind(data_pheno_run_shrinked$submitter_id.samples, age, stage_h, grade_h, race_white))
colnames(data_pheno_premerge)[1]<-"patient_id"

# survival
data_surv<-read.csv("TCGA-OV.survival.csv", h = T)
anyNA(data_surv)
length(which(data_surv=="", arr.ind = T))
data_surv_run<-data_surv[, c(1, 2, 4)]
colnames(data_surv_run)[1]<-"patient_id"

# merging data_gene_exp_run_t data_mirna_exp_t data_pheno_premerge data_surv_run
data_merge_1<-merge(data_surv_run, data_pheno_premerge, by = "patient_id")
data_merge_2<-merge(data_merge_1, data_mirna_exp_t, by = "patient_id")
data_merge<-merge(data_merge_2, data_gene_exp_run_t, by = "patient_id")
anyNA(data_merge)
table(data_merge$stage_h)

write.csv(data_merge, file = "data_merge.csv", row.names = F)

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
  survey<-c("OS", "race_white", "stage_h", "grade_h")
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
          if(sum(as.numeric(divided_data[[section]][, blobs])) > (nrow(divided_data[[section]]) - 5) | sum(as.numeric(divided_data[[section]][, blobs])) < 5){
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
    if(try_count >= 500 & carry_on == TRUE){
      print("Warning: reached 500 times maximum tryout limit, division still not desirable.")
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

out_dir<-"D:/DL/TCGA data/TCGA-GDC-OV/data_split"
all.equal(colnames(data_merge)[1889:ncol(data_merge)],rownames(pathway_mask_general))
third_pipeline_data_division(out_dir, data_merge, pathway_mask_general, 4, 8)
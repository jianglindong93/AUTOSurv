setwd("D:/AUTOSurv/ICGC_data/Caldas2007/data_processing_related")

# exp data
exp_data<-read.csv("naderi2007Exp_genomicMatrix.csv", h = T)
num_na_row<-0
for(i in 1:nrow(exp_data)){
  if(anyNA(exp_data[i,-1])){
    num_na_row<-num_na_row+1
  }
}
round(num_na_row/nrow(exp_data)*100, 2)

get_na_rows<-function(data_raw){
  removable_index<-c()
  num_na_entry<-c()
  for(i in 1:nrow(data_raw)){
    if(anyNA(as.numeric(data_raw[i,-1]))){
      removable_index<-c(removable_index, i)
      num_na_entry<-c(num_na_entry, sum(is.na(data_raw[i,-1])))
    }
  }
  cache_na_info<-list()
  cache_na_info[["removable_index"]]<-removable_index
  cache_na_info[["num_na_entry"]]<-num_na_entry
  return(cache_na_info)
}

anyNA(exp_data)
cache_na_info<-get_na_rows(exp_data)
rm_index1<-cache_na_info[["removable_index"]]
exp_data_compl<-exp_data[-rm_index1,]
anyNA(exp_data_compl)
min(cache_na_info[["num_na_entry"]])
max(cache_na_info[["num_na_entry"]])

# get pathway information from DAVID
setwd("D:/AUTOSurv/ICGC_data/Caldas2007/data_processing_related/gene_name_batches")
divide_index<-seq(1, nrow(exp_data_compl), 6000)
gene_id_divided<-list()
probe_list<-exp_data_compl$probe
for(i in 1:length(divide_index)){
  batch_name<-paste("gene_id_batch", i, sep = "_")
  if(i<length(divide_index)){
    gene_id_divided[[batch_name]]<-probe_list[c(divide_index[i]:(divide_index[i+1]-1))]
  }
  else{
    gene_id_divided[[batch_name]]<-probe_list[c(divide_index[i]:nrow(exp_data_compl))]
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

setwd("D:/AUTOSurv/ICGC_data/Caldas2007/data_processing_related/pathway_batches")
pathway_batch_cache<-list()
pathway_batch_names<-paste("pathway_batch", seq(1,length(gene_id_divided),1), sep = "_")
for(i in pathway_batch_names){
  pathway_batch_cache[[i]]<-read.csv(paste(i, ".csv", sep = ""), h = T)
}

pathway_data_merged<-merge_cache(pathway_batch_cache)
length(pathway_data_merged$ID)
length(unique(pathway_data_merged$ID))

which(table(pathway_data_merged$ID)>1, arr.ind = T)
length(grep("A_23_P125107", pathway_batch_cache[[3]]$ID, fixed = T))
grep("A_23_P125107", pathway_batch_cache[[3]]$ID)
pathway_batch_cache[[3]][820,]
pathway_batch_cache[[3]][821,]

# get pathway information from Reactome
write.csv(data.frame(probe = probe_list), file = "gene_list.csv", row.names = F)

# clean up pathway info
rep_id_names<-rownames(which(table(pathway_data_merged$ID)>1, arr.ind = T))
rep_freq<-c()
for(i in rep_id_names){
  rep_freq<-c(rep_freq, length(grep(i, pathway_data_merged$ID, fixed = T)))
}
max(rep_freq)
sum(rep_freq)-length(rep_id_names)

sort_rep_ids<-function(pathway_info){
  rep_id<-rownames(which(table(pathway_info$ID)>1, arr.ind = T))
  keep_ind<-c()
  check_id<-list()
  ind_cache<-list()
  rep_index<-c()
  for(i in rep_id){
    rep_index<-c(rep_index, grep(i, pathway_info$ID, fixed = T))
    temp_subpathway<-pathway_info[grep(i, pathway_info$ID, fixed = T),]
    if(length(unique(temp_subpathway[,"REACTOME_PATHWAY"]))==1){
      if(length(unique(temp_subpathway[,"ID"]))!=1){
        print(i)
        print("Inconsistent ID in keep_ind.")
      }
      keep_ind<-c(keep_ind, grep(i, pathway_info$ID, fixed = T)[1])
    }else{
      check_id[[i]]<-temp_subpathway
    }
  }
  ind_cache[["keep_ind"]]<-keep_ind
  ind_cache[["check_id"]]<-check_id
  ind_cache[["rep_index"]]<-rep_index
  return(ind_cache)
}

sorted_ids<-sort_rep_ids(pathway_data_merged)
View(sorted_ids[["check_id"]][[22]])
pathway_data_merged$ID[as.numeric(rownames(sorted_ids[["check_id"]][[4]]))]
keep_ind2<-c(3552, 460, 5716, 2920, 1199, 3534, 5773, 2917, 5687, 1082, 3169, 2938, 49, 902, 1349, 1756, 5497, 4127, 5074)
reformed_pathway_fraction<-data.frame(row_index = rownames(sorted_ids[["check_id"]][[2]]),
                                      ID = c("A_23_P118435", "A_23_P102937"),
                                      Gene.Name = c("small ubiquitin like modifier 2(SUMO2)", "small ubiquitin like modifier 3(SUMO3)"),
                                      Species = sorted_ids[["check_id"]][[2]]$Species,
                                      REACTOME_PATHWAY = sorted_ids[["check_id"]][[2]]$REACTOME_PATHWAY)
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = c(rownames(sorted_ids[["check_id"]][[3]])[2], rownames(sorted_ids[["check_id"]][[3]])[1]),
                                                                       ID = c("A_23_P13033", "A_23_P150255"),
                                                                       Gene.Name = c("RNA binding motif protein 4(RBM4)", "RNA binding motif protein 14(RBM14)"),
                                                                       Species = sorted_ids[["check_id"]][[3]]$Species,
                                                                       REACTOME_PATHWAY = c(sorted_ids[["check_id"]][[3]]$REACTOME_PATHWAY[2], sorted_ids[["check_id"]][[3]]$REACTOME_PATHWAY[1])))
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = rownames(sorted_ids[["check_id"]][[5]])[1:2],
                                                                       ID = c("A_23_P152804", "A_23_P141405"),
                                                                       Gene.Name = sorted_ids[["check_id"]][[5]]$Gene.Name[1:2],
                                                                       Species = sorted_ids[["check_id"]][[5]]$Species[1:2],
                                                                       REACTOME_PATHWAY = sorted_ids[["check_id"]][[5]]$REACTOME_PATHWAY[1:2]))
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = rownames(sorted_ids[["check_id"]][[7]]),
                                                                       ID = c("A_23_P150609", "A_23_P1981"),
                                                                       Gene.Name = sorted_ids[["check_id"]][[7]]$Gene.Name,
                                                                       Species = sorted_ids[["check_id"]][[7]]$Species,
                                                                       REACTOME_PATHWAY = sorted_ids[["check_id"]][[7]]$REACTOME_PATHWAY))
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = rownames(sorted_ids[["check_id"]][[15]]),
                                                                       ID = c("A_23_P79510", "A_23_P56709"),
                                                                       Gene.Name = sorted_ids[["check_id"]][[15]]$Gene.Name,
                                                                       Species = sorted_ids[["check_id"]][[15]]$Species,
                                                                       REACTOME_PATHWAY = sorted_ids[["check_id"]][[15]]$REACTOME_PATHWAY))
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = rownames(sorted_ids[["check_id"]][[16]])[c(1,3)],
                                                                       ID = c("A_23_P58359", "A_23_P81158"),
                                                                       Gene.Name = sorted_ids[["check_id"]][[16]]$Gene.Name[c(1,3)],
                                                                       Species = sorted_ids[["check_id"]][[16]]$Species[c(1,3)],
                                                                       REACTOME_PATHWAY = sorted_ids[["check_id"]][[16]]$REACTOME_PATHWAY[c(1,3)]))
reformed_pathway_fraction<-rbind(reformed_pathway_fraction, data.frame(row_index = rownames(sorted_ids[["check_id"]][[21]])[c(9,2,12,5,10,13,14)],
                                                                       ID = c("A_23_P93258", "A_23_P111037", "A_23_P70445", "A_23_P93282", "A_23_P133814", "A_23_P8004", "A_23_P42198"),
                                                                       Gene.Name = sorted_ids[["check_id"]][[21]]$Gene.Name[c(9,2,12,5,10,13,14)],
                                                                       Species = sorted_ids[["check_id"]][[21]]$Species[c(9,2,12,5,10,13,14)],
                                                                       REACTOME_PATHWAY = sorted_ids[["check_id"]][[21]]$REACTOME_PATHWAY[c(9,2,12,5,10,13,14)]))

fine_id<-seq(1, nrow(pathway_data_merged))[-sorted_ids[["rep_index"]]]
keep_inds<-c(sorted_ids[["keep_ind"]], keep_ind2, as.numeric(reformed_pathway_fraction$row_index))
cleaned_pathway_info<-pathway_data_merged
cleaned_pathway_info[as.numeric(reformed_pathway_fraction$row_index),]<-reformed_pathway_fraction[,-1]
keep_ids_overall<-sort(c(fine_id, keep_inds))
cleaned_pathway_info<-cleaned_pathway_info[keep_ids_overall,]
length(unique(cleaned_pathway_info$ID))
length(grep(",", cleaned_pathway_info$ID, fixed = T))
check_ids2<-cleaned_pathway_info$ID[grep(",", cleaned_pathway_info$ID, fixed = T)]
cleaned_pathway_info_run<-cleaned_pathway_info
for(i in check_ids2){
  cleaned_pathway_info_run$ID[grep(i, cleaned_pathway_info_run$ID, fixed = T)]<-unlist(strsplit(i, split = ",", fixed = T))[1]
}
length(grep(",", cleaned_pathway_info_run$ID, fixed = T))
length(grep(" ", cleaned_pathway_info_run$ID, fixed = T))

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

reactome_pathway_mask<-get_path_mask_r(cleaned_pathway_info_run)
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

data_gene_pathway<-exp_data_compl[match(rownames(reactome_pathway_mask), exp_data_compl$probe),]
all.equal(data_gene_pathway$probe, rownames(reactome_pathway_mask))
rownames(data_gene_pathway)<-seq(1,nrow(data_gene_pathway))
colnames(data_gene_pathway)[1]<-"gene_id"

setwd("D:/AUTOSurv/ICGC_data/Caldas2007/data_processing_related/gene_name_batches")
gene_names_shrinked<-data_gene_pathway$gene_id
write.csv(data.frame(gene_id = gene_names_shrinked), file = "gene_names_shrinked.csv", row.names = F)

# get chromosomal info
chromosome_info<-read.csv("chromosome_info.csv", h = T)
length(unique(chromosome_info$ID))
unique(chromosome_info$CHROMOSOME)

gene_info_match<-cleaned_pathway_info_run
gene_info_match<-gene_info_match[match(data_gene_pathway$gene_id, gene_info_match$ID),]
all.equal(gene_info_match$ID, data_gene_pathway$gene_id)
sum(gene_info_match$Gene.Name %in% chromosome_info$Gene.Name)

length(unique(chromosome_info$Gene.Name))
chromosome_info_shrinked<-chromosome_info[match(gene_info_match$Gene.Name, chromosome_info$Gene.Name),]
all.equal(chromosome_info_shrinked$Gene.Name, gene_info_match$Gene.Name)
chromosome_info_shrinked$ID<-gene_info_match$ID
multi_chromo_check<-c()
for(i in 1:nrow(chromosome_info_shrinked)){
  if(length(unlist(strsplit(chromosome_info_shrinked$CHROMOSOME[i], split = ",", fixed = T)))!=1){
    multi_chromo_check<-c(multi_chromo_check, i)
  }
}
length(multi_chromo_check)

removable_gene_chromo<-which(chromosome_info_shrinked$CHROMOSOME == "X|Y," | chromosome_info_shrinked$CHROMOSOME == "Y,", arr.ind = T)
chromosome_info_shrinked2<-chromosome_info_shrinked[-removable_gene_chromo,]

unique(chromosome_info_shrinked2$CHROMOSOME)
length(unique(chromosome_info_shrinked2$CHROMOSOME))

reactome_pathway_mask_shrinked<-reactome_pathway_mask[match(chromosome_info_shrinked2$ID, rownames(reactome_pathway_mask)),]
all.equal(rownames(reactome_pathway_mask_shrinked), chromosome_info_shrinked2$ID)
removable_pathway_index_2<-which(colSums(reactome_pathway_mask_shrinked)==0, arr.ind = T)
length(removable_pathway_index_2)
min(colSums(reactome_pathway_mask_shrinked))
setwd("D:/AUTOSurv/ICGC_data/Caldas2007")
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
all.equal(as.numeric(data_exp_pathway_t$A_23_P168909), as.numeric(data_gene_pathway_shrinked[grep("A_23_P168909", data_gene_pathway_shrinked$gene_id), -1]))
for(i in 1:nrow(data_exp_pathway_t)){
  data_exp_pathway_t$patient_id[i]<-unlist(strsplit(data_exp_pathway_t$patient_id[i], split = "X", fixed = T))[2]
}
write.csv(data_exp_pathway_t, file = "caldas2007_exp_data_run.csv", row.names = F)

# phenotype data
setwd("D:/AUTOSurv/ICGC_data/Caldas2007/data_processing_related")
patient_info<-read.csv("phenotype_data_raw.csv", h = T)

colnames(patient_info)
selected_var<-c("patient.No.", "Meno", "Size", "Grade", "Stage", "Survival", "Dead", "AGE")
selected_phenotype_info<-patient_info[, selected_var]
selected_phenotype_info<-selected_phenotype_info[-25,]
sum(selected_phenotype_info$patient.No. %in% data_exp_pathway_t$patient_id)
nrow(selected_phenotype_info)
#pick_singular<-data_exp_pathway_t$patient_id[!data_exp_pathway_t$patient_id %in% selected_phenotype_info$patient.No.]
grep("2196", selected_phenotype_info$patient.No., fixed = T)
selected_phenotype_info<-selected_phenotype_info[-62,]
anyNA(selected_phenotype_info)

phenotype_data_run<-data.frame(patient_id = selected_phenotype_info$patient.No.,
                               OS = ifelse(selected_phenotype_info$Dead == 1, 1, 0),
                               OS.time = selected_phenotype_info$Survival,
                               age = selected_phenotype_info$AGE,
                               post_meno = ifelse(selected_phenotype_info$Meno == 2, 1, 0),
                               size = selected_phenotype_info$Size,
                               stage_i = ifelse(selected_phenotype_info$Stage == 1, 1, 0),
                               grade_i = ifelse(selected_phenotype_info$Grade == 1, 1, 0),
                               grade_ii = ifelse(selected_phenotype_info$Grade == 2, 1, 0))


# merge data
merged_data<-merge(phenotype_data_run, data_exp_pathway_t, by = "patient_id")
anyNA(merged_data)
length(grep("A_23", colnames(merged_data), fixed = T))
grep("A_23", colnames(merged_data), fixed = T)[1]
colnames(merged_data)[10]
all.equal(as.numeric(merged_data[grep("2073", merged_data$patient_id, fixed = T), 10:ncol(merged_data)]), 
          as.numeric(data_exp_pathway_t[grep("2073", data_exp_pathway_t$patient_id), -1]))
setwd("D:/AUTOSurv/ICGC_data/Caldas2007")
write.csv(merged_data, file = "caldas2007_merged_data_run.csv", row.names = F)

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
  survey<-c("OS", "post_meno", "stage_i", "grade_i", "grade_ii")
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
          if(sum(divided_data[[section]][, blobs]) > (nrow(divided_data[[section]]) - 5) | sum(divided_data[[section]][, blobs]) < 5){
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
  data_cache<-list()
  data_cache[["data_train_std"]]<-data_train_std
  data_cache[["data_valid_std"]]<-data_valid_std
  data_cache[["pathway_mask_std"]]<-pathway_mask
  if("constant_omics_index" %in% names(data_std)){
    data_cache[["exist_constant_omics"]]<-TRUE
  }else{
    data_cache[["exist_constant_omics"]]<-FALSE
  }
  return(data_cache)
}

third_pipeline_data_division<-function(out_dir, data_raw, pathway_mask, age_index = 4, omics_index = 10){
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
  
  gene_start<-min(grep("A_23_", colnames(data_train_overall), fixed = T))
  gene_end<-max(grep("A_23_", colnames(data_train_overall), fixed = T))
  
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
  write.csv(data_test_overall, file="data_test_overall.csv", row.names=F)
  write.csv(pathway_mask_run, file="pathway_mask.csv")
  
  setwd(paste(out_dir, "train_test_split", sep = "/"))
  if(!("minmax_normalized" %in% list.files())){
    dir.create("minmax_normalized")
  }
  setwd(paste(paste(out_dir, "train_test_split", sep = "/"), "minmax_normalized", sep = "/"))
  write.csv(data_train_minmax_overall, file="data_train_minmax_overall.csv", row.names=F)
  write.csv(data_test_minmax_overall, file="data_test_minmax_overall.csv", row.names=F)
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
    write.csv(data_valid_tune, file=paste(paste("data_valid", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(pathway_mask_run, file="pathway_mask.csv")
    
    setwd(paste(out_dir, folder, sep = "/"))
    if(!("minmax_normalized" %in% list.files())){
      dir.create("minmax_normalized")
    }
    setwd(paste(paste(out_dir, folder, sep = "/"), "minmax_normalized", sep = "/"))
    write.csv(data_train_minmax_tune, file=paste(paste("data_train_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
    write.csv(data_valid_minmax_tune, file=paste(paste("data_valid_minmax", folder, sep="_"),".csv",sep=""), row.names=F)
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

#out_dir<-"D:/AUTOSurv/ICGC_data/Caldas2007/divided_data"
#set.seed(2)
out_dir<-"D:/AUTOSurv/ICGC_data/Caldas2007/divided_data_new/div_2"
all.equal(colnames(merged_data)[10:5708], rownames(reactome_pathway_mask_shrinked))
third_pipeline_data_division(out_dir, merged_data, reactome_pathway_mask_shrinked, 4, 10)

# Uni-variate Cox-PH analysis
setwd("D:/AUTOSurv/ICGC_data/Caldas2007/uni_coxph_analysis")
age_os_data<-merged_data[,c("OS", "OS.time", "age")]
write.csv(age_os_data, file = "age_os_data.csv", row.names = F)

postmeno_os_data<-merged_data[,c("OS", "OS.time", "post_meno")]
write.csv(postmeno_os_data, file = "postmeno_os_data.csv", row.names = F)

size_os_data<-merged_data[,c("OS", "OS.time", "size")]
write.csv(size_os_data, file = "size_os_data.csv", row.names = F)

stage_os_data<-merged_data[,c("OS", "OS.time", "stage_i")]
write.csv(stage_os_data, file = "stage_os_data.csv", row.names = F)

grade_os_data<-merged_data[,c("OS", "OS.time", "grade_i", "grade_ii")]
write.csv(grade_os_data, file = "grade_os_data.csv", row.names = F)








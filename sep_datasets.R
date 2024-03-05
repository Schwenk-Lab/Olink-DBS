# Data separation and merging of replicates

# Function for merging duplicates, triplicates and separating the datasets for datasets 1 and 2
sep_datasets <- function(npx, sinfo) {
  
  # Get only population samples, remove bridging samples in the bridged dataset
  sinfo_temp <- sinfo %>% filter((Analysis_date == "2020-06-23" & str_detect(Analysis_group, "UMAP")) |
                                 (Analysis_date == "2020-11-13" & str_detect(Analysis_group, "SPK")))
  # (before merging, store protein column order to restore later)
  prot_order <- colnames(npx)
  npx_temp <- npx[sinfo_temp$sample_id, ] %>%
    # Merge replicates with mean value by making long format data and grouping
    as.data.frame() %>% rownames_to_column("sample_id") %>%
    # Make names (individuals) to group by
    mutate(sample_id = str_remove(sample_id, "_[:alpha:]$")) %>%
    pivot_longer(cols = -sample_id, names_to = "binder", values_to = "npx") %>%
    group_by(sample_id, binder) %>%
    summarise(npx = mean(npx), .groups = "keep") %>%
    pivot_wider(id_cols = "sample_id", names_from = "binder", values_from = "npx") %>%
    column_to_rownames("sample_id") %>%
    select(!!prot_order)

  # Trim off remaints of the replicate names in the sample names
  sinfo_temp$sample_id <- str_remove(sinfo_temp$sample_id, "_[:alpha:]$")
  sinfo_temp$subj.id <- str_remove(sinfo_temp$subj.id, "_[:alpha:]$")
  
  # Separate the datasets (distinct used to remove duplicate rows that were from the replicates after removing the letters at the end of the sample names)
  sinfo_umap <- sinfo_temp %>% distinct(sample_id, .keep_all = T) %>% filter(Analysis_date == "2020-06-23") %>% as.data.frame() %>% `rownames<-`(.$sample_id)
  npx_umap <- npx_temp[sinfo_umap$sample_id, ]
  sinfo_spk <- sinfo_temp %>% distinct(sample_id, .keep_all = T) %>% filter(Analysis_date == "2020-11-13") %>% as.data.frame() %>% `rownames<-`(.$sample_id)
  npx_spk <- npx_temp[sinfo_spk$sample_id, ]
  
  return(list("npx_umap" = npx_umap, "npx_spk" = npx_spk,
              "sinfo_umap" = sinfo_umap, "sinfo_spk" = sinfo_spk))
}

# Function for merging replicates in dataset 3
merge_repl <- function(dat_list) {
  npx_out <- dat_list$npx %>%
    # Make individual-specific (not sample-specific) IDs to group by
    mutate(sample_id = str_remove(sample_id, "_B")) %>%
    # Make long for easy grouping
    pivot_longer(cols = -sample_id) %>%
    group_by(sample_id, name) %>%
    # Merge replicates by taking mean
    summarise(value = mean(value), .groups = "keep") %>%
    # Convert to wide format again
    pivot_wider(id_cols = sample_id)
  
  sinfo_out <- dat_list$sinfo %>% filter(sample_id %in% npx_out$sample_id) %>% arrange(sample_id)
  
  dat_out <- dat_list
  dat_out$npx <- npx_out
  dat_out$sinfo <- sinfo_out
  
  return(dat_out)
}

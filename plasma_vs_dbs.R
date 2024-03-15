#                   Title: Comparison of plasma and DBS eluates from 12 donors

# description: Matched dbs and plasma samples from 12 donors. Data for 92 proteins using the Olink Cardiovascular III panel
# creator: Annika Bendes
######################################################################
# load packages
library(tidyverse)
library(reshape2)
library(scales)
library(gridExtra)
require(ggrepel)
library(RColorBrewer)
require(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)


## Load gene name information:
file_name <- c("gene_name_key.csv")
seperator <- ","
gene_names_table = read.csv(file_name, sep=seperator, dec=",", skip=0, stringsAsFactors=F, header=T)

# load data with sample, protein and npx information 
input_data <- "olink_data.RData"
load(input_data)

select_batch <- "cbc03"

# specify data set
dat <- cbc03_olink_02 # specidy data to use (there are two Olink normalization methods),

binder <- dat[["binder"]] # get binder data
sinfo <- dat[["sinfo"]] # get sample info
npx <- dat[["npx"]] # get npx

txt <- paste0("Olink, ", unique(binder$Comment)) # add noramlization type

# only keep samples with 3 sample types (paired sample types)
select_ind <- names(which(table(sinfo$individual) == 3)) 

df <- npx[rownames(npx) %in% sinfo$id[sinfo$individual %in% select_ind], ]
df <- data.frame("id" = rownames(df), df)

# make one data frame with sample and npx info
df <- left_join(df, sinfo, by="id")

plot_data <- 
  df %>% 
  select(-id, -qc_warning, -sample_type) #%>% 

plot_data <- reshape(plot_data, idvar = "individual",
                     timevar = c("subtype"),
                     direction = "wide")

# Select the proteins (binders) to include in analysis:
selected_binders <- binder$id # include all binders

# extract spearman correlation values
correlation_df <- data.frame("binder_id" = selected_binders)

for(ab in selected_binders){
  
  # cDBS vs Plasma
  correlation_df$rho_cDBS_plasma[correlation_df$binder_id %in% ab] <-
    round(cor((plot_data[, paste0(ab, ".cDBS")]),
              (plot_data[, paste0(ab, ".Plasma")]),
              method="spearman", use="complete.obs"), 2)
}

##################### Extract p-values and fold change         ######################

# Create a data frame:
Pvalue_df <- data.frame("binder_id" = selected_binders)

for(ab in selected_binders){

  ################ Add the fold change values to the pvalue table:
  fold_change <- median(plot_data[, paste0(ab, ".cDBS")]) - median(plot_data[, paste0(ab, ".Plasma")])
  
  Pvalue_df$dist_cDBS_plasma[Pvalue_df$binder_id %in% ab] <- fold_change
  
  
  ############ Add paired t-test p.value value:
  # Compute paired t-test
  p.value_t <- t.test(plot_data[, paste0(ab, ".cDBS")],
  										plot_data[, paste0(ab, ".Plasma")], paired = TRUE)$p.value
  	
  Pvalue_df$p_ttest[Pvalue_df$binder_id %in% ab] <- p.value_t
  
}


############## FDR correciton:

Pvalue_df$FDR_ttest <- p.adjust(Pvalue_df$p_ttest, method="fdr", n=length(Pvalue_df[, 1]))

############ Color threshold:
for(ab in selected_binders){
	
	 FDR_value <- Pvalue_df$FDR_ttest[Pvalue_df$binder_id %in% ab]
	 log2_change <- Pvalue_df$dist_cDBS_plasma[Pvalue_df$binder_id %in% ab]

	if (FDR_value < 0.01 && log2_change >= 1) {
		Pvalue_df$threshold_FDR[Pvalue_df$binder_id %in% ab] <- "FDR and log2(FC)"
		Pvalue_df$dev_analyte_FDR[Pvalue_df$binder_id %in% ab] <- gene_names_table$Gene[grep(as.character(ab), gene_names_table$Assay)] #as.character(ab)
		
	} else if ( FDR_value < 0.01 && log2_change <= -1) {
		Pvalue_df$threshold_FDR[Pvalue_df$binder_id %in% ab] <- "FDR and log2(FC)"
		Pvalue_df$dev_analyte_FDR[Pvalue_df$binder_id %in% ab] <- gene_names_table$Gene[grep(as.character(ab), gene_names_table$Assay)] #as.character(ab)
		
	} else if ( FDR_value < 0.01 && log2_change > -1) {
		Pvalue_df$threshold_FDR[Pvalue_df$binder_id %in% ab] <- "FDR"
		
	} else if ( FDR_value < 0.01 && log2_change < 1) {
		Pvalue_df$threshold_FDR[Pvalue_df$binder_id %in% ab] <- "FDR"
		
	} else {
		Pvalue_df$threshold_FDR[Pvalue_df$binder_id %in% ab] <- "0"
		Pvalue_df$dev_analyte_FDR[Pvalue_df$binder_id %in% ab] <- ""
	}
	
}


#############								Volcanoplot
theme_set(theme_bw(base_size = 12))

volcano_plot_FDR <- ggplot(data = Pvalue_df, 
                        aes(x = dist_cDBS_plasma, 
                            y = -log10(FDR_ttest), 
                            colour=threshold_FDR)) +
  geom_point(size=3.5) +
  geom_text_repel(aes(label = dev_analyte_FDR), colour = "black", size = 8) +
  labs(x = "Δ NPX(DBS-plasma)", y = "-log10(FDR-value)") +
  geom_hline(yintercept = c(-log10(0.01)), linetype="dotted", color = c("black"), lwd=2) + 
  geom_vline(xintercept = c(-1, 1), linetype="dotted", lwd=2) + 
  xlim(-8,8) +
  theme_classic() +
  theme(axis.line = element_line(size = 1) , axis.text = element_text(size = 27, colour = "black"), axis.title = element_text(size = 30, face = "bold"), 
        plot.title = element_text(size = 18), legend.title = element_text(size = 14), legend.position = "none", axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm")) + 
  scale_colour_brewer(palette = "Set2")


##############      Histogram of the spearman rho correlation:

p_histo <- ggplot(correlation_df, aes(x=rho_cDBS_plasma)) + 
	geom_histogram(bins=17, color="black", fill="grey80") +
	labs(x = "Correlation [rho]", y = "Count [N]") +
	geom_vline(xintercept = 0, linetype="dotted", lwd = 2) +
  theme_classic() +
  xlim(-1,1) + 
  theme(axis.line = element_line(size = 1) , axis.text = element_text(size = 27, colour = "black"), axis.title = element_text(size = 30, face = "bold"), plot.title = element_text(size = 18), legend.title = element_text(size = 16), legend.position = "none", axis.ticks.length.x = unit(0.3, "cm"), axis.ticks.length.y = unit(0.3, "cm")) + scale_colour_brewer(palette = "Set2")

###                   Investigate IQR per protein and sample type

# Reshape the data frame to long format
plot_data_long <- melt(plot_data)

plot_data_long <- plot_data_long %>%
mutate(sample_type = case_when(
  grepl("\\.vDBS$", variable) ~ "vDBS",
  grepl("\\.cDBS$", variable) ~ "cDBS",
  grepl("\\.Plasma$", variable) ~ "Plasma",
  TRUE ~ NA_character_
)) %>%
  mutate(target = sub("\\.(vDBS|cDBS|Plasma)$", "", variable))


# Calculate IQR per protein per sample type
iqr_per_protein_sample <- plot_data_long %>%
  group_by(target, sample_type) %>%
  mutate(iqr = IQR(value)) %>%
  distinct(target, sample_type, iqr) %>%
  ungroup()

# Create an empty data frame to store the results
df_iqr <- data.frame(target = character(),
                     sample_type = character(),
                     iqr = numeric(),
                     stringsAsFactors = FALSE)

# Get unique combinations of target, sample_type, and individual
unique_combinations <- unique(plot_data_long[, c("target", "sample_type")])

# Loop through each unique combination and calculate IQR
for (i in seq_len(nrow(unique_combinations))) {
  target <- unique_combinations[i, "target"]
  sample_type <- unique_combinations[i, "sample_type"]
  
  # Subset the data for the current combination
  subset_data <- plot_data_long[plot_data_long$target == target &
                                  plot_data_long$sample_type == sample_type, ]
  
  # Calculate the IQR for the subset
  iqr_value <- IQR(subset_data$value, na.rm = TRUE)
  
  # Create a data frame for the current result
  result <- data.frame(target = target,
                       sample_type = sample_type,
                       iqr = iqr_value,
                       stringsAsFactors = FALSE)
  
  # Append the result to the df_iqr data frame
  df_iqr <- bind_rows(df_iqr, result)
}

# Filter rows with sample_type other than "vDBS"
df_filtered <- df_iqr[df_iqr$sample_type != "vDBS", ]

# Create a boxplot to compare IQR values between cDBS and Plasma sample types, and calculate p-value using t.test
ggplot(df_filtered, aes(x = sample_type, y = iqr, fill = sample_type)) +
  geom_boxplot() +
  labs(x = "Sample Type", y = "IQR Value", fill = "Sample Type") +
  ggtitle("Comparison of IQR Values between cDBS and Plasma") +
  theme_classic() +
  stat_compare_means(method = "t.test", label = "p.format", 
                     comparisons = list(c("cDBS", "Plasma"))) +
  theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(size = 18),legend.position = "none")



## Calculate delta IQR:
df_delta_iqr <- df_filtered %>%
  pivot_wider(names_from = sample_type, values_from = iqr) %>%
  mutate(delta_iqr = Plasma - cDBS)

max_index <- which.max(df_delta_iqr$delta_iqr)
min_index <- which.min(df_delta_iqr$delta_iqr)

delta_iqr_plot <- ggplot(df_delta_iqr, aes(x= reorder(target, delta_iqr), y = delta_iqr) ) +
  geom_point() +
  labs(x = "Target", y = "Delta IQR (Plasma - DBS)") +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey", size=1) +
  geom_text(
    data = df_delta_iqr[c(max_index, min_index), ],
    aes(label = target),
    color = "black",
  nudge_x = c(-3, 3) ,
  nudge_y = c(-0.05, 0.05) ,
  size = 4
  ) +
  theme(axis.line = element_line(size = 1, colour = "black"), axis.text = element_text(size = 0, angle = 90), axis.title = element_text(size = 16), plot.title = element_text(size = 18))


# master_cov_olink.R
#
# description:  Master script for pipelining analysis of covid19 Olink DBS data
#               Sources/Renders other scripts needed for analysis
#         
#               Panels
#               - Olink Target 96 Cardiometabolic(v.3603)
#               - Olink Target 96 Cardiovascular III(v.6113)
#               - Olink Target 96 Metabolism(v.3403)
#
#
# required_packages: rmarkdown
#
# creator: Tea Dodig-Crnkovic
# contributors: Tea Dodig-Crnkovic, Leo Dahl
######################################################################

# clean up workspace
rm(list=ls())

# packages
library(rmarkdown)

# set working directory
setwd("...")

# load_cov_npx.R contain functions for loading and transforming data

####### Read in and summarize data #######

# summary of covid19 DBS samples
render(input = "summary_covid_dbs_samples.Rmd",
       output_file = paste0("../results/sample_summary/",
                            format(Sys.time(),"%Y-%m-%d_%H%M%S"),
                            "_sample_summary.html"))

####### Quality control #######

# sample CV
render(input = "cv_samples.Rmd",
       output_file = paste0("../results/cv/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_sample_cv.html"))

# sample correlation
render(input = "sample_correlation_dbs.Rmd", 
       output_file = paste0("../results/correlation_plots/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_correlation_plots.html"))

# Binder IQR
render(input = "binder_iqr.Rmd",
       output_file = paste0("../results/binder_iqr/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_binder_iqr.html"))

# Outlier check, based on Olink QC plot of median vs IQR
render(input = "outlier_median_iqr.Rmd",
       output_file = paste0("../results/outliers/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_outliers_median_iqr.html"))

# Global overview + outlier check
render(input = "global_overview.Rmd",
       output_file = paste0("../results/global_overview/", 
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_global_overview.html"))

####### Association tests #######

# Lasso regression
render(input = "lasso_regression.Rmd", 
       output_file = paste0("../results/lasso_regression/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_lasso_regression.html"))

# Logistic regression against metadata and serostatus, with dataset 3 included
render(input = "association_tests_all_datasets.Rmd",
       output_file = paste0("../results/association_tests_lm_kw/",
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_association_logreg_all.html"))

# HeVa with all datasets
render(input = "heva_analysis.Rmd", 
       output_file = paste0("../results/heva/", 
                            format(Sys.time(), "%Y-%m-%d_%H%M%S"),
                            "_heva_analysis.html"))

# Correlation between antigens and olink binders
render(input = "population_antigen_relation.Rmd", 
			 output_file = paste0("../results/population_antigen_relation/",
			 										 format(Sys.time(), "%Y-%m-%d_%H%M%S"),
			 										 "_population_antigen_relation.html"))

####### Longitudinal analysis #######

# Longitudinal analysis of one individual
render(input = "longitudinal.Rmd",
			 output_file = paste0("../results/longitudinal/",
			 										 format(Sys.time(), "%Y-%m-%d_%H%M%S"),
			 										 "_longitudinal.html"))

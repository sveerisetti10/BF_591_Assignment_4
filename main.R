#Load Packages
library(DESeq2)
library(tidyverse)

#Sriramteja Veerisetti 
#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x m) tibble with a 'gene' column followed by
#' sample names as column names. 
#' 
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#' 
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(verse_counts){
  read_tsv(verse_counts)
  return()
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x m) tibble of raw read counts
#'
#' @return tibble: a (n x m) tibble of raw reads with genes that have 
#' zero variance across samples removed
#' 
#' @note (g >= n)
#' 
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_count) {
  #we want to omit the first column, which is why we do -c(1)
  #for the margin, we only want the rows to be manipulated, so we say 1
  #for the FUN we want to calculate the variance of the dataframe so var = variance
  gene_vars <- apply(verse_count[,-c(1)], c(1), var)
  #we do not want our verse_count containing any rows that have only 0 values. != means "does not equal" 
  return(verse_counts[gene_vars!=0,])
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point 
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(timepoint_replicate) {
  timepoint <- "[A-Z][a-z,0-9]"
  analysis <- stringr::str_extract(timepoint_replicate, timepoint)
  
  return(analysis)
}


#' Grab sample replicate number from sample name
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#' 
#' @note you may choose to return numeric values instead of strings here
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(timepoint_replicate) {
  str_split(timepoint_replicate, "_")[[1]][2] %>%
  return()

}


#' Generate sample-level metadata from sample names. 
#' 
#' Will include columns named "sample", "timepoint", and "replicate" that store 
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample", 
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint" 
#' stores sample time points; and "replicate" stores sample replicate
#' 
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  replication <- sapply(sample_names, sample_replicate)
  timepointt <- sapply(sample_names, timepoint_from_sample)
  
  return(tibble(Sample = sample_names, 
                Timepoint = timepointt,
                Replicate = replication))
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x m) tibble of raw read counts.
#'
#' @return named vector: numeric vector of read totals from each sample
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  count_data[colnames(count_data)!='gene'] %>%
    summarize(dplyr::across(.cols = everything(), ~sum(., is.na(.), 0))) %>%
    return()
}

#' Normalize raw count data to counts per million WITH pseudocounts using the 
#' following formula:
#'     (count + 1) / ((sample_library_size/10^6) + 1)
#'
#' @param count_data tibble: a (n x m) matrix of raw read counts.
#'
#' @return tibble: a (n x m) matrix with read count normalized to counts
#' per million
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) { 
  count_data %>%
    #we have to unlist the function to get it from a list to a vector 
    #otherwise this does not work 
    mutate((.[c(-1)] + 1) / (unlist(get_library_size (.[c(-1)]) /1000000)) +1) %>%
    return ()
  
}


#' Normalize raw count data using DESeq2
#'
#'
#' @param count_data tibble: a (n x m) matrix of raw reads
#' @param meta_data tibble: sample-level information tibble containing time point,
#' sample, and replicate information.
#' @param design_formula formula: formula of comparision of interest
#'
#' @return tibble: a (n x m) tibble of DESeq2 normalized count data.
#'
#' @example ' `deseq_normalize(count_data, meta_data, ~ timepoint)`

deseq_normalize <- function(count_data, meta_data, design_formula) {
  count_select <- select(count_data, -c(gene))
  DEsequece <- DESeqDataSetFromMatrix(
    countData = count_select,
    colData = meta_data,
    design = ~1) 
  
  DEsequece <- estimateSizeFactors(DEsequece)
  Dnormalizedcounts <- as_tibble(counts(DEsequece, normalized = TRUE)) %>%
    mutate(gene = count_data$gene) %>%
    relocate(gene) %>%
  
  return()
}

#make sure to load in DESeq2
#count_select <- select(count_data)
#DEsequece <- DESeqDataSetFromMatrix(
#countData = count_select,
#colData = meta_data,
#design = design_formula)

#DEsequece <- estimateSizeFactors(DEsequece)
#Dnormalizedcounts <- as_tibble(counts(DEsequece, normalized = TRUE)) %>%
#mutate(gene = count_data$gene) %>%
#relocate(gene) %>%

#return()


#' Perform and plot PCA using processed data.
#' 
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs. 
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(no_header, sample_meta, title="Raw Counts PCA") {
  comp <- prcomp(t(no_header),
                 center = FALSE, 
                 scale = TRUE)
  plot <- sample_meta 
  plot$PC1 <- comp$x[ ,1]
  plot$PC2 <- comp$x[ ,2]
  pca_variance <- comp$sdev^2
  pca_variance_percentage <- pca_variance/sum(pca_variance)
  pca_plot <- ggplot(plot, aes(x = PC1, y = PC2, col = Timepoint)) +
    geom_point() +
    labs(x = round(pca_variance_percentage[1]) * 100,
         y = round(pca_variance_percentage[2]) * 100,
         title = "Raw Counts PCA") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  
  return(pca_plot)
}


#' Plot gene count distributions for each sample using boxplots.
#' 
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(no_header, scale_y_axis=FALSE, title="Sample Distributions (Log10(Raw Counts))") {
  
  pivot <- pivot_longer(no_header, names_to = 'sample', 
                        cols = colnames(no_header), values_to = 'counts') %>%
    mutate(sample = factor(sample, levels = colnames(no_header)))
  
  distribution_plot <- ggplot(pivot, aes(x = sample, y = counts, col = sample)) +
    geom_boxplot() + 
    labs(title = "Sample Distributions (Log10(Raw Counts))") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  
  if (scale_y_axis) {
    pivot <- filter(pivot, counts != 0) 
    distribution_plot <- distribution_plot + ggplot2::scale_y_log10()
     }
  
  return(distribution_plot)
    
  }  
    
# return(verse_counts[gene_vars!=0,])

#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(no_header, scale_y_axis=FALSE, title="Variance vs. Mean (Log10(Raw Counts))") {
  
  xaxis_mean <- apply(no_header, 1, mean)
  yaxis_variance <- apply(no_header, 1, var)
  plot_tibble <- tibble::tibble(mean = xaxis_mean, variance = yaxis_variance)
  variance_mean_plot <- ggplot(plot_tibble, aes(x = mean, y = variance)) +
    geom_point() +
    geom_line() +
    labs(x = "mean", y = "Variance", title = "Variance vs. Mean (Log10(Raw Counts))") +  
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  
  if (scale_y_axis) {
    variance_mean_plot <- variance_mean_plot + ggplot2::scale_y_log10()
  }
  
  return(variance_mean_plot)
}

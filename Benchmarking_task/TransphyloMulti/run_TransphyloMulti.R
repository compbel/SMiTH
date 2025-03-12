# Author: Maryam KafiKang
# Created: 2024
# Description: This R script runs the TransPhyloMulti analysis to infer transmission networks from
#             phylogenetic trees. It processes tree data to compute transmission probabilities and
#             intermediate matrices for detailed analysis. The workflow includes:
#             1. Tree Processing:
#                - Reads and modifies Newick format phylogenetic trees
#                - Maps host labels for analysis
#             2. Transmission Inference:
#                - Infers transmission trees using TransPhyloMulti
#                - Computes transmission probability matrices
#                - Calculates intermediate matrices for additional insights
#             3. Output Generation:
#                - Saves transmission matrices and intermediate results
#                - Records execution time for performance evaluation
# 
# Input:
#     - Newick format phylogenetic trees
#     - Command line arguments for shape parameter and sample index
# 
# Output:
#     - Transmission probability matrices
#     - Intermediate matrices for detailed analysis
#     - Execution time statistics

#devtools::install_github("DrJCarson/TransPhyloMulti")
.libPaths("~/rlibs")
library('TransPhyloMulti')
library(ape)
library(stringr)
source("MyUtils.R")
source("calculate_w.R")


computeMatTDistM = function(record,burnin=0.5)
{
  #Remove burnin
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  t1=extractTTreeM(record[[1]]$ctree)
  n=length(unique(t1$obs[, 2])) #Number of sampled individuals
  mat=matrix(0,n,n)
  colnames(mat)<-1:n
  rownames(mat)<-1:n
  
  for (i in 1:length(record))
  {
    ttree=extractTTreeM(record[[i]]$ctree)$ttree
    for (a in 2:n) for (b in 1:a) {
      aa=a
      bb=b
      count=0
      while (aa!=bb) {
        if (ttree[aa,1]>ttree[bb,1]) aa=ttree[aa,3] else bb=ttree[bb,3]
        count=count+1
      }
      mat[a,b]=mat[a,b]+count/length(record)
      mat[b,a]=mat[a,b]
    }
  }
  return(mat)
}

main <- function(tree, host, shape, scale){
  newick_str_modified <- replace_underscores_in_file(tree)
  
  dated_tree <- read.tree(text = newick_str_modified)
  
  start_time <- Sys.time()
  pt=ptreeFromPhyloM(dated_tree,lubridate::decimal_date(as.Date('2022/7/1')))
  #plot(pt)
  
  #res <- inferTTreeM(pt, w.shape = shape, w.scale = scale, obs.end = lubridate::decimal_date(as.Date('2024/7/1')))
  
  res <- inferTTreeM(pt, w.shape = shape, w.scale = scale ,obs.end = lubridate::decimal_date(as.Date('2025/11/1')))
  end_time <- Sys.time()
  duration <- as.numeric(end_time - start_time, units = "secs")
  #plot(res)
  
  # Print the results
  #print(res)
  
  #plot(res[[length(res)]]$ctree)
  # Matrix of transmission probabilities from each host (row) to any other (column).
  mat <- computeMatWIWM(res)
  
  mat_intermediate <- computeMatTDistM(res)
  print(mat_intermediate)
  
  return(list(
    transmission = mat,
    mat_intermediate = mat_intermediate,
    running_time = duration
  ))
}


input_dir <- "tree"
output_dir <- "TransphyloMultiResultMat/"

args <- commandArgs(trailingOnly = TRUE)
shape <- as.numeric(args[1])
i <- as.numeric(args[2])


tryCatch({
  data <- paste0(input_dir,"/PrunedTree_E1_",i,".nwk")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (file.exists(input_dir)) {
    
    new_tree <- map_host_labels_tree(data)
    tree <- new_tree$new_tree_string
    host <- new_tree$host_mapping
    host_number <- length(host)
    
    #tree_new <- readLines(data, warn = FALSE)
    
    
    # Extract all unique host labels
    #host_labels <- str_extract_all(tree_string, "\\d+_N\\d+")[[1]]
    
    if (host_number>= 5 & host_number <= 30 ) {
      print(paste("#host:", host_number, "#test:", i))
      
      ###########################################################################################################
      #####this section assign w$shape and w$scale based on the this formula : w$shaep *w$scale = 1/infectionRate
      ###########################################################################################################
      target_mean <- 20 #(1/infection rate) = (1/ 0.05)
      scale <- target_mean / shape
      print(paste("shape:", shape, "scale:", scale))
      output_base <- paste0(output_dir, "Shape_", shape, "_Scale_", scale,"/wiwMatrix")
      output_base_inter <- paste0(output_dir, "Shape_", shape, "_Scale_", scale, "/intermediateMatrix")
      output_base_dir <- paste0(output_base, "/FAVITES_output_expSI_contemp_T200_N100_E1_", i, "_transphylo.txt")
      output_base_intermediate_mat <- paste0(output_base_inter, "/FAVITES_output_expSI_contemp_T200_N100_E1_", i, "_interMatrix.txt")

      if (!file.exists(output_base_dir)) {
      tryCatch(
        {
          if (!dir.exists(output_base)) {
            dir.create(output_base, recursive = TRUE)
          }
          if (!dir.exists(output_base_inter)) {
            dir.create(output_base_inter, recursive = TRUE)
          }

          result <- main(tree, host, shape, scale)
          reverse_host <- setNames(names(host), unlist(host))
          mat <- result$transmission
          # Rename rows and columns based on the host mapping
          rownames(mat) <- reverse_host[as.character(1:nrow(mat))]
          colnames(mat) <- reverse_host[as.character(1:ncol(mat))]
          
          mat_inter <- result$mat_intermediate
          # Rename rows and columns based on the host mapping
          rownames(mat_inter) <- reverse_host[as.character(1:nrow(mat_inter))]
          colnames(mat_inter) <- reverse_host[as.character(1:ncol(mat_inter))]
         
          write.table(mat, file = output_base_dir, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
          write.table(mat_inter, file = output_base_intermediate_mat, row.names = TRUE, col.names = NA, quote = FALSE)
          write(sprintf("Running time: %.4f seconds", result$running_time), file = output_base_dir, append = TRUE)
          }, error = function(e) {
          print(paste("Error in assigning w using loop in the test", i, ":", conditionMessage(e)))
        })
      }
    }
  }
}, error = function(e) {
  print(paste("Error in test", i, ":", conditionMessage(e)))
})





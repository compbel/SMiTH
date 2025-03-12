map_host_labels_tree <- function(input_file_path) {
  # Read the tree string from the input file
  tree_string <- readLines(input_file_path, warn = FALSE)
  tree_string <- paste(tree_string, collapse = " ")
  
  # Extract all unique host labels
  host_labels <- str_extract_all(tree_string, "\\d+_N\\d+")[[1]]
  unique_hosts <- unique(sub("_N\\d+", "", host_labels))
  
  # Create a mapping from original host labels to sequential numbers starting from 1
  host_mapping <- setNames(seq_along(unique_hosts), unique_hosts)
  
  # Create a function to replace the host labels in the tree string
  replace_labels <- function(label) {
    host <- sub("_N\\d+", "", label)
    new_host <- host_mapping[host]
    return(sub(host, new_host, label))
  }
  
  # Apply the replacement function to all labels in the tree string
  new_tree_string <- tree_string
  for (label in host_labels) {
    new_tree_string <- str_replace_all(new_tree_string, label, replace_labels(label))
  }
  
  # Create a dictionary for the mapping
  mapping_dict <- as.list(host_mapping)
  
  
  return(list(new_tree_string = new_tree_string, host_mapping = mapping_dict))
}


extract_max_node_value_from_file <- function(file_path) {
  # Read the phylogeny tree from the file
  tree_string <- readLines(file_path)
  
  # Collapse the lines into a single string (in case the file has multiple lines)
  tree_string <- paste(tree_string, collapse = "")
  
  # Use a regular expression to match the numerical values after the colon
  matches <- gregexpr(":[0-9]+\\.[0-9]+", tree_string)
  values <- regmatches(tree_string, matches)
  
  # Extract the numerical values and convert them to numeric type
  numeric_values <- as.numeric(sub(":", "", unlist(values)))
  
  # Return the maximum value
  max_value <- max(numeric_values)
  
  return(numeric_values)
}


replace_underscores_in_file <- function(newick_str) {
  # Read the Newick file
  
  newick_str <- paste(newick_str, collapse = "")
  #gsub("(_\\d+):", ":", "10_N1_100:")
  #gsub("([a-zA-Z0-9]+)_([a-zA-Z0-9]+)", "\\1.\\2", "10_N1:")
  # Replace underscores with periods
  newick_str_modified <- gsub("_", ".", newick_str)
  
  return(newick_str_modified)
}


find_transmission <- function(mat){
  result <- list()
  for (i in 1:nrow(mat)){
    val <- max(mat[i, ])
    col <- which(mat[i, ] == val)
    result[[i]] <- c(row = i, column = col, max_value = val)
  }
  return(result)
}

calc_prim_obs <- function(ptree) {
  
  host <- ptree$host
  ptree <- ptree$ptree
  
  prim_times <- numeric(max(host))
  
  for (i in 1:max(host)) {
    
    prim_times[i] <- min(ptree[which(host == i), 1])
    
  }
  
  return(prim_times)
  
}


find_transmission <- function(mat, host) {
  # Create a reverse mapping from numbers to host labels
  reverse_host <- setNames(names(host), unlist(host))
  
  # Initialize vectors to store the results
  rows <- c()
  columns <- c()
  max_values <- c()
  
  for (i in 1:nrow(mat)) {
    val <- max(mat[i, ])
    col <- which(mat[i, ] == val)
    
    row_host <- reverse_host[as.character(i)]
    col_host <- reverse_host[as.character(col)]
    
    rows <- c(rows, row_host)
    columns <- c(columns, col_host)
    max_values <- c(max_values, val)
  }
  
  # Create a dataframe from the vectors
  result_df <- data.frame(
    row = rows,
    column = columns,
    max_value = max_values,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

whole_matrix <- function(mat, host) {
  # Create a reverse mapping from numbers to host labels
  reverse_host <- setNames(names(host), unlist(host))
  
  # Initialize vectors to store the results
  rows <- c()
  columns <- c()
  values <- c()
  
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      row_host <- reverse_host[as.character(i)]
      col_host <- reverse_host[as.character(j)]
      
      rows <- c(rows, row_host)
      columns <- c(columns, col_host)
      values <- c(values, mat[i, j])
    }
  }
  
  # Create a dataframe from the vectors
  result_df <- data.frame(
    row = rows,
    column = columns,
    value = values,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}


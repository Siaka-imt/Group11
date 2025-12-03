# data_utils.R - Real Datasets Version
# Handles downloading, parsing, and matrix conversion for real-world data

load_data <- function(dataset_name) {
  
  # ==========================================
  # 1. Real Dataset: MovieLens 100k (Sparse)
  # ==========================================
  if (dataset_name == "MovieLens 100k (Real)") {
    
    # Define local filename
    ml_file <- "ml-100k-u.data"
    
    # Check if file exists, if not, download from official source
    if (!file.exists(ml_file)) {
      message("Downloading MovieLens 100k dataset...")
      tryCatch({
        download.file("https://files.grouplens.org/datasets/movielens/ml-100k/u.data", 
                      destfile = ml_file, method = "auto", quiet = TRUE)
      }, error = function(e) {
        stop("Failed to download MovieLens data. Please check your internet connection.")
      })
    }
    
    # Read data (Format: User ID | Item ID | Rating | Timestamp)
    # This imports as a data frame (Long format)
    raw_data <- read.table(ml_file, sep = "\t", header = FALSE, 
                           col.names = c("user_id", "item_id", "rating", "timestamp"))
    
    # ==== Data Transformation: Long Format to Sparse Matrix (User x Item) ====
    # Logic: Use xtabs to pivot data.
    # Optimization: We take the top 200 users and 300 movies to ensure 
    # the pure R implementation runs smoothly during the demo.
    
    top_users <- unique(raw_data$user_id)[1:200]
    top_items <- unique(raw_data$item_id)[1:300]
    
    subset_data <- raw_data[raw_data$user_id %in% top_users & 
                              raw_data$item_id %in% top_items, ]
    
    # Convert to matrix (Missing ratings become 0)
    X_sparse <- xtabs(rating ~ user_id + item_id, data = subset_data)
    X <- as.matrix(X_sparse)
    
    # Enforce double precision for numerical stability
    storage.mode(X) <- "double"
    
    return(list(
      X = X,
      title = "MovieLens 100k (Subset)",
      description = paste("Real user ratings from MovieLens.",
                          "\nDimensions:", nrow(X), "Users x", ncol(X), "Movies.",
                          "\nSparsity:", round(sum(X==0)/length(X)*100, 1), "% zeros.",
                          "\n(Subsampled for performance efficiency)")
    ))
    
    # ==========================================
    # 2. Real Dataset: Volcano Topography (Dense)
    # ==========================================
  } else if (dataset_name == "Volcano Topography (Real)") {
    
    # Load R's built-in Volcano dataset
    data(volcano)
    
    # Volcano is an 87x61 matrix
    X <- as.matrix(volcano)
    
    # Normalize to 0-100 range for better NMF compatibility (NMF requires non-negative)
    X <- X - min(X) 
    
    # Enforce double precision
    storage.mode(X) <- "double"
    
    return(list(
      X = X,
      title = "Maunga Whau Volcano",
      description = paste("Topographic information on Maunga Whau volcano.",
                          "\nDimensions:", nrow(X), "x", ncol(X), "grid points.",
                          "\nType: Dense, Non-negative matrix.",
                          "\n(Ideal for demonstrating image/surface compression)")
    ))
  }
}


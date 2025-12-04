# data_utils.R - Real & Large Scale Version
# Includes MovieLens, Volcano, and a Large Synthetic Benchmark

load_data <- function(dataset_name) {
  
  # ==========================================
  # 1. Real Dataset: MovieLens 100k (Sparse)
  # ==========================================
  if (dataset_name == "MovieLens 100k (Real)") {
    ml_file <- "ml-100k-u.data"
    if (!file.exists(ml_file)) {
      message("Downloading MovieLens 100k dataset...")
      tryCatch({
        download.file("https://files.grouplens.org/datasets/movielens/ml-100k/u.data", 
                      destfile = ml_file, method = "auto", quiet = TRUE)
      }, error = function(e) {
        stop("Failed to download MovieLens data. Please check your internet connection.")
      })
    }
    raw_data <- read.table(ml_file, sep = "\t", header = FALSE, 
                           col.names = c("user_id", "item_id", "rating", "timestamp"))
    top_users <- unique(raw_data$user_id)[1:200]
    top_items <- unique(raw_data$item_id)[1:300]
    subset_data <- raw_data[raw_data$user_id %in% top_users & 
                              raw_data$item_id %in% top_items, ]
    X_sparse <- xtabs(rating ~ user_id + item_id, data = subset_data)
    X <- as.matrix(X_sparse)
    storage.mode(X) <- "double"
    
    return(list(
      X = X,
      title = "MovieLens 100k (Subset)",
      description = paste("Real user ratings from MovieLens.",
                          "\nDimensions:", nrow(X), "Users x", ncol(X), "Movies.",
                          "\nSparsity:", round(sum(X==0)/length(X)*100, 1), "% zeros.")
    ))
    
    # ==========================================
    # 2. Real Dataset: Volcano Topography (Dense)
    # ==========================================
  } else if (dataset_name == "Volcano Topography (Real)") {
    data(volcano)
    X <- as.matrix(volcano)
    X <- X - min(X) 
    storage.mode(X) <- "double"
    
    return(list(
      X = X,
      title = "Maunga Whau Volcano",
      description = paste("Topographic information on Maunga Whau volcano.",
                          "\nDimensions:", nrow(X), "x", ncol(X), "grid points.",
                          "\nType: Dense, Non-negative matrix.")
    ))
    
    # ==========================================
    # 3. Large Scale Benchmark (Synthetic)
    # ==========================================
  } else if (dataset_name == "Large Synthetic (Benchmark)") {
    
    # set large scale：2000 x 2000 to show significant speed difference
    n <- 2000
    m <- 2000
    true_rank <- 50
    
    message("Generating large synthetic matrix (2000x2000)...")
    
    # generate low rank structure A = U * V^T
    U <- matrix(rnorm(n * true_rank), nrow = n)
    V <- matrix(rnorm(m * true_rank), nrow = m)
    
    # add signal decay to make it more real
    S <- diag(seq(100, 1, length.out = true_rank))
    
    # core signal
    Signal <- U %*% S %*% t(V)
    
    # add noise
    Noise <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n)
    
    X <- Signal + Noise
    
    # ensure non-negative
    X <- X - min(X) + 0.1 
    
    storage.mode(X) <- "double"
    
    return(list(
      X = X,
      title = "Large Scale Synthetic Benchmark",
      description = paste("A large dense matrix constructed for performance testing.",
                          "\nDimensions: 2000 x 2000.",
                          "\nTrue Rank: 50 + Gaussian Noise.",
                          "\nDesigned to highlight the speed of Randomized SVD.")
    ))
  }
}


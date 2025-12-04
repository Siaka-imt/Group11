# algorithms.R - Fixed Return Values for Compatibility

# Helper function: Calculate RMSE
calc_rmse <- function(X, X_recon) {
  sqrt(mean((X - X_recon)^2))
}

# Helper function: Calculate relative Frobenius error
calc_rel_error <- function(X, X_recon) {
  norm(X - X_recon, type = "F") / norm(X, type = "F")
}

# ==========================================
# 1. SVD Algorithms
# ==========================================

algo_svd_normal <- function(X, k) {
  if(!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  s <- svd(X, nu = k, nv = k)
  
  # Compute reconstruction
  D <- diag(s$d[1:k], nrow = k, ncol = k)
  recon <- s$u %*% D %*% t(s$v)
  
  list(recon = recon, components = s$u, name = "Deterministic SVD")
}

algo_svd_randomized <- function(X, k, n_oversamples = 10, n_iter = 2) {
  if(!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); m <- ncol(X)
  p <- min(k + n_oversamples, m)
  
  Omega <- matrix(rnorm(m * p), nrow = m)
  Y <- X %*% Omega
  
  for (j in 1:n_iter) {
    Y <- X %*% crossprod(X, Y)
  }
  
  qr_res <- qr(Y)
  Q <- qr.Q(qr_res)
  B <- crossprod(Q, X)
  s <- svd(B, nu = k, nv = k)
  U <- Q %*% s$u
  
  recon <- U %*% diag(s$d[1:k], nrow = k) %*% t(s$v)
  
  list(recon = recon, components = U, name = "Randomized SVD")
}

# ==========================================
# 2. NMF Algorithms
# ==========================================

# Core update logic (Multiplicative Updates)
run_nmf_core <- function(X, W, H, max_iter = 200, tol = 1e-4) {
  errors <- c()
  for (i in 1:max_iter) {
    WH <- W %*% H
    H <- H * (crossprod(W, X)) / (crossprod(W, WH) + 1e-9)
    WH <- W %*% H
    W <- W * (X %*% t(H)) / (WH %*% t(H) + 1e-9)
    
    if (i %% 10 == 0) {
      e <- norm(X - W %*% H, type = "F")
      errors <- c(errors, e)
      if (length(errors) > 1 && abs(errors[length(errors)] - errors[length(errors)-1]) < tol) break
    }
  }
  list(W = W, H = H, errors = errors)
}

algo_nmf_normal <- function(X, k) {
  X <- abs(X)
  n <- nrow(X); m <- ncol(X)
  W_init <- matrix(runif(n*k), nrow=n)
  H_init <- matrix(runif(k*m), nrow=k)
  
  res <- run_nmf_core(X, W_init, H_init)
  
  list(recon = res$W %*% res$H, components = res$W, name = "Standard NMF (Random Init)")
}

# ==========================================
# Randomized HALS NMF (Paper Implementation) [FIXED]
# ==========================================
algo_nmf_randomized <- function(X, k, max_iter = 100, tol = 1e-4) {
  
  # 1. Parameter Setup
  p <- 10             
  l <- k + p          
  n <- ncol(X)
  m <- nrow(X)
  
  # 2. Randomized Compression
  Omega <- matrix(rnorm(n * l), nrow = n, ncol = l)
  Y <- X %*% Omega
  
  q <- 2
  for(i in 1:q) {
    Y <- X %*% crossprod(X, Y)
  }
  
  qr_res <- qr(Y)
  Q <- qr.Q(qr_res) 
  B <- crossprod(Q, X)
  
  # 3. Initialization
  W_tilde <- matrix(runif(l * k), nrow = l, ncol = k)
  H <- matrix(runif(k * n), nrow = k, ncol = n)
  
  # 4. HALS Iterations
  errors <- c()
  eps <- 1e-16
  
  for(iter in 1:max_iter) {
    
    # Update H
    V <- crossprod(W_tilde, W_tilde)
    Cross <- crossprod(W_tilde, B)
    for (j in 1:k) {
      H[j, ] <- pmax(0, H[j, ] + (Cross[j, ] - drop(V[j, , drop=FALSE] %*% H)) / (V[j, j] + eps))
    }
    
    # Update W_tilde
    T_mat <- B %*% t(H)
    S <- H %*% t(H)
    
    for (j in 1:k) {
      W_tilde[, j] <- W_tilde[, j] + (T_mat[, j] - drop(W_tilde %*% S[, j])) / (S[j, j] + eps)
      
      # Projection to High Dimension to enforce non-negativity
      w_high <- Q %*% W_tilde[, j, drop=FALSE]
      w_high[w_high < 0] <- 0
      W_tilde[, j] <- crossprod(Q, w_high)
    }
    
    # Simple convergence check
    if (iter %% 5 == 0) {
      curr_error <- norm(B - W_tilde %*% H, "F")
      errors <- c(errors, curr_error)
      if(length(errors) > 1 && abs(diff(tail(errors, 2)))/tail(errors, 1) < tol) break
    }
  }
  
  # 5. Final Result Construction
  W_final <- Q %*% W_tilde
  W_final[W_final < 0] <- 0 
  
  # ==== FIXED HERE: Calculate reconstruction and assign components ====
  recon <- W_final %*% H
  
  list(
    recon = recon,           # Required for plotting and error calc
    components = W_final,    # Required for 'Deep Dive' plots
    W = W_final, 
    H = H, 
    iter = iter, 
    name = "Randomized HALS (Paper Impl.)"
  )
}

# ==========================================
# 3. CUR Matrix Decomposition
# ==========================================

compute_leverage_scores <- function(X, k) {
  s <- svd(X, nu = k, nv = k)
  lev_col <- rowSums(s$v[, 1:k]^2) / k
  lev_row <- rowSums(s$u[, 1:k]^2) / k
  list(col = lev_col, row = lev_row)
}

algo_cur_normal <- function(X, k) {
  c_num <- k + 5
  r_num <- k + 5
  c_num <- min(c_num, ncol(X)); r_num <- min(r_num, nrow(X))
  
  probs <- compute_leverage_scores(X, k)
  
  col_idx <- order(probs$col, decreasing = TRUE)[1:c_num]
  row_idx <- order(probs$row, decreasing = TRUE)[1:r_num]
  
  C <- X[, col_idx, drop = FALSE]
  R <- X[row_idx, , drop = FALSE]
  
  pinv <- function(M) {
    s <- svd(M); tol <- 1e-5
    d_inv <- 1 / s$d; d_inv[s$d < tol] <- 0
    s$v %*% diag(d_inv, nrow=length(d_inv)) %*% t(s$u)
  }
  
  U_core <- pinv(C) %*% X %*% pinv(R)
  recon <- C %*% U_core %*% R
  
  list(recon = recon, components = C, name = "Deterministic CUR (Top Leverage)")
}

# Randomized: Randomized CUR (Probability sampling based on leverage scores)
algo_cur_randomized <- function(X, k) {
  c_num <- k + 5
  r_num <- k + 5
  c_num <- min(c_num, ncol(X)); r_num <- min(r_num, nrow(X))
  
  probs <- compute_leverage_scores(X, k)
  
  # Probabilistic sampling
  col_idx <- sample(1:ncol(X), c_num, prob = probs$col, replace = TRUE)
  row_idx <- sample(1:nrow(X), r_num, prob = probs$row, replace = TRUE)
  
  C <- X[, col_idx, drop = FALSE]
  R <- X[row_idx, , drop = FALSE]
  
  pinv <- function(M) {
    s <- svd(M); tol <- 1e-5
    d_inv <- 1 / s$d; d_inv[s$d < tol] <- 0
    s$v %*% diag(d_inv, nrow=length(d_inv)) %*% t(s$u)
  }
  
  U_core <- pinv(C) %*% X %*% pinv(R)
  recon <- C %*% U_core %*% R
  
  list(recon = recon, components = C, name = "Randomized CUR (Weighted Sampling)")
}


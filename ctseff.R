ctseff <- function(a0, learner, y, a, x, gamma = 1, nsplits = 5,
                   double_split = TRUE, ...) {
  
  params <- list(...)
  
  n <- length(y)
  ng <- length(gamma)
  
  tau.u <-  gamma / (1 + gamma)
  tau.l <- 1 / (1 +  gamma)
  
  c1.l <-  gamma - 1 / gamma
  c2.l <- 1 / gamma
  c1.u <- - c1.l
  c2.u <- 1 / c2.l
  
  s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  
  est.l <- est.u <- replicate(length(learner), 
                              array(NA, dim = c(length(a0), 3, nsplits, ng),
                                    dimnames = list(NULL, 
                                                    c("est", "ci.lo", "ci.hi"),
                                                    paste0("split", 1:nsplits),
                                                    gamma)),
                              simplify = FALSE)
  names(est.l) <- names(est.u) <- learner
  pseudo.u.folds <- pseudo.l.folds <- a.te.folds <- 
    vector(mode = "list", length = nsplits)
  pihat_folds <- vector(mode = "list", length = nsplits)
  names(est.u) <- names(est.l) <- learner
  ps <- params[["ps"]]
  ps.method <- params[["ps.method"]]
  
  sl.lib <- params[["sl.lib"]]
  
  if(any(learner  == "dr")) {
    
    mu.method <- params[["mu.method"]]
    v.method <- params[["v.method"]]
    q.method <- params[["q.method"]]
    drl.method <- params[["drl.method"]]
    
    mu <- params[["mu"]]
    q <- params[["q"]]
    v <- params[["v"]]
    drl <- params[["drl"]]
    
  }
  
  if(length(unique(a)) < 100) a.vals <- unique(sort(a))
  else a.vals <- seq(min(a), max(a), length.out = 100)
  
  
  for(k in 1:nsplits) {
    
    if(nsplits == 1) {
      tr.idx <- te.idx <- rep(TRUE, n)
    } else {
      tr.idx <- s != k
      te.idx <- s == k
    }

    n.tr <- sum(tr.idx)
    n.te <- sum(te.idx)
    pihat <- pahat <- v.mhat <- rep(NA, n.te)
    
    x.tr <- x[tr.idx, , drop = FALSE]
    x.te <- x[te.idx, , drop = FALSE]
    
    y.tr <- y[tr.idx]
    y.te <- y[te.idx]
    
    a.tr <- a[tr.idx]
    a.te <- a[te.idx]
    
    a.fctr <- factor(a.tr, ordered = TRUE)
    a0.fctr <- factor(a.vals, ordered = TRUE)
    pihat_fit <- nnet::multinom(a ~ ., data = cbind(a = a.fctr, x.tr), maxit = 500)
    pihat_mat_tr <- predict(pihat_fit, newdata = x.tr, type = "prob")
    pihat_mat <- predict(pihat_fit, newdata = x.te, type = "prob")
    pahat_tmp <- apply(pihat_mat_tr, 2, mean)
    
    for(u in 1:n.te) {
      pihat[u] <- pihat_mat[u, which(a.te[u] == colnames(pihat_mat))]
      pahat[u] <- pahat_tmp[which(a.te[u] == names(pahat_tmp))]
    }

    what <- pahat / pihat
    
    if(double_split) {
      s.tr <- sample(rep(1:3, ceiling(n.tr / 3))[1:n.tr])
      tr.idx1 <- s.tr == 1
      tr.idx2 <- s.tr == 2
      tr.idx3 <- s.tr == 3
    } else {
      tr.idx1 <- tr.idx2 <- tr.idx3 <- tr.idx
    }
 
    x.tr1 <- x[tr.idx1, , drop = FALSE]
    x.tr2 <- x[tr.idx2, , drop = FALSE]
    
    y.tr1 <- y[tr.idx1]
    y.tr2 <- y[tr.idx2]
    
    a.tr1 <- a[tr.idx1]
    a.tr2 <- a[tr.idx2]
    
    if(any(learner == "dr")) {
      
      tau.vals <- sort(unique(c(tau.l, tau.u)))
      qhat <- q(y.tr, a.tr, x.tr, new.x = x.te, new.a = a.te, 
                tau = tau.vals)
      
      qhat.mu <- q(y = y.tr1, a = a.tr1, x = x.tr1, new.x = x.tr2, 
                   new.a = a.tr2, tau = tau.vals)
      
      qhat.u <- qhat[, which(tau.vals >= 0.5), drop = FALSE]
      qhat.l <- qhat[,  which(tau.vals <= 0.5), drop = FALSE]
      qhat.mu.u <- qhat.mu[, which(tau.vals >= 0.5), drop = FALSE]
      qhat.mu.l <- qhat.mu[, which(tau.vals <= 0.5), drop = FALSE]
      qhat.l <- qhat.l[, order(ng:1), drop = FALSE]
      qhat.mu.l <- qhat.mu.l[, order(ng:1), drop = FALSE]
      
      if(double_split) {
        fit_v <- v(y = c(y.tr1, y.tr2), a = c(a.tr1, a.tr2), 
                   x = rbind(x.tr1, x.tr2))
      
      } else {
        fit_v <- v(y = c(y.tr), a = c(a.tr), x = x.tr)
        
      }
      
      vhat.te <- predict(fit_v, newdata = cbind(a = a.te, x.te))
      
      ax.te.vals <- expand.grid(a.vals, x.te)
      ax.u.vals <- expand.grid(1:length(a.te), 1:length(a.te))
      ax.u.vals <- ax.u.vals[ax.u.vals[, 1] != ax.u.vals[, 2], ]
      vhat.u.vals <- predict(fit_v, 
                             newdata = cbind(a = a.te[ax.u.vals[, 1]], 
                                             x.te[ax.u.vals[, 2], , 
                                                  drop = FALSE]))

      y.mu.u <- apply(qhat.mu.u, 2, function(x) I(y.tr2 <= x)) 
      y.mu.l <- apply(qhat.mu.l, 2, function(x) I(y.tr2 <= x)) 
      
      fit_mu.u <- fit_mu.l <- vector(mode = "list", length = ng)
      y.mu.u.pseudo <- t(t(y.mu.u) * y.tr2)
      y.mu.l.pseudo <- t(t(y.mu.l) * y.tr2)
      muhat.te.u <- muhat.te.l <- if.mu.mhat.u.a <- if.mu.mhat.l.a <- 
        if.mu.mhat.u.x <- if.mu.mhat.l.x <- matrix(NA, ncol = ng, nrow = n.te)
      mu.mhat.u <- mu.mhat.l <- matrix(NA, ncol = ng, nrow = nrow(ax.u.vals))
      
      mu.form <- paste0(paste0("y ~ ", paste0(disc.x, collapse = " + ")),
                        " + a + dmeduc + dfeduc + dmage + dfage + ",
                        " dfage_zero + nprevist + disllb + dlivord")
      
      for(u in 1:ng) {
        fit_mu.u[[u]] <- mu(y = y.mu.u.pseudo[, u], a = a.tr2, x = x.tr2)
        fit_mu.l[[u]] <- mu(y = y.mu.l.pseudo[, u], a = a.tr2, x = x.tr2)
      }
      for(u in 1:ng) {
        
        muhat.te.u[, u] <- predict(fit_mu.u[[u]], 
                                   newdata = cbind(a = a.te, x.te))
        muhat.te.l[, u] <- predict(fit_mu.l[[u]], 
                                   newdata = cbind(a = a.te, x.te))
          
        mu.mhat.u[, u] <- predict(fit_mu.u[[u]], 
                                  newdata = cbind(a = a.te[ax.u.vals[, 1]], 
                                                  x.te[ax.u.vals[, 2], , drop = FALSE]))
        mu.mhat.l[, u] <- predict(fit_mu.l[[u]], 
                                  newdata = cbind(a = a.te[ax.u.vals[, 1]], 
                                                  x.te[ax.u.vals[, 2], , drop = FALSE]))
      }
      
    }
    varphi.u <- varphi.l <- 
      if.moment.u <- if.moment.l <- matrix(NA, nrow = n.te, ncol = ng)
    pseudo.u <- pseudo.l <- matrix(NA, nrow = nrow(ax.u.vals), ncol = ng)
    
    for(j in 1:ng) {
      
      varphi.u[, j] <- c1.u[j] * what * 
        (qhat.u[, j] * (tau.u[j] - I(y.te <= qhat.u[, j])) + 
           y.te * I(y.te <= qhat.u[, j]) - muhat.te.u[, j]) +
        c2.u[j] * what * (y.te - vhat.te) 
      
      varphi.l[, j] <- c1.l[j] * what * 
        (qhat.l[, j] * (tau.l[j] - I(y.te <= qhat.l[, j])) + 
           y.te * I(y.te <= qhat.l[, j]) - muhat.te.l[, j]) +
        c2.l[j] * what * (y.te - vhat.te)
      
      pseudo.u[, j] <- c(varphi.u[ax.u.vals[, 1], j] + 
                           c1.u[j] * mu.mhat.u[, j] + c2.u[j] * vhat.u.vals)
      
      pseudo.l[, j] <- c(varphi.l[ax.u.vals[, 1], j] +
                           c1.l[j] * mu.mhat.l[, j] + c2.l[j] * vhat.u.vals)
    }
    a.reg <- a.te[ax.u.vals[, 1]]
    for(j in 1:ng) {
      for(alg in 1:length(learner)){
        if(ng > 1) {
          if(any(learner == "dr")) {
    
            est.u[[alg]][, , k, j] <- drl[[alg]](y = pseudo.u[, j], a = a.reg, new.a = a0)
            est.l[[alg]][, , k, j] <- drl[[alg]](y = pseudo.l[, j], a = a.reg, new.a = a0)
          
          }
        }
      } 
    }
    pseudo.u.folds[[k]] <- pseudo.u
    pseudo.l.folds[[k]] <- pseudo.l
    a.te.folds[[k]] <- a.te[ax.u.vals[, 1]]
  }
  ret.u <- sapply(1:length(learner), function(x) apply(est.u[[x]], c(1, 2, 4), mean),
                  simplify = FALSE)
  ret.l <- sapply(1:length(learner), function(x) apply(est.l[[x]], c(1, 2, 4), mean),
                  simplify = FALSE)
  ret <- list(est.u = ret.u, est.l = ret.l, fold.est.u = est.u, 
              fold.est.l = est.l, pseudo.u.folds = pseudo.u.folds, 
              pseudo.l.folds = pseudo.l.folds)
  return(ret)
}


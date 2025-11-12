# =========================================================
# MISE À JOUR : ajout du test manquant "Best Subset (Cp)"
# (Le reste de votre pipeline est inchangé. On active Cp
#  dans run_once() et on l’affiche dans les résumés du §5.)
# =========================================================

# =========================================================
# 0) SETUP DES PACKAGES
# =========================================================
pkgs <- c("MASS","leaps")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, function(p) require(p, character.only = TRUE)))

# =========================================================
# 1) GENERATION DES DONNEES (2 prédicteurs, corrélation contrôlée)
# =========================================================
gen_data <- function(n, beta0 = 1, beta1 = 2, beta2 = -1.5,
                     sigma2 = 1, rho = 0, n_test = n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  X <- MASS::mvrnorm(n = n + n_test, mu = c(0,0), Sigma = Sigma)
  colnames(X) <- c("X1","X2")
  eps <- rnorm(n + n_test, sd = sqrt(sigma2))
  y <- beta0 + beta1*X[,"X1"] + beta2*X[,"X2"] + eps
  
  idx_tr <- 1:n
  idx_te <- (n+1):(n+n_test)
  list(
    Xtr = X[idx_tr, , drop = FALSE],
    ytr = y[idx_tr],
    Xte = X[idx_te, , drop = FALSE],
    yte = y[idx_te],
    params = c(n = n, beta0 = beta0, beta1 = beta1, beta2 = beta2,
               sigma2 = sigma2, rho = rho)
  )
}

# =========================================================
# 2) METHODES DE SELECTION
# =========================================================

# 2.1) Sélection par test de Student (p < alpha)
student_select <- function(X, y, alpha = 0.05, intercept = TRUE, adjust = NULL) {
  df <- data.frame(y = y, X1 = X[,1], X2 = X[,2])
  form <- if (intercept) y ~ X1 + X2 else y ~ X1 + X2 - 1
  fit <- lm(form, data = df)
  sm  <- summary(fit)
  pvs <- sm$coefficients[, "Pr(>|t|)"]
  if (intercept) pX <- pvs[setdiff(names(pvs), "(Intercept)")] else pX <- pvs
  if (!is.null(adjust)) pX <- p.adjust(pX, method = adjust)
  sel <- names(pX)[pX < alpha]
  list(selected = sel, fit = fit, pvals = pX)
}

# 2.2) Recherche exhaustive via leaps::regsubsets + refit (AIC/BIC/Cp)
best_subset_select <- function(X, y, crit = c("AIC","BIC","Cp"), intercept = TRUE) {
  crit <- match.arg(crit)
  df <- data.frame(y = y, X1 = X[,1], X2 = X[,2])
  form <- if (intercept) y ~ X1 + X2 else y ~ X1 + X2 - 1
  
  rs <- regsubsets(form, data = df, nvmax = 2, method = "seqrep")
  ss <- summary(rs)
  
  n <- nrow(df)
  fits <- vector("list", nrow(ss$which))
  crit_vals <- rep(Inf, nrow(ss$which))
  
  fit_full <- lm(form, data = df)
  RSS_full <- sum(residuals(fit_full)^2)
  k_full   <- length(coef(fit_full))
  sigma2_full <- RSS_full / (n - k_full)
  
  for (m in seq_len(nrow(ss$which))) {
    inc <- ss$which[m, ]
    nm  <- names(inc)[inc]
    if (intercept) {
      nm <- setdiff(nm, "(Intercept)")
      f <- if (length(nm) == 0) y ~ 1 else reformulate(nm, response = "y")
    } else {
      f <- if (length(nm) == 0) y ~ -1 else reformulate(nm, response = "y", intercept = FALSE)
    }
    fitm <- lm(f, data = df)
    fits[[m]] <- fitm
    RSS <- sum(residuals(fitm)^2)
    k   <- length(coef(fitm))
    if (crit == "Cp")      crit_vals[m] <- RSS / sigma2_full - (n - 2*k)
    else if (crit == "AIC") crit_vals[m] <- n * log(RSS / n) + 2 * k
    else                    crit_vals[m] <- n * log(RSS / n) + log(n) * k
  }
  idx <- which.min(crit_vals)
  best <- fits[[idx]]
  sel  <- setdiff(names(coef(best)), "(Intercept)")
  list(selected = sel, fit = best, crit_value = crit_vals[idx])
}

# 2.3) Stepwise (AIC ou BIC)
stepwise_select <- function(X, y, criterion = c("AIC","BIC"), intercept = TRUE) {
  criterion <- match.arg(criterion)
  df <- data.frame(y = y, X1 = X[,1], X2 = X[,2])
  form_full <- if (intercept) y ~ X1 + X2 else y ~ X1 + X2 - 1
  fit0 <- lm(form_full, data = df)
  if (criterion == "AIC") {
    st <- MASS::stepAIC(fit0, direction = "both", trace = FALSE)
  } else {
    n  <- nrow(df)
    st <- step(fit0, direction = "both", trace = FALSE, k = log(n))
  }
  sel <- setdiff(names(coef(st)), "(Intercept)")
  list(selected = sel, fit = st)
}

# =========================================================
# 3) EVALUATION (RMSE test + taille du modèle)
# =========================================================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

eval_method <- function(name, fit_obj, Xte, yte) {
  if (is.null(fit_obj) || length(coef(fit_obj)) == 0L) {
    yhat <- rep(0, nrow(Xte))
  } else {
    nd <- data.frame(X1 = Xte[, 1], X2 = Xte[, 2])
    yhat <- as.numeric(predict(fit_obj, newdata = nd))
  }
  data.frame(
    method    = name,
    k         = if (is.null(fit_obj)) 0L else length(coef(fit_obj)),
    test_RMSE = rmse(yte, yhat),
    stringsAsFactors = FALSE
  )
}

# =========================================================
# 4) PIPELINE D’UN ESSAI (toutes méthodes) — ACTIVE Cp
# =========================================================
run_once <- function(n = 150, beta0 = 1, beta1 = 2, beta2 = -1.5,
                     sigma2 = 1, rho = 0, seed = 123, n_test = 150,
                     alpha = 0.05, adjust = NULL, intercept = TRUE) {
  dat <- gen_data(n, beta0, beta1, beta2, sigma2, rho, n_test, seed)
  
  # Student
  stu <- student_select(dat$Xtr, dat$ytr, alpha = alpha, intercept = intercept, adjust = adjust)
  f_stu <- if (length(stu$selected) == 0) {
    if (intercept) y ~ 1 else y ~ -1
  } else {
    reformulate(stu$selected, response = "y", intercept = intercept)
  }
  fit_stu <- lm(f_stu, data = data.frame(y = dat$ytr, dat$Xtr))
  
  # Best subset (AIC/BIC/Cp)  ← Cp ACTIVÉ
  bs_aic <- best_subset_select(dat$Xtr, dat$ytr, crit = "AIC", intercept = intercept)
  bs_bic <- best_subset_select(dat$Xtr, dat$ytr, crit = "BIC", intercept = intercept)
  bs_cp  <- best_subset_select(dat$Xtr, dat$ytr, crit = "Cp",  intercept = intercept)
  
  # Stepwise (AIC/BIC)
  sw_aic <- stepwise_select(dat$Xtr, dat$ytr, criterion = "AIC", intercept = intercept)
  sw_bic <- stepwise_select(dat$Xtr, dat$ytr, criterion = "BIC", intercept = intercept)
  
  # Tableau des résultats
  res <- rbind(
    eval_method("Student(alpha=0.05)", fit_stu,       dat$Xte, dat$yte),
    eval_method("BestSubset(AIC)",     bs_aic$fit,    dat$Xte, dat$yte),
    eval_method("BestSubset(BIC)",     bs_bic$fit,    dat$Xte, dat$yte),
    eval_method("BestSubset(Cp)",      bs_cp$fit,     dat$Xte, dat$yte),   # ← AJOUTÉ
    eval_method("Stepwise(AIC)",       sw_aic$fit,    dat$Xte, dat$yte),
    eval_method("Stepwise(BIC)",       sw_bic$fit,    dat$Xte, dat$yte)
  )
  
  sel_to_str <- function(s) if (length(s)==0) "NONE" else paste(s, collapse = "+")
  res$selected <- c(
    sel_to_str(stu$selected),
    sel_to_str(bs_aic$selected),
    sel_to_str(bs_bic$selected),
    sel_to_str(bs_cp$selected),     # ← AJOUTÉ
    sel_to_str(sw_aic$selected),
    sel_to_str(sw_bic$selected)
  )
  
  as.data.frame(cbind(res, t(as.data.frame(dat$params))), stringsAsFactors = FALSE)
}

# =========================================================
# 5) EXEMPLES RAPIDES (avec Cp dans les sorties)
# =========================================================
# a) Effet de n
res_n_small <- run_once(n = 80,  beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0, seed=1)
res_n_big   <- run_once(n = 400, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0, seed=1)
rbind(res_n_small, res_n_big)

# b) Effet du bruit (sigma2)
res_sig_low <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=0.25, rho=0, seed=2)
res_sig_hi  <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=2.0,  rho=0, seed=2)
rbind(res_sig_low, res_sig_hi)

# c) Effet de la taille d'effet (|beta2|)
res_b2_weak <- run_once(n=150, beta0=1, beta1=2, beta2=-0.4, sigma2=1, rho=0, seed=3)
res_b2_strg <- run_once(n=150, beta0=1, beta1=2, beta2=-2.5, sigma2=1, rho=0, seed=3)
rbind(res_b2_weak, res_b2_strg)

# c.2) Effet de la taille d'effet (|beta1|)
res_b2_weak <- run_once(n=150, beta0=1, beta1=0.25, beta2=2.5, sigma2=1, rho=0, seed=3)
res_b2_strg <- run_once(n=150, beta0=1, beta1=2, beta2=-2.5, sigma2=1, rho=0, seed=3)
rbind(res_b2_weak, res_b2_strg)

# d) Effet de la colinéarité (rho)
res_rho0   <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0.0, seed=4)
res_rho08  <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0.8, seed=4)
rbind(res_rho0, res_rho08)


# e) Effet de n avec (rho)
res_rho0   <- run_once(n=40, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0.8, seed=4)
res_rho08  <- run_once(n=400, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0.8, seed=4)
rbind(res_rho0, res_rho08)


# f) Effet du bruit (sigma2) avec (rho)
res_rho0   <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=0.25, rho=0.8, seed=4)
res_rho08  <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=2, rho=0.8, seed=4)
rbind(res_rho0, res_rho08)




# g) Effet de la taille d'effet (|beta1|) avec (rho)
res_rho0   <- run_once(n=150, beta0=1, beta1=0.25, beta2=-1.5, sigma2=1, rho=0.8, seed=4)
res_rho08  <- run_once(n=150, beta0=1, beta1=2, beta2=-1.5, sigma2=1, rho=0.8, seed=4)
rbind(res_rho0, res_rho08)


# h) Effet de la taille d'effet (|beta1|) avec (rho)
res_rho0   <- run_once(n=150, beta0=1, beta1=2, beta2=-0.4, sigma2=1, rho=0.8, seed=4)
res_rho08  <- run_once(n=150, beta0=1, beta1=2, beta2=-2.5, sigma2=1, rho=0.8, seed=4)
rbind(res_rho0, res_rho08)







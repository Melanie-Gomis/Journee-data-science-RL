## =========================================
## 0) Chargement des packages nécessaires
## =========================================
pkgs <- c("MASS","leaps")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, require, character.only = TRUE))

## =========================================
## 1) Générateur de données n x p
## =========================================
generate.lm <- function(n.train, p, p0, sigma2, n.test = 10*n.train, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- n.train + n.test
  train <- 1:n.train
  test  <- (n.train+1):n
  beta <- numeric(p)
  S.star <- sample(1:p, p0)
  beta[S.star] <- runif(p0, 1, 2) * sample(c(-1,1), p0, replace=TRUE)
  X <- matrix(rnorm(n*p), n, p)
  bruit <- rnorm(n, sd = sqrt(sigma2))  # sigma2 = variance du bruit
  y <- as.vector(X %*% beta + bruit)
  list(y = y, X = X, beta = beta, sigma2 = sigma2, train = train, test = test)
}

## =========================================
## 2) Fonctions utilitaires
## =========================================
.ensure_colnames <- function(X) {
  X <- as.data.frame(X, check.names = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  X
}
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

support_metrics <- function(pred, true) {
  TP <- sum(pred == 1 & true == 1)
  FP <- sum(pred == 1 & true == 0)
  FN <- sum(pred == 0 & true == 1)
  precision <- ifelse(TP + FP > 0, TP/(TP+FP), NA)
  rappel    <- ifelse(TP + FN > 0, TP/(TP+FN), NA)
  F1        <- ifelse(precision + rappel > 0, 2*precision*rappel/(precision+rappel), NA)
  c(precision = precision, rappel = rappel, F1 = F1)
}

## =========================================
## 3) Méthodes de sélection -> coefficients (p), zéros si non sélectionnés
## =========================================
betas_student <- function(X, y, alpha = 0.05) {
  X <- .ensure_colnames(X); vars <- colnames(X); p <- ncol(X)
  df <- data.frame(y = y, X)
  fit_full <- lm(y ~ ., data = df)
  pv <- summary(fit_full)$coefficients[-1, 4]
  sel_vars <- vars[which(pv < alpha)]
  betas <- setNames(rep(0, p), vars)
  if (!length(sel_vars)) return(betas)
  fit_sel <- lm(reformulate(sel_vars, response = "y"), data = df)
  bhat <- coef(fit_sel); bhat <- bhat[names(bhat) != "(Intercept)"]
  betas[names(bhat)] <- unname(bhat)
  betas
}

betas_stepwise <- function(X, y, criterion = c("AIC","BIC")) {
  X <- .ensure_colnames(X); vars <- colnames(X); p <- ncol(X)
  df <- data.frame(y = y, X)
  criterion <- match.arg(criterion)
  fit_full <- lm(y ~ ., data = df)
  fit_sw <- if (criterion == "AIC") {
    MASS::stepAIC(fit_full, direction = "both", trace = FALSE)
  } else {
    step(fit_full, direction = "both", k = log(nrow(df)), trace = FALSE)
  }
  bhat <- coef(fit_sw); bhat <- bhat[names(bhat) != "(Intercept)"]
  betas <- setNames(rep(0, p), vars)
  if (length(bhat)) betas[names(bhat)] <- unname(bhat)
  betas
}

betas_bestsubset <- function(X, y, criterion = c("AIC","BIC"), nvmax = 10, method = "seqrep") {
  X <- .ensure_colnames(X); vars <- colnames(X); p <- ncol(X)
  df <- data.frame(y = y, X)
  criterion <- match.arg(criterion)
  nvmax <- min(nvmax, p)
  rs <- leaps::regsubsets(y ~ ., data = df, nvmax = nvmax, method = method)
  ss <- summary(rs)
  n <- nrow(df)
  best_val <- Inf; best_vars <- character(0)
  for (m in seq_len(nrow(ss$which))) {
    inc <- ss$which[m, ]
    cand <- setdiff(names(inc)[inc], "(Intercept)")
    if (!length(cand)) next  # éviter le modèle vide
    fitm <- lm(reformulate(cand, response = "y"), data = df)
    RSS <- sum(residuals(fitm)^2); k <- length(coef(fitm))
    val <- if (criterion == "AIC") n*log(RSS/n) + 2*k else n*log(RSS/n) + log(n)*k
    if (val < best_val) { best_val <- val; best_vars <- cand }
  }
  betas <- setNames(rep(0, p), vars)
  if (!length(best_vars)) return(betas)
  fit_best <- lm(reformulate(best_vars, response = "y"), data = df)
  bhat <- coef(fit_best); bhat <- bhat[names(bhat) != "(Intercept)"]
  if (length(bhat)) betas[names(bhat)] <- unname(bhat)
  betas
}

## =========================================
## 4) Pipeline complet : renvoie tous les objets avec tailles n et p
## =========================================
run_selection_pipeline <- function(n.train = 100, p = 50, p0 = 10, sigma2 = 1,
                                   n.test = 100, nvmax = 10, seed = 42,
                                   standardiser = TRUE) {
  data <- generate.lm(n.train, p, p0, sigma2, n.test, seed)
  Xtr <- data$X[data$train, , drop = FALSE]
  ytr <- data$y[data$train]
  Xte <- data$X[data$test,  , drop = FALSE]
  yte <- data$y[data$test]
  colnames(Xtr) <- paste0("X", seq_len(ncol(Xtr)))
  colnames(Xte) <- colnames(Xtr)
  
  # (optionnel) standardiser les colonnes de X pour plus de stabilité
  if (standardiser) {
    mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv == 0] <- 1
    Xtr <- scale(Xtr, center = mu, scale = sdv)
    Xte <- scale(Xte, center = mu, scale = sdv)
  }
  
  # Coefficients (p x 5)
  b_student   <- betas_student(Xtr, ytr, alpha = 0.05)
  b_step_aic  <- betas_stepwise(Xtr, ytr, "AIC")
  b_step_bic  <- betas_stepwise(Xtr, ytr, "BIC")
  b_best_aic  <- betas_bestsubset(Xtr, ytr, "AIC", nvmax = nvmax)
  b_best_bic  <- betas_bestsubset(Xtr, ytr, "BIC", nvmax = nvmax)
  
  B <- cbind(Student = b_student,
             StepAIC = b_step_aic,
             StepBIC = b_step_bic,
             BestAIC = b_best_aic,
             BestBIC = b_best_bic)        # dimension p x 5
  
  # Sélection binaire (p x 5)
  S <- (B != 0) * 1L
  rownames(S) <- rownames(B) <- colnames(Xtr)
  
  # Prédictions sur le jeu de test (n.test x 5)
  yhat <- matrix(NA_real_, nrow = n.test, ncol = 5,
                 dimnames = list(NULL, c("Student","StepAIC","StepBIC","BestAIC","BestBIC")))
  # Fonction auxiliaire pour prédire à partir des coefficients (sans intercept)
  predict_with_betas <- function(X, betas) as.vector(X %*% betas)
  
  yhat[, "Student"] <- predict_with_betas(Xte, b_student)
  yhat[, "StepAIC"] <- predict_with_betas(Xte, b_step_aic)
  yhat[, "StepBIC"] <- predict_with_betas(Xte, b_step_bic)
  yhat[, "BestAIC"] <- predict_with_betas(Xte, b_best_aic)
  yhat[, "BestBIC"] <- predict_with_betas(Xte, b_best_bic)
  
  # RMSE sur le jeu de test pour chaque méthode
  rmse_vec <- apply(yhat, 2, function(col) rmse(yte, col))
  
  # Qualité de récupération du support (comparée au vrai beta)
  support_true <- as.integer(data$beta != 0)
  metrics <- rbind(
    Student = support_metrics(S[, "Student"], support_true),
    StepAIC = support_metrics(S[, "StepAIC"], support_true),
    StepBIC = support_metrics(S[, "StepBIC"], support_true),
    BestAIC = support_metrics(S[, "BestAIC"], support_true),
    BestBIC = support_metrics(S[, "BestBIC"], support_true)
  )
  
  list(
    # tailles explicites
    Xtr = Xtr,               # n.train x p
    ytr = ytr,               # n.train
    Xte = Xte,               # n.test x p
    yte = yte,               # n.test
    beta_true = data$beta,   # p
    support_true = support_true, # p
    B = B,                   # p x 5 (coefficients ; 0 si non sélectionnés)
    S = S,                   # p x 5 (sélection binaire)
    Yhat_test = yhat,        # n.test x 5
    RMSE_test = rmse_vec,    # 5 valeurs
    support_metrics = metrics # 5 x 3
  )
}

## =========================================
## 5) Exemple d’utilisation : tailles n, p contrôlées
## =========================================
res <- run_selection_pipeline(
  n.train = 100, p = 50, p0 = 10, sigma2 = 1,
  n.test = 100, nvmax = 12, seed = 123, standardiser = TRUE
)

# Dimensions attendues
dim(res$Xtr)       # n.train x p
length(res$ytr)    # n.train
dim(res$Xte)       # n.test x p
length(res$yte)    # n.test
dim(res$B)         # p x 5
dim(res$S)         # p x 5
dim(res$Yhat_test) # n.test x 5

# Aperçu rapide
head(res$B)

matrice <- res$B

res$RMSE_test
round(res$support_metrics, 3)


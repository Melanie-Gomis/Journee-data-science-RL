source("scripts/00_setup.R")
source("scripts/01_generate_data.R")
source("scripts/04_selection_des_variables.R")
# saveRDS(cas_a_explorer,"data/cas_a_explorer.rds")

cat("\nPerformance des modèles")

perf_safe <- function(X_test, y_test, beta, beta.star) {
  # -- forcer les bons types
  X_test   <- as.matrix(X_test)
  y_test   <- as.numeric(y_test)
  beta     <- as.numeric(beta)
  beta.star<- as.numeric(beta.star)
  
  pX <- ncol(X_test)
  pB <- length(beta)
  pT <- length(beta.star)
  
  # -- sécurités de dimension
  if (pB != pX) {
    # on tente de « réparer » : tronquer ou zéro-remplir
    if (pB > pX) {
      beta <- beta[seq_len(pX)]
    } else {
      beta <- c(beta, rep(0, pX - pB))
    }
  }
  if (length(beta.star) != pX) {
    # si le vrai beta n'a pas la bonne taille, on ajuste pareil
    if (length(beta.star) > pX) {
      beta.star <- beta.star[seq_len(pX)]
    } else {
      beta.star <- c(beta.star, rep(0, pX - length(beta.star)))
    }
  }
  
  # -- métriques de support
  nzero     <- which(beta != 0)
  zero      <- which(beta == 0)
  true.nzero<- which(beta.star != 0)
  true.zero <- which(beta.star == 0)
  
  TP <- sum(nzero %in% true.nzero)
  TN <- sum(zero  %in% true.zero)
  FP <- sum(nzero %in% true.zero)
  FN <- sum(zero  %in% true.nzero)
  
  recall      <- if ((TP + FN) > 0) TP/(TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN/(TN + FP) else NA_real_
  precision   <- if ((TP + FP) > 0) TP/(TP + FP) else NA_real_
  
  rmse <- sqrt(mean((beta - beta.star)^2, na.rm = TRUE))
  
  # -- prédiction test (conformable garanti)
  yhat_test <- as.numeric(X_test %*% beta)
  rerr <- sqrt(mean((y_test - yhat_test)^2))
  
  res <- round(c(precision, recall, specificity, rmse, rerr), 4)
  res[is.nan(res)] <- 0
  names(res) <- c("precision","recall","specificity","rmse","prediction")
  res
}

# 1) calculer la list-col des performances
cas_a_explorer <- cas_a_explorer %>%
  mutate(
    perf = pmap(
      list(data, beta_final),
      function(dat, beta_hat) {
        Xte <- dat$X[dat$test, , drop = FALSE]
        yte <- dat$y[dat$test]
        perf_safe(X_test = Xte, y_test = yte, beta = beta_hat, beta.star = dat$beta)
      }
    )
  )

# 2) extraire la métrique demandée par 'comparaison'
lib_to_metric <- c(
  "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train" = "rmse",
  "Erreur moyenne de prédiction calculée sur l’ensemble test"     = "prediction",
  "Précision du support estimé"                                   = "precision",
  "Sensibilité du support estimé"                                 = "recall",
  "Spécificité du support estimé"                                 = "specificity"
)


cas_a_explorer <- cas_a_explorer %>%
  mutate(
    metric_name = recode(comparaison, !!!lib_to_metric),
    score = map2_dbl(perf, metric_name, ~ .x[[.y]])
  ) 

cas_a_explorer_simpli <- cas_a_explorer %>%
  dplyr::select(-data,-beta_final,metric_name)


saveRDS(cas_a_explorer_simpli,"data/cas_a_explorer_simpli.rds")

source("scripts/00_setup.R")
source("scripts/01_generate_data.R")
source("scripts/02_fonction_test.R")

library(furrr)

# Petit routeur vers la bonne fonction de sélection
.fit_selector <- function(selection, X, Y) {
  if (selection == "Test de Student de signification des coefficients") {
    return(testStudent(X, Y))
  }
  if (selection == "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)") {
    return(bestsubset(X, Y))
  }
  if (selection == "Algorithme de sélection pas à pas (stepwise)") {
    return(stepwise(X, Y))
  }
  stop("Méthode de sélection inconnue: ", selection)
}

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

lib_to_metric <- c(
  "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train" = "rmse",
  "Erreur moyenne de prédiction calculée sur l’ensemble test"     = "prediction",
  "Précision du support estimé"                                   = "precision",
  "Sensibilité du support estimé"                                 = "recall",
  "Spécificité du support estimé"                                 = "specificity"
)


cas_a_explorer <- expand.grid(
  type_donnees = c("Indépendantes","Dépendance longitudinale"),
  p = 10,
  n = c(50, 100),                        # ici: interprété comme n.train
  sigma2 = c(1, 200),
  correlation_longitudinale = c(0.3, 0.8),
  # selection = c(
  #   "Test de Student de signification des coefficients",
  #   "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
  #   "Algorithme de sélection pas à pas (stepwise)"
  # ),
  # comparaison = c(
  #   "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train",
  #   "Erreur moyenne de prédiction calculée sur l’ensemble test",
  #   "Précision du support estimé",
  #   "Sensibilité du support estimé",
  #   "Spécificité du support estimé"
  # ),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>% 
  mutate(rapport_observation_variable = n / p) 


# Effet taille ----
no_cores <- availableCores() - 1
plan(multicore, workers = no_cores)

cas_a_explorer_taille <- expand.grid(
  type_donnees = c("Indépendantes"),
  p = 10,
  n = c(seq(50,100,by=50),seq(200,500,by=100),seq(750,1000,by=250)),# ici: interprété comme n.train
  sigma2 = 1,
  # selection = c(
  #   "Test de Student de signification des coefficients",
  #   "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
  #   "Algorithme de sélection pas à pas (stepwise)"
  # ),
  # comparaison = c(
  #   "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train",
  #   "Erreur moyenne de prédiction calculée sur l’ensemble test",
  #   "Précision du support estimé",
  #   "Sensibilité du support estimé",
  #   "Spécificité du support estimé"
  # ),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>% 
  mutate(rapport_observation_variable = n / p) 


cas_a_explorer_taille_simpli <-cas_a_explorer_taille %>%
  mutate(
    correlation_longitudinale=0,
    perf =  future_pmap(
      # pmap(
      pick(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable),
      function(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable) {
        p0 <- max(1, floor(p / 5))   # ex. 20% de variables actives
        if (type_donnees == "Indépendantes") {
          data_temp <- generate.lm(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rapport_observation_variable = rapport_observation_variable
          )
        } else {
          data_temp <- generate.lm.longitudinal(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rho = correlation_longitudinale,
            rapport_observation_variable = rapport_observation_variable
          )
        }
        
        Xtr <- data_temp$X[data_temp$train, , drop = FALSE]
        Ytr <- data_temp$y[data_temp$train]
        
        selection = c(
          "Test de Student de signification des coefficients",
          "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
          "Algorithme de sélection pas à pas (stepwise)"
        )
        list_res <- map_df(selection,function(selection){
          temp_beta <- .fit_selector(selection, Xtr, Ytr)
          
          Xte <- data_temp$X[data_temp$test, , drop = FALSE]
          yte <- data_temp$y[data_temp$test]
          res <- perf_safe(X_test = Xte, y_test = yte, beta = temp_beta, beta.star = data_temp$beta)
          res_tbl <- tibble(
            score = res,
            metric=c("precision","recall","specificity","rmse","prediction")
          )
          return(res_tbl)
        })
        
        list_res <- list_res %>% 
          mutate(selection=rep(selection,each=5))
        return(list_res)  
      })
  ) %>% 
  unnest()  



graph_cas_a_explorer_taille <- cas_a_explorer_taille_simpli %>% 
  mutate(
    metric_name = factor(metric, levels=lib_to_metric,labels= str_wrap(names(lib_to_metric),width=25))
    # score = map2_dbl(perf, metric_name, ~ .x[[.y]])
  ) %>% 
  filter(metric%in% c("recall","rmse","prediction")) %>%
  ggplot(aes(x = n,y=score,color=selection))+
  # geom_line()+
  geom_smooth(se = F)+
  facet_wrap("metric_name" ,scales="free")+
  theme(
    legend.position = "bottom"
    )+
  labs(
    title = "Comparaison des methodes de selection de variables\nen fonction de l'évolution de la taille"
  )



# Effet correlation longitudinale ----
no_cores <- availableCores() - 1
plan(multicore, workers = no_cores)

cas_a_explorer_dependance_longitudinale <- expand.grid(
  type_donnees = c("Dépendance longitudinale"),
  p = 10,
  n = 100,# ici: interprété comme n.train
  sigma2 = 1,
  correlation_longitudinale = seq(0,1,0.1),
  # selection = c(
  #   "Test de Student de signification des coefficients",
  #   "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
  #   "Algorithme de sélection pas à pas (stepwise)"
  # ),
  # comparaison = c(
  #   "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train",
  #   "Erreur moyenne de prédiction calculée sur l’ensemble test",
  #   "Précision du support estimé",
  #   "Sensibilité du support estimé",
  #   "Spécificité du support estimé"
  # ),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>% 
  mutate(rapport_observation_variable = n / p) 


cas_a_explorer_dependance_longitudinale_simpli <-cas_a_explorer_dependance_longitudinale %>%
  mutate(
    perf =  future_pmap(
      # pmap(
      pick(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable),
      function(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable) {
        p0 <- max(1, floor(p / 5))   # ex. 20% de variables actives
        if (type_donnees == "Indépendantes") {
          data_temp <- generate.lm(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rapport_observation_variable = rapport_observation_variable
          )
        } else {
          data_temp <- generate.lm.longitudinal(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rho = correlation_longitudinale,
            rapport_observation_variable = rapport_observation_variable
          )
        }
        
        Xtr <- data_temp$X[data_temp$train, , drop = FALSE]
        Ytr <- data_temp$y[data_temp$train]
        
        selection = c(
          "Test de Student de signification des coefficients",
          "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
          "Algorithme de sélection pas à pas (stepwise)"
        )
        list_res <- map_df(selection,function(selection){
          temp_beta <- .fit_selector(selection, Xtr, Ytr)
          
          Xte <- data_temp$X[data_temp$test, , drop = FALSE]
          yte <- data_temp$y[data_temp$test]
          res <- perf_safe(X_test = Xte, y_test = yte, beta = temp_beta, beta.star = data_temp$beta)
          res_tbl <- tibble(
            score = res,
            metric=c("precision","recall","specificity","rmse","prediction")
          )
          return(res_tbl)
        })
        
        list_res <- list_res %>% 
          mutate(selection=rep(selection,each=5))
        return(list_res)  
      })
  ) %>% 
  unnest()  



graph_cas_a_explorer_dependance_longitudinale <- cas_a_explorer_dependance_longitudinale_simpli %>% 
  mutate(
    metric_name = factor(metric, levels=lib_to_metric,labels= str_wrap(names(lib_to_metric),width=25))
  ) %>% 
  filter(metric%in% c("recall","rmse","prediction")) %>% 
  ggplot(aes(x =correlation_longitudinale,y=score,color=selection))+
  # geom_line()+
  geom_smooth(se = F)+
  facet_wrap("metric_name" ,scales="free")+
  theme(
    legend.position = "bottom"
  )+
  labs(
    title = "Comparaison des methodes de selection de variables\nen fonction de l'évolution de la correlation longitudinale"
  )


# Effet sigma2 ----
no_cores <- availableCores() - 1
plan(multicore, workers = no_cores)

cas_a_explorer_sigma <- expand.grid(
  type_donnees = c("Indépendantes"),
  p = 10,
  n = 100,# ici: interprété comme n.train
  sigma2 = unique(c(seq(1,10,by=1),seq(20,50,by=10),seq(50,100,by=50),seq(200,500,by=100),seq(750,1000,by=250))),
  # selection = c(
  #   "Test de Student de signification des coefficients",
  #   "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
  #   "Algorithme de sélection pas à pas (stepwise)"
  # ),
  # comparaison = c(
  #   "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train",
  #   "Erreur moyenne de prédiction calculée sur l’ensemble test",
  #   "Précision du support estimé",
  #   "Sensibilité du support estimé",
  #   "Spécificité du support estimé"
  # ),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>% 
  mutate(rapport_observation_variable = n / p) 


cas_a_explorer_sigma_simpli <-cas_a_explorer_sigma %>%
  mutate(
    correlation_longitudinale=0,
    perf =  future_pmap(
      # pmap(
      pick(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable),
      function(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable) {
        p0 <- max(1, floor(p / 5))   # ex. 20% de variables actives
        if (type_donnees == "Indépendantes") {
          data_temp <- generate.lm(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rapport_observation_variable = rapport_observation_variable
          )
        } else {
          data_temp <- generate.lm.longitudinal(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rho = correlation_longitudinale,
            rapport_observation_variable = rapport_observation_variable
          )
        }
        
        Xtr <- data_temp$X[data_temp$train, , drop = FALSE]
        Ytr <- data_temp$y[data_temp$train]
        
        selection = c(
          "Test de Student de signification des coefficients",
          "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
          "Algorithme de sélection pas à pas (stepwise)"
        )
        list_res <- map_df(selection,function(selection){
          temp_beta <- .fit_selector(selection, Xtr, Ytr)
          
          Xte <- data_temp$X[data_temp$test, , drop = FALSE]
          yte <- data_temp$y[data_temp$test]
          res <- perf_safe(X_test = Xte, y_test = yte, beta = temp_beta, beta.star = data_temp$beta)
          res_tbl <- tibble(
            score = res,
            metric=c("precision","recall","specificity","rmse","prediction")
          )
          return(res_tbl)
        })
        
        list_res <- list_res %>% 
          mutate(selection=rep(selection,each=5))
        return(list_res)  
      })
  ) %>% 
  unnest()  



graph_cas_a_explorer_sigma <- cas_a_explorer_sigma_simpli %>% 
  mutate(
    metric_name = factor(metric, levels=lib_to_metric,labels= str_wrap(names(lib_to_metric),width=25))
  ) %>% 
  filter(metric%in% c("recall","rmse","prediction")) %>%
  ggplot(aes(x = sigma2,y=score,color=selection))+
  # geom_line()+
  geom_smooth(se = F)+
  facet_wrap("metric_name" ,scales="free")+
  theme(
    legend.position = "bottom"
  )+
  labs(
    title = "Comparaison des methodes de selection de variables\nen fonction de l'évolution de sigma2"
  )




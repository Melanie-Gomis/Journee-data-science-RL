# --- Chargement de la bibliothèque ---
library(MASS)
library(leaps)
cat("\nFonction selection de modèles")

stepwise <- function(X, Y) {
  
  # 1. Préparation des données
  # S'assurer que X est une matrice et obtenir le nombre de prédicteurs
  if (!is.matrix(X)) X <- as.matrix(X)
  p <- ncol(X)
  
  # Assigner des noms aux colonnes de X si elles n'en ont pas
  # (essentiel pour que lm() et stepAIC() fonctionnent correctement)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  # Combiner en un data.frame pour les fonctions lm()
  dat <- as.data.frame(cbind(Y, X))
  
  
  # 2. Définir les modèles
  # Modèle complet (point de départ)
  full_model <- lm(Y ~ ., data = dat)
  # Modèle nul (limite inférieure de la recherche)
  null_model <- lm(Y ~ 1, data = dat)
  
  
  # 3. Exécuter la sélection Stepwise (AIC)
  # direction = "both" effectue une sélection avant et arrière
  # k = 2 correspond au critère AIC (c'est la valeur par défaut)
  # trace = FALSE pour ne pas afficher les étapes dans la console
  step_model <- stepAIC(
    full_model,
    scope = list(lower = null_model, upper = full_model),
    direction = "both",
    k = 2, # k=2 pour AIC
    trace = FALSE
  )
  
  
  # 4. Extraire les coefficients du modèle final
  selected_coeffs <- coef(step_model)
  
  # Retirer l'intercept, car il ne fait pas partie des prédicteurs X
  selected_coeffs <- selected_coeffs[names(selected_coeffs) != "(Intercept)"]
  
  
  # 5. Créer le vecteur beta final (de taille p)
  # Obtenir les noms de TOUS les prédicteurs initiaux
  all_predictor_names <- colnames(X)
  
  # Créer un vecteur de 0 de la bonne longueur (p)
  beta_final <- rep(0, p)
  # Assigner les noms de tous les prédicteurs à ce vecteur
  names(beta_final) <- all_predictor_names
  
  # 6. Remplir le vecteur avec les coefficients sélectionnés
  # En utilisant les noms, on s'assure que les coefficients
  # sont assignés aux bonnes variables.
  beta_final[names(selected_coeffs)] <- selected_coeffs
  
  return(beta_final)
}


bestsubset <- function(X, Y) {
  
  # 1. Préparation des données
  if (!is.matrix(X)) X <- as.matrix(X)
  p <- ncol(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  # regsubsets() fonctionne mieux avec un data.frame
  dat <- as.data.frame(cbind(Y, X))
  all_predictor_names <- colnames(X)
  
  # 2. Exécuter la recherche exhaustive (Best Subset)
  # (Comme mentionné dans la section 1 de votre document)
  best_subset_model <- regsubsets(
    Y ~ ., 
    data = dat, 
    nvmax = p, 
    method = "exhaustive"
  )
  
  bs_summary <- summary(best_subset_model)
  
  # 3. Trouver le meilleur modèle basé sur le BIC
  # (Critère AIC/BIC mentionné dans la section 2.4 de votre document)
  best_model_index_bic <- which.min(bs_summary$bic)
  
  # 4. Obtenir les variables de ce modèle
  selected_vars_logical <- bs_summary$outmat[best_model_index_bic, ] == "*"
  selected_vars_names <- names(selected_vars_logical[selected_vars_logical])
  
  # Retirer l'intercept
  selected_vars_names <- selected_vars_names[selected_vars_names != "(Intercept)"]
  
  # 5. Créer le vecteur beta final (de taille p)
  beta_final <- rep(0, p)
  names(beta_final) <- all_predictor_names
  
  # 6. Remplir le vecteur
  if (length(selected_vars_names) > 0) {
    
    # Il faut ré-estimer le modèle linéaire OLS (moindres carrés ordinaires)
    # avec SEULEMENT les variables sélectionnées pour obtenir les coefficients.
    
    # Créer la formule pour le modèle final
    formula_best <- as.formula(
      paste("Y ~", paste(selected_vars_names, collapse = " + "))
    )
    
    # Ajuster le modèle OLS final
    final_model <- lm(formula_best, data = dat)
    
    # Extraire les coefficients (en ignorant l'intercept)
    final_coeffs <- coef(final_model)
    final_coeffs <- final_coeffs[names(final_coeffs) != "(Intercept)"]
    
    # Assigner les coefficients au vecteur final
    beta_final[names(final_coeffs)] <- final_coeffs
  }
  
  return(beta_final)
}


testStudent <- function(X, Y, p_value_threshold = 0.05) {
  
  # 1. Préparation des données
  if (!is.matrix(X)) X <- as.matrix(X)
  p <- ncol(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", 1:p)
  }
  
  dat <- as.data.frame(cbind(Y, X))
  all_predictor_names <- colnames(X)
  
  # 2. Ajuster le modèle complet pour obtenir les p-valeurs
  full_model <- lm(Y ~ ., data = dat)
  model_summary <- summary(full_model)
  
  # 3. Extraire les p-valeurs (coefficients[-1,...] pour ignorer l'intercept)
  p_values <- model_summary$coefficients[-1, "Pr(>|t|)"]
  
  # 4. Sélectionner les variables significatives
  selected_vars_names <- names(p_values[p_values < p_value_threshold])
  
  # 5. Créer le vecteur beta final (de taille p)
  beta_final <- rep(0, p)
  names(beta_final) <- all_predictor_names
  
  # 6. Remplir le vecteur
  # Si au moins une variable a été sélectionnée...
  if (length(selected_vars_names) > 0) {
    
    # ...on ré-ajuste un modèle OLS avec SEULEMENT ces variables
    # pour obtenir les bonnes valeurs de coefficients.
    
    # Créer la formule pour le modèle final
    formula_final <- as.formula(
      paste("Y ~", paste(selected_vars_names, collapse = " + "))
    )
    
    # Ajuster le modèle OLS final
    final_model <- lm(formula_final, data = dat)
    
    # Extraire les coefficients (en ignorant l'intercept)
    final_coeffs <- coef(final_model)
    final_coeffs <- final_coeffs[names(final_coeffs) != "(Intercept)"]
    
    # Assigner les coefficients au vecteur final
    beta_final[names(final_coeffs)] <- final_coeffs
  }
  
  return(beta_final)
}





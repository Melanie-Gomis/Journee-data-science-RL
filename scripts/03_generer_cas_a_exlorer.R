cat("\nFonctions générer cas à explorer")
cas_a_explorer <- expand.grid(
  type_donnees = c("Indépendantes","Dépendance longitudinale"),
  p = 10,
  n = c(50, 1000),                        # ici: interprété comme n.train
  sigma2 = c(1, 200),
  correlation_longitudinale = c(0.3, 0.8),
  selection = c(
    "Test de Student de signification des coefficients",
    "Recherche exhaustive (Best subset avec l’algorithme Leaps and bound)",
    "Algorithme de sélection pas à pas (stepwise)"
  ),
  comparaison = c(
    "Erreur quadratique moyenne de β^ (RMSE) sur l’ensemble train",
    "Erreur moyenne de prédiction calculée sur l’ensemble test",
    "Précision du support estimé",
    "Sensibilité du support estimé",
    "Spécificité du support estimé"
  ),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  as_tibble() %>% 
  mutate(rapport_observation_variable = n / p) 

dataset_a_explorer <- cas_a_explorer %>% 
  distinct(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable)

dataset_a_explorer <- dataset_a_explorer %>%
  # on crée une list-col 'data' en appelant le bon générateur par ligne
  mutate(
    data = pmap(
      pick(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable),
      function(type_donnees, p, n, sigma2, correlation_longitudinale, rapport_observation_variable) {
        p0 <- max(1, floor(p / 5))   # ex. 20% de variables actives
        if (type_donnees == "Indépendantes") {
          generate.lm(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rapport_observation_variable = rapport_observation_variable
          )
        } else {
          generate.lm.longitudinal(
            n.train = n, p = p, p0 = p0, sigma2 = sigma2,
            rho = correlation_longitudinale,
            rapport_observation_variable = rapport_observation_variable
          )
        }
      }
    )
  )


cas_a_explorer <- cas_a_explorer %>% 
  left_join(
    dataset_a_explorer,
    by=c("type_donnees", "p", "n", "sigma2", "correlation_longitudinale", "rapport_observation_variable")
  )

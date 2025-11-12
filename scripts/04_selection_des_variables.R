source("scripts/02_fonction_test.R")
source("scripts/03_generer_cas_a_exlorer.R")
cat("\nSelection des variables")

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


# Ajout de la colonne beta_final
 

cas_a_explorer <- cas_a_explorer %>%
  mutate(
    beta_final = pmap(
      list(selection, data),
      function(selection, dat) {
        # Extraction du train
        Xtr <- dat$X[dat$train, , drop = FALSE]
        Ytr <- dat$y[dat$train]
        # Application de la bonne méthode
        # .fit_selector(selection, Xtr, Ytr)
      }
    )
  )

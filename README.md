# Demi-journée Data Science : Evaluation des méthodes de sélection de variable (Protocole de simulation)


## Objectifs
- Problématique :

Un problème essentiel de la régression linéaire multiple est le choix des variables explicatives à conserver.
Sélectionner un sous-ensemble de variables permet notamment de réaliser un compromis entre biais et variance.

Cela représente également une solution intéressante en cas de multi-colinéarité.

Le choix de modèle dépend des objectifs de la régression (description des données, estimation des paramètres, prévision de nouvelles valeurs) et de la connaissance des données.

## Organisation du dépôt
```
Journee-data-science-RL/
├── data/                    # jeux de données (non versionnés si sensibles)
├── scripts/                 # code R par étapes
├── rapport.Rmd               # rapport reproductible
├── README.md
└── .gitignore
```

## Scripts
- `scripts/00_setup.R` : charge/installe les packages, fixe le seed (123), options knitr.
- `scripts/01_generer.R` : Générer les données.
- `scripts/02_selection.R` : la recherche exhaustive best subsets avec critère Cp, AIC et BIC et la régression stepwise avec critère AIC et BIC.
- `scripts/03_comparaison.R` : RMSE, $\hat{\varepsilon}$, précisions, sensibilité, spécificité.
- `scripts/04_simulation.R` :  Planning de simulations


## Exécution rapide
1) Ouvrir `report.Rmd` dans RStudio → Knit.
2) (option) Lancer les scripts dans l’ordre depuis `scripts/`.

## Auteurs
- Mélanie: (rôle) —
- Cédric: (rôle) —
- Maximiliano (rôle) —
- Enzo (rôle) —
- Arnaud (rôle) —

library(mvtnorm)

set.seed(123)  # reproductibilité
cat("\nFonctions generate data")
generate.lm <- function(n.train, p, p0, sigma2,rapport_observation_variable=10 ,n.test=rapport_observation_variable*n.train) {
  ## génération ensemble test/apprentissage
  ## Rq : Le test peut être grand ici, car il n'impacte pas la taille de l'échantillon d'apprentissage
  ## (Plus le test sera grand, plus les modèles seront évalués avec précision.)
  
  n <- n.train + n.test
  train <- 1:n.train
  test <- (n.train+1):n
  
  ## vecteur des vrais coefficients
  beta <- numeric(p)
  S.star <- sample(1:p, p0) ## vrai support S*
  beta[S.star] <- runif(p0, 1, 2)*sample(c(-1,1),p0, replace=T)
  ## bruit blanc gaussien de variance sigma2
  noise <- rnorm(n, sd=sigma2)
  ## prédicteurs sans structure particulière
  X <- matrix(rnorm(n*p), n, p)
  ## vecteur d'observation des réponses
  y <- ## vecteur des vrais coefficients
    beta <- numeric(p)
  S.star <- sample(1:p, p0) ## vrai support S*
  beta[S.star] <- runif(p0, 1, 2)*sample(c(-1,1),p0, replace=T)
  ## bruit blanc gaussien de variance sigma2
  noise <- rnorm(n, sd=sigma2)
  ## prédicteurs sans structure particulière
  X <- matrix(rnorm(n*p), n, p)
  ## vecteur d'observation des réponses
  y <-  X%*%beta + noise
  
  list(y=y, X=X, beta=beta, sigma2=sigma2, train=train, test=test)
}


generate.lm.longitudinal <- function(n.train, p, p0, sigma2,
                                     Ttime = 5,         # nombre de temps par sujet
                                     rho = 0.6,         # corrélation AR(1)
                                     tau2 = 0.5,        # variance effet aléatoire d'intercept
                                     rapport_observation_variable = 10,
                                     n.test = rapport_observation_variable * n.train,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- dimensions
  n.subj <- n.train + n.test          # nb de sujets
  N      <- n.subj * Ttime            # nb total d'observations (lignes)
  
  # --- index train/test (sur les LIGNES, comme dans votre version plate)
  id    <- rep(1:n.subj, each = Ttime)
  time  <- rep(0:(Ttime-1), times = n.subj)
  train.subj <- 1:n.train
  test.subj  <- (n.train + 1):n.subj
  train <- which(id %in% train.subj)
  test  <- which(id %in% test.subj)
  
  # --- vrais coefficients et support
  beta   <- numeric(p)
  S.star <- sample(1:p, p0)
  beta[S.star] <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  
  # --- design : prédicteurs sans structure particulière (pas d'intercept fixe)
  X <- matrix(rnorm(N * p), nrow = N, ncol = p)
  
  # --- structure de dépendance longitudinale
  # Erreurs résiduelles AR(1) intra-sujet
  SigmaT <- sigma2 * toeplitz(rho^(0:(Ttime-1)))
  # Effet aléatoire d'intercept par sujet (b_i ~ N(0, tau2))
  
  # --- génération de y
  y <- numeric(N)
  for (i in 1:n.subj) {
    rows_i <- ((i - 1) * Ttime + 1):(i * Ttime)
    Xi     <- X[rows_i, , drop = FALSE]
    
    b_i    <- rnorm(1, mean = 0, sd = sqrt(tau2))
    eps_i  <- as.numeric(rmvnorm(1, mean = rep(0, Ttime), sigma = SigmaT))
    
    y[rows_i] <- as.numeric(Xi %*% beta) + b_i + eps_i
  }
  
  # --- sortie (mêmes éléments que votre fonction + id/time utiles)

  list(
    y = y,
    X = X,
    beta = beta,
    sigma2 = sigma2,
    train = train,   # indices de lignes "train"
    test  = test   # indices de lignes "test"
  )
}


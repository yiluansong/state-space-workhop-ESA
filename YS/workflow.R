source("./install_workshop.R")

source("YS/nab.R")
df_daymet
df

df_join <- df_daymet %>%
  left_join(df)

df_binned <- df_join %>%
  gather(key = "var", value = "value", -date, -station) %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(station, year, var) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  spread(key = "var", value = "value") %>%
  ungroup() %>%
  drop_na() %>%
  select(station, year, Acer, Betula, Fraxinus, Morus, Populus, Quercus, Ulmus, everything())


df_pollen_matrix <- df_binned %>% # Select taxa
  select(Acer:Ulmus) %>%
  as.matrix()

df_daymet_matrix_scaled <- df_binned %>% # select covariates
  select(`prcp..mm.day.`, `tmax..deg.c.`) %>%
  as.matrix() %>%
  scale() # Scale covariates

sample_idx <- which(rowSums(df_pollen_matrix) != 0)

# Set-up parameters
p <- ncol(df_daymet_matrix_scaled) + 1 # Number of independent variables plus intercept
n <- ncol(df_pollen_matrix) # number of taxa
V.fixed.glmm <- diag(n)
diag(V.fixed.glmm) <- NA
V.fixed.glmm[1] <- 1
B.fixed.glmm <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), p, n) # reference taxa [,1] are set to 0
B.start.glmm <- matrix(c(rep(0, p), rep(.01, (n - 1) * p)), p, n) # reference taxa [,1] are set to 0

# fit glmm
start_time <- Sys.time()
glmm_mod <- mnGLMM(
  Y = df_pollen_matrix[sample_idx, ],
  X = df_daymet_matrix_scaled[sample_idx, , drop = F],
  B.start = B.start.glmm,
  B.fixed = B.fixed.glmm,
  V.fixed = V.fixed.glmm
)
end_time <- Sys.time()
end_time - start_time
summary(glmm_mod)

# B.start etc
B0.start.mnTS <- glmm_mod$B[1, , drop = F]
B.start.mnTS <- glmm_mod$B[2:p, , drop = F]

sigma.start.mnTS <- glmm_mod$sigma

V.fixed.mnTS <- matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed.mnTS[1] <- 1

V.start.mnTS <- V.fixed.mnTS
V.start.mnTS <- glmm_mod$V

B.fixed.mnTS <- matrix(NA, p - 1, n)
B.fixed.mnTS[, 1] <- 0
B0.fixed.mnTS <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.mnTS <- .5 * diag(n)
C.fixed.mnTS <- C.start.mnTS
C.fixed.mnTS[C.fixed.mnTS != 0] <- NA

# Set-up C with interactions between Fagus and Quercus
C.start.int.mnTS <- matrix(0.001, n, n)
C.start.int.mnTS <- C.start.int.mnTS + (.5 - .001) * diag(n)
C.fixed.int.mnTS <- C.start.int.mnTS
C.fixed.int.mnTS[C.fixed.int.mnTS != 0] <- NA

start_time <- Sys.time()
mnTS_mod <- mnTS(
  Y = df_pollen_matrix[sample_idx, ],
  X = df_daymet_matrix_scaled, Tsample = sample_idx,
  B0.start = B0.start.mnTS, B0.fixed = B0.fixed.mnTS,
  B.start = B.start.mnTS, B.fixed = B.fixed.mnTS,
  C.start = C.start.mnTS, C.fixed = C.fixed.mnTS,
  V.start = V.start.mnTS, V.fixed = V.fixed.mnTS,
  dispersion.fixed = 1, maxit.optim = 1e+2
)
# maxit.optim is the max number of iterations the optimiser will complete before stopping.
# increase maxit.optim if the model needs a lot of time to fit.
end_time <- Sys.time()
end_time - start_time
coef(mnTS_mod)

start_time <- Sys.time()
mnTS_mod_int <- mnTS(
  Y = df_pollen_matrix[sample_idx, ],
  X = df_daymet_matrix_scaled, Tsample = sample_idx,
  B0.start = mnTS_mod$B0, B0.fixed = B0.fixed.mnTS,
  B.start = mnTS_mod$B, B.fixed = B.fixed.mnTS,
  C.start = mnTS_mod$C, C.fixed = C.fixed.int.mnTS,
  V.start = mnTS_mod$V, V.fixed = V.fixed.mnTS,
  dispersion.fixed = 1, maxit.optim = 1e+2
)
end_time <- Sys.time()
end_time - start_time
coef(mnTS_mod_int)

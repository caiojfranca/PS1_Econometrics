library(data.table)
library(foreign)
library(dplyr)
library(plm)

url <- "https://raw.githubusercontent.com/caiojfranca/PS1_Econometrics/main/Problem%20Set%201/cornwell.dta"

data <- read.dta(url)

data <- data[c("county", "year", "lcrmrte", "lprbarr", "lprbconv", "lprbpris", "lavgsen", "lpolpc")]

data <- data %>%
  arrange(county, year)

## --------- ##
## A)
## --------- ##

data_within <- data %>%
  group_by(county) %>% 
  mutate(
    mean_lcrmrte = mean(lcrmrte),
    mean_lprbarr = mean(lprbarr),
    mean_lprbconv= mean(lprbconv),
    mean_lprbpris= mean(lprbpris),
    mean_lavgsen = mean(lavgsen),
    mean_lpolpc  = mean(lpolpc),
    
    lcrmrte_within = lcrmrte - mean_lcrmrte,
    lprbarr_within  = lprbarr  - mean_lprbarr,
    lprbconv_within = lprbconv - mean_lprbconv,
    lprbpris_within = lprbpris - mean_lprbpris,
    lavgsen_within  = lavgsen  - mean_lavgsen,
    lpolpc_within   = lpolpc   - mean_lpolpc
  ) %>%
  ungroup()

Y_tilde <- as.matrix(data_within$lcrmrte_within)

X_tilde <- data_within %>%
  select(lprbarr_within, lprbconv_within, lprbpris_within, lavgsen_within, lpolpc_within) %>%
  as.matrix()                                     

beta_FE <- solve(t(X_tilde) %*% X_tilde) %*% (t(X_tilde) %*% Y_tilde)
beta_FE

## --------- ##
## B)
## --------- ##

N    <- n_distinct(data$county)
Time <- n_distinct(data$year)     
K    <- ncol(X_tilde)              

u_hat <- as.vector(Y_tilde - X_tilde %*% beta_FE)

sigma2_e_hat <- sum(u_hat^2) / (N * (Time - 1) - K)
sigma2_e_hat

data_between <- data %>%
  group_by(county) %>%
  summarise(
    Y_bar        = mean(lcrmrte),
    lprbarr_bar  = mean(lprbarr),
    lprbconv_bar = mean(lprbconv),
    lprbpris_bar = mean(lprbpris),
    lavgsen_bar  = mean(lavgsen),
    lpolpc_bar   = mean(lpolpc),
    .groups = "drop"
  )

Y_bar <- as.matrix(data_between$Y_bar)

X_bar <- data_between %>%
  select(lprbarr_bar, lprbconv_bar, lprbpris_bar, lavgsen_bar, lpolpc_bar) %>%
  mutate(const = 1) %>%
  as.matrix()                       

beta_BE <- solve(t(X_bar) %*% X_bar) %*% (t(X_bar) %*% Y_bar)

e_be <- as.vector(Y_bar - X_bar %*% beta_BE)

K_BE <- ncol(X_bar)

sigma2_alpha_hat <- (sum(e_be^2) / (N - K_BE)) - (1 / Time) * sigma2_e_hat
sigma2_alpha_hat

lambda_hat <- 1 - (sqrt(sigma2_e_hat) / sqrt(Time * sigma2_alpha_hat + sigma2_e_hat))
lambda_hat

data_RE <- data_within %>%
  mutate(
    const = 1,
    
    lcrmrte_RE = lcrmrte - lambda_hat * mean_lcrmrte,
    
    lprbarr_RE  = lprbarr  - lambda_hat * mean_lprbarr,
    lprbconv_RE = lprbconv - lambda_hat * mean_lprbconv,
    lprbpris_RE = lprbpris - lambda_hat * mean_lprbpris,
    lavgsen_RE  = lavgsen  - lambda_hat * mean_lavgsen,
    lpolpc_RE   = lpolpc   - lambda_hat * mean_lpolpc
  )

Y_RE <- as.matrix(data_RE$lcrmrte_RE)

X_RE <- data_RE %>%
  select(const, lprbarr_RE, lprbconv_RE, lprbpris_RE, lavgsen_RE, lpolpc_RE) %>%
  as.matrix()

K_RE <- ncol(X_RE)   

beta_RE <- solve(t(X_RE) %*% X_RE) %*% (t(X_RE) %*% Y_RE)

beta_RE_short <- beta_RE[-1, , drop = FALSE]
beta_RE_short

## --------- ##
## C)
## --------- ##

Var_FE <- sigma2_e_hat * solve(t(X_tilde) %*% X_tilde)

v_hat <- as.vector(Y_RE - X_RE %*% beta_RE)

sigma2_RE_hat <- sum(v_hat^2) / ((N * Time) - K_RE) 

Var_RE <- sigma2_RE_hat * solve(t(X_RE) %*% X_RE)

Var_RE_short <- Var_RE[-1, -1, drop = FALSE]

beta_diff <- beta_FE - beta_RE_short

Var_diff <- Var_FE - Var_RE_short

H_stat <- abs(as.numeric(t(beta_diff) %*% solve(Var_diff) %*% beta_diff))
p_value <- 1 - pchisq(H_stat, df = nrow(beta_diff))

H_stat
p_value

## --------- ##
## D)
## --------- ##

data_diff <- data %>%
  arrange(county, year) %>%            
  group_by(county) %>%
  mutate(
    d_crmrte   = lcrmrte   - dplyr::lag(lcrmrte),
    d_lprbarr  = lprbarr  - dplyr::lag(lprbarr),
    d_lprbconv = lprbconv - dplyr::lag(lprbconv),
    d_lprbpris = lprbpris - dplyr::lag(lprbpris),
    d_lavgsen  = lavgsen  - dplyr::lag(lavgsen),
    d_lpolpc   = lpolpc   - dplyr::lag(lpolpc)
  ) %>%
  filter(dplyr::row_number() > 1) %>%
  ungroup()

Y_D <- as.matrix(data_diff$d_crmrte)

X_D <- data_diff %>%
  select(d_lprbarr, d_lprbconv, d_lprbpris, d_lavgsen, d_lpolpc) %>%
  mutate(const = 1) %>%
  as.matrix()

beta_D <- solve(t(X_D) %*% X_D) %*% (t(X_D) %*% Y_D)
beta_D_short <- beta_D[-nrow(beta_D), , drop = FALSE]
beta_D_short


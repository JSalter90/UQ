# h/t Daniel Williamson, UNRISK training, April 26

# ---- EBM ----
simulate_ebm <- function(years, F, lambda = 1.2, C = 10, T0 = 0, dt = 1) {
  n <- length(years)
  T <- numeric(n)
  T[1] <- T0
  for (i in 2:n) T[i] <- T[i-1] + (dt / C) * (F[i-1] - lambda * T[i-1])
  T
}

# ---- CO2 forcing ----
forcing_co2_of_C <- function(C, C0 = 280) 5.35 * log(C / C0)

co2_path_knots <- function(years,
                           knots_year = c(1850, 1980, 2020, 2100),
                           knots_C    = c(280, 338, 415, 560)) {
  approx(x = knots_year, y = knots_C, xout = years, rule = 2)$y
}

# ---- Volcano forcing: fixed events (known history) ----
forcing_volcano_from_events <- function(years, events, tau_v = 1.2) {
  n <- length(years)
  Fv <- rep(0, n)
  if (is.null(events) || nrow(events) == 0) return(Fv)
  for (j in seq_len(nrow(events))) {
    t0 <- match(events$year[j], years)
    if (is.na(t0)) next
    A <- events$amp[j]
    k <- 0:(n - t0)
    Fv[t0:n] <- Fv[t0:n] + A * exp(-k / tau_v)
  }
  Fv
}

# ---- Volcano forcing: random future ----
forcing_volcano_random_future <- function(years, start_year = 2020,
                                          rate = 0.08, amp_mean = -2.0, amp_sd = 0.6,
                                          tau_v = 1.2) {
  n <- length(years)
  Fv <- rep(0, n)
  idx_start <- which(years >= start_year)
  if (length(idx_start) == 0) return(Fv)
  i0 <- idx_start[1]
  events <- i0:n
  hit <- events[runif(length(events)) < rate]
  if (length(hit) == 0) return(Fv)
  
  for (t0 in hit) {
    A <- rnorm(1, mean = amp_mean, sd = amp_sd)
    k <- 0:(n - t0)
    Fv[t0:n] <- Fv[t0:n] + A * exp(-k / tau_v)
  }
  Fv
}

# Developing a history + future
years_all  <- 1850:2100 # past
years_proj <- 2020:2100  # future

C0_ref <- 280
C_t <- co2_path_knots(years_all)
F_co2_all <- forcing_co2_of_C(C_t, C0 = C0_ref)

events_fixed <- data.frame(
  year = c(1883, 1963, 1982, 1991),
  amp  = c(-1.0, -1.0, -1.2, -2.5)
)
F_volc_fixed_all <- forcing_volcano_from_events(years_all, events_fixed, tau_v = 1.2)

F_base_all <- F_co2_all + F_volc_fixed_all

set.seed(1)

F_volc_future <- forcing_volcano_random_future(years_all, start_year = 2020, rate = 0.08)
F_one <- F_base_all + F_volc_future

T_one_all <- simulate_ebm(years_all, F_one, lambda = 1.2, C = 10, T0 = 0)

#idx_proj <- years_all %in% years_proj
idx_proj <- years_all %in% years_all
df_one <- data.frame(
  year = years_all[idx_proj],
  forcing = F_one[idx_proj],
  temp = T_one_all[idx_proj]
)

df_one_long <- rbind(
  data.frame(year = df_one$year, panel = "Forcing  F(t)  (W/m^2)", value = df_one$forcing),
  data.frame(year = df_one$year, panel = "Temperature  T(t)  (K since 1850)", value = df_one$temp)
)

# p_two_panel <- ggplot(df_one_long, aes(year, value)) +
#   geom_hline(yintercept = 0, linetype = 3, linewidth = 0.4) +
#   geom_line(linewidth = 0.9) +
#   facet_wrap(~panel, ncol = 1, scales = "free_y") +
#   labs(x = "Year", y = NULL) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor = element_blank(),
#         strip.text = element_text(face = "bold"))

#p_two_panel

set.seed(2)

idx_proj <- years_all %in% years_proj

M <- 2000
k_show <- 60

r_lognorm <- function(n, meanlog, sdlog) exp(rnorm(n, meanlog, sdlog))
lambda_draw <- r_lognorm(M, meanlog = log(1.2), sdlog = 0.25)
C_draw      <- r_lognorm(M, meanlog = log(10),  sdlog = 0.35)

T_proj_mat <- matrix(NA_real_, nrow = length(years_proj), ncol = M)

for (m in 1:M) {
  F_future_m <- forcing_volcano_random_future(years_all, start_year = 2020, rate = 0.08)
  F_m <- F_base_all + F_future_m
  T_m_all <- simulate_ebm(years_all, F_m, lambda = lambda_draw[m], C = C_draw[m], T0 = 0)
  T_proj_mat[, m] <- T_m_all[idx_proj]
}

qs <- apply(T_proj_mat, 1, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
df_q <- data.frame(year = years_proj, lo = qs[1,], mid = qs[2,], hi = qs[3,])

set.seed(3)
idx <- sample(seq_len(M), min(k_show, M))
df_spag <- data.frame(
  year = rep(years_proj, times = length(idx)),
  draw = rep(seq_along(idx), each = length(years_proj)),
  T = as.vector(T_proj_mat[, idx])
)

# p_traj <- ggplot() +
#   geom_line(data = df_spag, aes(year, T, group = draw), alpha = 0.25) +
#   geom_ribbon(data = df_q, aes(year, ymin = lo, ymax = hi), alpha = 0.20) +
#   geom_line(data = df_q, aes(year, mid), linewidth = 1.05) +
#   labs(x = "Year", y = "Temperature  T(t)  (K since 1850)") +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor = element_blank())
# 
# p_traj

thr <- 2
T_end <- T_proj_mat[nrow(T_proj_mat), ]

p_hat <- mean(T_end > thr)
se_hat <- sqrt(p_hat * (1 - p_hat) / length(T_end))
lab <- sprintf("P(T(2100) > %.0f) ≈ %.3f\nMC SE ≈ %.3f", thr, p_hat, se_hat)

h <- hist(T_end, breaks = 40, plot = FALSE)
df_hist <- data.frame(
  left = h$breaks[-length(h$breaks)],
  right = h$breaks[-1],
  count = h$counts
)
df_hist$mid <- 0.5 * (df_hist$left + df_hist$right)
df_hist$exceed <- df_hist$mid > thr

x_min <- min(df_hist$left); x_max <- max(df_hist$right)
x_lab <- x_min + 0.05 * (x_max - x_min)
y_top <- max(df_hist$count)

# p_hist <- ggplot(df_hist) +
#   geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = count, alpha = exceed),
#             colour = "white") +
#   scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.85), guide = "none") +
#   geom_vline(xintercept = thr, linetype = 2, linewidth = 0.9) +
#   annotate("label", x = 3.4, y = 0.95 * y_top, hjust = 0, vjust = 1,
#            label = lab, size = 4) +
#   labs(x = "T(2100) (K since 1850)", y = "Count") +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor = element_blank())
# 
# p_hist

p_exc <- rowMeans(T_proj_mat > thr, na.rm = TRUE)
df_exc <- data.frame(year = years_proj, p = p_exc)

# p_exceed <- ggplot(df_exc, aes(year, p)) +
#   geom_line(linewidth = 1.05) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "Year", y = sprintf("P(T(t) > %.0f)", thr)) +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor = element_blank())
# 
# p_exceed

F2x <- 5.35 * log(2)
ECS <- F2x / lambda_draw

df_ecs <- data.frame(ECS = ECS)
q <- quantile(ECS, c(0.05, 0.5, 0.95))
lab_ecs <- sprintf("Prior ECS (K)\n5%%=%.2f  50%%=%.2f  95%%=%.2f", q[1], q[2], q[3])

h2 <- hist(ECS, breaks = 40, plot = FALSE)
dfh2 <- data.frame(left = h2$breaks[-length(h2$breaks)],
                   right = h2$breaks[-1],
                   count = h2$counts)
y2_top <- max(dfh2$count)
x2_min <- min(dfh2$left); x2_max <- max(dfh2$right)
x2_lab <- x2_min + 0.05*(x2_max - x2_min)

# p_ecs <- ggplot(dfh2) +
#   geom_rect(aes(xmin = left, xmax = right, ymin = 0, ymax = count),
#             colour = "white", alpha = 0.7) +
#   geom_vline(xintercept = q, linetype = c(2,1,2), linewidth = 0.9) +
#   annotate("label", x = 4.5, y = 0.95*y2_top, hjust = 0, vjust = 1,
#            label = lab_ecs, size = 4) +
#   labs(x = "ECS (K)", y = "Count") +
#   theme_minimal(base_size = 14) +
#   theme(panel.grid.minor = element_blank())
# 
# p_ecs





# Use same ranges for lambda, C as previously
mu_lambda <- log(1.2); sd_lambda <- 0.25
mu_C      <- log(10);  sd_C      <- 0.35
lambda_rng <- qlnorm(c(0.001, 0.999), mu_lambda, sd_lambda)
C_rng      <- qlnorm(c(0.001, 0.999), mu_C,      sd_C)

# map unit-cube -> physical parameters, linear in log-space (stable)
unit_to_theta <- function(Xu) {
  Xu <- as.matrix(Xu)
  lam <- exp(log(lambda_rng[1]) + Xu[,1] * (log(lambda_rng[2]) - log(lambda_rng[1])))
  CC  <- exp(log(C_rng[1])      + Xu[,2] * (log(C_rng[2])      - log(C_rng[1])))
  list(lambda = lam, C = CC)
}
# Reverse
theta_to_unit <- function(lambda, C) {
  u1 <- (log(lambda) - log(lambda_rng[1])) / (log(lambda_rng[2]) - log(lambda_rng[1]))
  u2 <- (log(C)      - log(C_rng[1]))      / (log(C_rng[2])      - log(C_rng[1]))
  cbind(u1, u2)
}

# Run simulator
ebm_output <- function(X, forcing = F_one){
  th <- unit_to_theta(X)
  out <- matrix(NA_real_, length(years_all), nrow(X))
  for (i in 1:nrow(X)){
    out[,i] <- simulate_ebm(years_all, forcing, lambda = th$lambda[i], C = th$C[i], T0 = 0)
  }
  return(out)
}

# Theoph 群体 PK, nlmixr2 FOCE-I
# 先跑 run_analysis.py (要用 results/model_fits.csv)
# M1 -> M5, 最后用 M5

suppressPackageStartupMessages({
  library(nlmixr2)
  library(ggplot2)
  library(patchwork)
})

set.seed(1110)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

BLUE <- "steelblue"; ORANGE <- "darkorange"; DARK <- "grey20"; GREY <- "grey85"
theme_set(theme_bw(base_size = 9))
save_plot <- function(p, file, w, h) ggsave(file.path("figures", file), p, width = w, height = h, dpi = 150, bg = "white")

cv_pct <- function(omega2) 100 * sqrt(exp(omega2) - 1)  # omega^2 (log scale) -> CV%


# data ----
raw <- read.csv("data/theoph.csv")
raw$dose_mg <- raw$Dose * raw$Wt  # dose in the file is mg/kg

subjects <- unique(raw[, c("Subject", "Wt", "dose_mg")])

# 每个人一条给药记录 (EVID = 1, 进 depot)
dose_rows <- data.frame(ID = subjects$Subject, TIME = 0, DV = NA, AMT = subjects$dose_mg,
                        EVID = 1, CMT = "depot", WT = subjects$Wt)

# 观测记录 (EVID = 0, CMT = center)。TIME = 0 的点去掉，和 python 一样：
# 有 3 个人给药前就有一点浓度，模型从 0 开始拟合不了
obs <- subset(raw, Time > 0)
obs_rows <- data.frame(ID = obs$Subject, TIME = obs$Time, DV = obs$conc, AMT = 0,
                       EVID = 0, CMT = "center", WT = obs$Wt)

dat <- rbind(dose_rows, obs_rows)
dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]


# models ----
# thetas on log scale; combined error sd = sqrt(add^2 + (prop*pred)^2)

m1_base <- function() {
  ini({
    tka <- log(1.5); label("Absorption rate constant ka (1/h)")
    tcl <- log(2.8); label("Apparent clearance CL/F (L/h)")
    tv  <- log(32);  label("Apparent volume V/F (L)")
    eta.ka ~ 0.5
    eta.cl ~ 0.1
    eta.v  ~ 0.05
    add.sd  <- 0.2; label("Additive residual error (mg/L)")
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# M2: 加吸收滞后 tlag（python 的残差前 2 h 有问题）。tlag 先不加 eta，
# 每个人 0.3 h 之前只有一个点
m2_lag <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("Apparent clearance CL/F (L/h)")
    tv    <- log(32);   label("Apparent volume V/F (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    eta.ka ~ 0.5
    eta.cl ~ 0.1
    eta.v  ~ 0.05
    add.sd  <- 0.2; label("Additive residual error (mg/L)")
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl)
    v    <- exp(tv  + eta.v)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# M3: 加体重，CL ~ WT^0.75, V ~ WT (70 kg)，指数固定
m3_lag_wt <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("CL/F for a 70 kg adult (L/h)")
    tv    <- log(32);   label("V/F for a 70 kg adult (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    eta.ka ~ 0.5
    eta.cl ~ 0.1
    eta.v  ~ 0.05
    add.sd  <- 0.2; label("Additive residual error (mg/L)")
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl) * (WT / 70)^0.75
    v    <- exp(tv  + eta.v)  * (WT / 70)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# 默认的 covMethod "r,s" 对有 tlag 的模型会退成 "s"，CL/F 的 RSE 变成 40-50%，所以用 "r"
fit_model <- function(model) {
  nlmixr2(model, dat, est = "focei",
          control = foceiControl(print = 0, covMethod = "r"),
          table = tableControl(cwres = TRUE))
}

# M4: 加 CL-V 的相关
# M5: M3 里 add.sd 的 RSE ~70%，M4 里变成 ~0.01，试试只用比例误差
# m4_corr <- m3_lag_wt |> ini(eta.cl + eta.v ~ c(0.1, 0.01, 0.05))
m4_corr <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("CL/F for a 70 kg adult (L/h)")
    tv    <- log(32);   label("V/F for a 70 kg adult (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    eta.ka ~ 0.5
    eta.cl + eta.v ~ c(0.1, 0.01, 0.05)
    add.sd  <- 0.2; label("Additive residual error (mg/L)")
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl) * (WT / 70)^0.75
    v    <- exp(tv  + eta.v)  * (WT / 70)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

m5_prop <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("CL/F for a 70 kg adult (L/h)")
    tv    <- log(32);   label("V/F for a 70 kg adult (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    eta.ka ~ 0.5
    eta.cl + eta.v ~ c(0.1, 0.01, 0.05)
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl) * (WT / 70)^0.75
    v    <- exp(tv  + eta.v)  * (WT / 70)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

# M5 去掉体重试一下。体重 dOFV 只有 -1.6，看看去掉以后参数变多少
m5_nowt <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("Apparent clearance CL/F (L/h)")
    tv    <- log(32);   label("Apparent volume V/F (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    eta.ka ~ 0.5
    eta.cl + eta.v ~ c(0.1, 0.01, 0.05)
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl)
    v    <- exp(tv  + eta.v)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

# CL 的体重指数不固定 0.75，让模型自己估 (V 还是固定 1)
m5_estwt <- function() {
  ini({
    tka   <- log(1.5);  label("Absorption rate constant ka (1/h)")
    tcl   <- log(2.8);  label("CL/F for a 70 kg adult (L/h)")
    tv    <- log(32);   label("V/F for a 70 kg adult (L)")
    ttlag <- log(0.15); label("Absorption lag time tlag (h)")
    dwt.cl <- 0.75;     label("WT exponent on CL")
    eta.ka ~ 0.5
    eta.cl + eta.v ~ c(0.1, 0.01, 0.05)
    prop.sd <- 0.1; label("Proportional residual error")
  })
  model({
    ka   <- exp(tka + eta.ka)
    cl   <- exp(tcl + eta.cl) * (WT / 70)^dwt.cl
    v    <- exp(tv  + eta.v)  * (WT / 70)
    tlag <- exp(ttlag)
    d/dt(depot)  <- -ka * depot
    alag(depot)  <- tlag
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

fits <- list(M1 = fit_model(m1_base),
             M2 = fit_model(m2_lag),
             M3 = fit_model(m3_lag_wt),
             M4 = fit_model(m4_corr),
             M5 = fit_model(m5_prop))
             # M5_noWT = fit_model(m5_nowt),    # 去掉体重 OFV +2.5，还是留着体重
             # M5_estWT = fit_model(m5_estwt))  # 指数估不出来 (0.58, RSE 126%)，固定 0.75

# OFV = -2LL，多 1 个参数要降 3.84 以上才算显著 (p < 0.05)。M5 是少了 1 个参数，OFV 升一点没关系
comparison <- data.frame(
  model = c("M1 1-cmt", "M2 + tlag", "M3 + WT", "M4 + corr(CL,V)", "M5 = M4, prop error only"),
  # , "M5 without WT", "M5, CL exponent estimated"),
  n_par = sapply(fits, function(f) attr(logLik(f), "df")),
  OFV   = sapply(fits, function(f) f$objf),
  AIC   = sapply(fits, AIC)
)
comparison$dOFV_vs_previous <- c(NA, diff(comparison$OFV))
write.csv(comparison, "results/nlme_model_comparison.csv", row.names = FALSE)

cat("\n=== Model comparison ===\n")
print(comparison, row.names = FALSE, digits = 4)
# print(fits$M5_estWT$parFixedDf["dwt.cl", ])

final <- fits$M5
final_name <- "M5"


# diagnostics ----
params <- final$parFixedDf
# nlmixr2 的 %RSE 是在估计的尺度上算的，log 参数 (tka, tcl, tv, ttlag) 这样不对
# (V/F 只有 ~1%)。log 尺度上的 SE 本身约等于原尺度的 RSE
log_par <- c("tka", "tcl", "tv", "ttlag")
params$RSE_natural_pct <- ifelse(rownames(params) %in% log_par, 100 * params$SE, params[["%RSE"]])
write.csv(cbind(parameter = rownames(params), params), "results/nlme_parameters.csv", row.names = FALSE)
cat("\n=== Final model (", final_name, ") parameters ===\n", sep = "")
print(params)

# shrinkage
etas <- final$eta
omega <- diag(final$omega)
shrinkage <- data.frame(
  eta       = names(omega),
  omega2    = omega,
  CV_pct    = cv_pct(omega),
  shrink_pct = 100 * (1 - sapply(names(omega), function(e) sd(etas[[e]])) / sqrt(omega))
)
write.csv(shrinkage, "results/nlme_shrinkage.csv", row.names = FALSE)
cat("\n=== Between-subject variability and shrinkage ===\n")
print(shrinkage, row.names = FALSE, digits = 3)

# Goodness-of-fit plots
d <- as.data.frame(final)
d <- d[!is.na(d$DV), ]
lim <- range(c(d$DV, d$PRED, d$IPRED))
gof_xy <- function(x, xlab) {
  ggplot(d, aes(.data[[x]], DV)) +
    geom_abline(colour = DARK, linetype = "dashed") +
    geom_point(colour = BLUE, alpha = 0.7, size = 1.3) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ORANGE, linewidth = 0.7) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = xlab, y = "Observed (mg/L)")
}
gof_res <- function(x, xlab) {
  ggplot(d, aes(.data[[x]], CWRES)) +
    geom_hline(yintercept = 0, colour = DARK, linetype = "dashed") +
    geom_hline(yintercept = c(-2, 2), colour = GREY) +
    geom_point(colour = BLUE, alpha = 0.7, size = 1.3) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ORANGE, linewidth = 0.7) +
    labs(x = xlab, y = "CWRES")
}
p_gof <- (gof_xy("PRED", "Population prediction (mg/L)") | gof_xy("IPRED", "Individual prediction (mg/L)")) /
  (gof_res("TIME", "Time (h)") | gof_res("PRED", "Population prediction (mg/L)")) +
  plot_annotation(title = paste0("GOF (", final_name, ")"))
save_plot(p_gof, "nlme_gof.png", 8, 7.5)

# 个体拟合图：PRED (eta = 0) 和 IPRED
ap <- augPred(final)
ap_lines <- subset(ap, ind != "Observed")
ap_obs <- subset(ap, ind == "Observed")
p_ind <- ggplot() +
  geom_line(data = ap_lines, aes(time, values, colour = ind, linetype = ind), linewidth = 0.6) +
  geom_point(data = ap_obs, aes(time, values), colour = DARK, size = 1.2) +
  facet_wrap(~ id, ncol = 4, labeller = labeller(id = function(x) paste("Subject", x))) +
  scale_colour_manual(values = c(Population = ORANGE, Individual = BLUE)) +
  scale_linetype_manual(values = c(Population = "dashed", Individual = "solid")) +
  labs(x = "Time (h)", y = "Concentration (mg/L)", colour = NULL, linetype = NULL,
       title = paste0("individual fits (", final_name, ")")) +
  theme(legend.position = "top")
save_plot(p_ind, "nlme_individual_fits.png", 9, 7)

# eta vs 体重，M2 和 M3 比一下
eta_wt <- function(fit, label) {
  e <- merge(fit$eta, subjects, by.x = "ID", by.y = "Subject")
  rbind(data.frame(model = label, eta = "eta.cl", WT = e$Wt, value = e$eta.cl),
        data.frame(model = label, eta = "eta.v",  WT = e$Wt, value = e$eta.v))
}
ew <- rbind(eta_wt(fits$M2, "M2 (no weight)"), eta_wt(fits$M3, "M3 (allometric weight)"))
p_eta <- ggplot(ew, aes(WT, value)) +
  geom_hline(yintercept = 0, colour = DARK, linetype = "dashed") +
  geom_point(colour = BLUE, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = ORANGE, linewidth = 0.7) +
  facet_grid(eta ~ model) +
  labs(x = "Body weight (kg)", y = "Individual eta", title = "eta vs WT")
save_plot(p_eta, "nlme_eta_vs_weight.png", 7, 5)

# VPC, 500 replicates
p_vpc <- vpcPlot(final, n = 500, show = list(obs_dv = TRUE), log_y = FALSE) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       title = paste0("VPC (", final_name, ")"))
save_plot(p_vpc, "nlme_vpc.png", 7, 4.5)


# dosing simulation ----
# 模拟 1000 个病人：体重从 12 个人里重抽，eta 按 omega 抽。没考虑参数的不确定性 (SE)
nSim <- 1000
th <- final$theta
eta_sim <- MASS::mvrnorm(nSim, mu = rep(0, nrow(final$omega)), Sigma = final$omega)
wt_sim <- sample(subjects$Wt, nSim, replace = TRUE)
patients <- data.frame(
  ka   = exp(th[["tka"]] + eta_sim[, "eta.ka"]),
  cl   = exp(th[["tcl"]] + eta_sim[, "eta.cl"]) * (wt_sim / 70)^0.75,
  v    = exp(th[["tv"]]  + eta_sim[, "eta.v"])  * (wt_sim / 70),
  tlag = exp(th[["ttlag"]])
)

sim_model <- rxode2({
  d/dt(depot)  <- -ka * depot
  alag(depot)  <- tlag
  d/dt(center) <-  ka * depot - cl / v * center
  cp <- center / v
})

regimens <- data.frame(label = c("300 mg every 8 h", "450 mg every 12 h"),
                       dose = c(300, 450), tau = c(8, 12))
t_end <- 96
sim_all <- do.call(rbind, lapply(seq_len(nrow(regimens)), function(i) {
  r <- regimens[i, ]
  ev <- et(amt = r$dose, ii = r$tau, addl = t_end / r$tau - 1, cmt = "depot") |>
    et(seq(0, t_end, by = 0.1))
  s <- as.data.frame(rxSolve(sim_model, params = patients, events = ev))
  s$regimen <- r$label
  s$tau <- r$tau
  s[, c("sim.id", "time", "cp", "regimen", "tau")]
}))

# 稳态 = 最后一个给药间隔
ss <- subset(sim_all, time >= t_end - tau)
per_patient <- aggregate(cp ~ sim.id + regimen, data = ss,
                         FUN = function(x) c(min = min(x), max = max(x)))
per_patient <- data.frame(regimen = per_patient$regimen,
                          cmin = per_patient$cp[, "min"], cmax = per_patient$cp[, "max"])

target <- function(pp, lo, hi) {
  do.call(rbind, lapply(split(pp, pp$regimen), function(g) data.frame(
    regimen = g$regimen[1], range = paste0(lo, "-", hi, " mg/L"),
    median_Cmin = median(g$cmin), median_Cmax = median(g$cmax),
    pct_whole_interval_in_range = 100 * mean(g$cmin >= lo & g$cmax <= hi),
    pct_trough_below = 100 * mean(g$cmin < lo),
    pct_peak_above   = 100 * mean(g$cmax > hi)
  )))
}
# regimen_table <- rbind(target(per_patient, 10, 20), target(per_patient, 5, 15))
regimen_table <- target(per_patient, 10, 20)
write.csv(regimen_table, "results/nlme_regimens.csv", row.names = FALSE)
cat("\n=== Steady-state target attainment, 1000 virtual patients ===\n")
print(regimen_table, row.names = FALSE, digits = 3)

band <- do.call(rbind, lapply(split(sim_all, list(sim_all$regimen, sim_all$time), drop = TRUE), function(g)
  data.frame(regimen = g$regimen[1], time = g$time[1],
             p05 = quantile(g$cp, 0.05), p50 = median(g$cp), p95 = quantile(g$cp, 0.95))))
p_reg <- ggplot(band, aes(time)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 20, fill = GREY, alpha = 0.6) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = BLUE, alpha = 0.25) +
  geom_line(aes(y = p50), colour = BLUE, linewidth = 0.7) +
  facet_wrap(~ regimen, ncol = 1) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       title = "1000 simulated patients, 900 mg/day")
save_plot(p_reg, "nlme_regimens.png", 8, 5.5)

cat("\nDone. Tables in results/nlme_*.csv, figures in figures/nlme_*.png\n")

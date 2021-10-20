#!/usr/bin/env Rscript

#' Analyze calibration and select best-fitting runs
#'
#' @description
#' imports model runs and finds the runs that fit empirical data best, for use
#' in future analyses
#'
#' @param TODO need to make function
#'
#' @returns
#' writes parameters for calibrated runs to file
#'

rm(list = ls())
# Load required packages
library(dplyr)
library(openxlsx)
library(abind)
library(ggplot2)
library(xlsx)

# Load scripts
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("parallel.R")
source("data_input.R")


cal_results <- NULL
# read in and bind calibration results
for (batch.ind in 1:10) {
  temp_results <- readRDS(
    paste0("calibration/CalibrationOutputs", batch.ind, ".rds")
  )
  cal_results <- rbind(cal_results, temp_results)
}
rm(temp_results)

# add the calibration parameters to the calibration results
calibration_params <- readRDS(paste0("calibration/Calib_par_table.rds"))
cal_results <- cbind(cal_results, calibration_params)

# read in workbook of targets
workbook <- openxlsx::loadWorkbook("Inputs/MasterTable.xlsx")
target <- openxlsx::read.xlsx(workbook, sheet = "Target")
target_data <- target$pe


# Calculate goodness of fit (gof)
for (ss in seq_len(nrow(cal_results))) {
  prediction <- cal_results[ss, c(
    "od_death16", "od_death17", "od_death18", "od_death19",
    "fx.death16", "fx.death17", "fx.death18", "fx.death19",
    "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19"
  )]
  gof <- 0

  # TODO needs to be less hardcoded at some point
  for (j in seq_len(length(target_data))) {
    if (j %in% c(5:8)) { # 5:8 is proportion that involves fx
      tmp <- abs(prediction[j] * 100 - target_data[j] * 100)
      adj <-  tmp / (target_data[j] * 100)
    } else {
      adj <- (abs(prediction[j] - target_data[j]) / target_data[j])
    }
    gof <- gof + adj / length(target_data)
  }
  cal_results[ss, "gof"] <- gof
}

# sort by goodness of fit and select top 100 fits
sorted_results <- cal_results[order(cal_results[, "gof"], decreasing = F), ]
cal_sample <- 100
cal_results <- sorted_results[1:cal_sample, ]
write.xlsx(cal_results,
  file = paste0("calibration/Calibrated_results.xlsx"),
  col.names = T, row.names = F
)


## save calibrated results as parameter lists (prepare for main analysis)##
# INPUT PARAMETERS
cal_param_names <- names(calibration_params)
cal_params <- list()


# for each selected calibration run, determine the run's parameters
for (cc in seq_len(nrow(cal_results))) {
  for (pp in seq_len(length(cal_param_names))) {
    cal_params[[cal_param_names[pp]]] <- cal_results[cc, cal_param_names[pp]]
  }
  # Overdose probability parameter matrix (per month)
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["rx", "subs"] <- params$od_rx_sub
  overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
  overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
  overdose_probs["NODU", "subs"] <- params$od_nodu_sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
  params$overdose_probs <- overdose_probs

  # Baseline mortality parameters, excluding overdose (per month)
  mortality_probs <- matrix(0, nrow = 2, ncol = length(mor_gp))
  rownames(mortality_probs) <- c("bg", "drug")
  colnames(mortality_probs) <- mor_gp
  mortality_probs["bg", ] <- params$mortality_base
  mortality_probs["drug", ] <- params$mortality_drug
  params$mortality_probs <- mortality_probs
  params$od_911_pub <- params$od_911_priv * params$od_911_pub_mul

  cal_params[[cc]] <- params
}

saveRDS(cal_params, file = paste0("calibration/CalibratedData.rds"))
saveRDS(cal_results[, "seed"], file = paste0("calibration/CalibratedSeed.rds"))


## Plot ##
par(mfrow = c(1, 3))
# REVIEWED md = median, change to mean
mean_oddeath <- apply(
  cal_results[, c("od_death16", "od_death17", "od_death18", "od_death19")],
  2,
  mean
)
ymax <- max(
  cal_results[, c("od_death16", "od_death17", "od_death18", "od_death19")]
)
plot(
  x = 2016:2019,
  target_data[1:4],
  col = "black",
  pch = 18,
  xlab = "Year",
  ylab = "Number of opioid overdose deaths",
  cex = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3,
  xaxt = "n",
  ylim = c(0, 500),
  frame.plot = FALSE
)
for (i in 1:cal_sample) {
  res_name <- c("od_death16", "od_death17", "od_death18", "od_death19")
  lines(
    x = 2016:2019,
    y = cal_results[i, res_name],
    col = adjustcolor("indianred1", alpha = 0.2),
    lwd = 3
  )
}
lines(x = 2016:2019, y = mean_oddeath, col = "red", lwd = 4)
points(x = 2016:2019,
       target_data[1:4],
       col = "black",
       pch = 16,
       cex = 1.2,
       cex.axis = 0.95
      )
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)

abline(h = 0)
legend("top",
  legend = c("Target", "Model: mean", "Model: simulation"),
  col = c("black", "red", adjustcolor("indianred1", alpha = 0.2)),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = 3,
  bty = "n",
  pt.cex = 1.2,
  cex = 1.1,
  text.col = "black",
  horiz = T
)

death_fx_names <- c("fx.death16", "fx.death17", "fx.death18", "fx.death19")
mean_fxdeath <- apply(cal_results[, death_fx_names], 2, median)
# REVIEWED: change from hardcoded ylim; mean instead of median for everything
plot(
  x = 2016:2019,
  target_data[5:8],
  col = "black",
  pch = 18,
  xlab = "Year",
  ylab = "Proportion of OOD deaths with fentanyl present",
  cex = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3,
  xaxt = "n",
  yaxt = "n",
  ylim = c(0, 1),
  frame.plot = FALSE
)

for (i in 1:cal_sample) {
  lines(x = 2016:2019,
        y = cal_results[i,
                    c("fx.death16", "fx.death17", "fx.death18", "fx.death19")],
        col = adjustcolor("indianred1", alpha = 0.2), lwd = 2
      )
}

lines(x = 2016:2019, y = mean_fxdeath, col = "red", lwd = 3)
points(x = 2016:2019,
       target_data[5:8],
       col = "black",
       pch = 16,
       cex = 1.2,
       cex.axis = 0.95
      )

axis(1,
     at = 2016:2019,
     pos = 0,
     lwd.ticks = 0,
     cex.axis = 1.2
    )

axis(2,
     at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
     labels = c("0", "20%", "40%", "60%", "80%", "100%"),
     cex.axis = 1.2
    )

abline(h = 0)


md_edvisits <- apply(
  cal_results[, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")],
  2,
  median
)
plot(
  x = 2016:2019,
  target_data[9:12],
  col = "black",
  pch = 18,
  xlab = "Year",
  ylab = "Number of ED visists for opioid overdose",
  cex = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3,
  xaxt = "n",
  ylim = c(0, 2500),
  frame.plot = FALSE
)

for (i in 1:cal_sample) {
  lines(
    x = 2016:2019,
    y = cal_results[i, c(
        "ed.visit16",
        "ed.visit17",
        "ed.visit18",
        "ed.visit19"
      )],
    col = adjustcolor("indianred1", alpha = 0.2),
    lwd = 2
  )
}
lines(x = 2016:2019, y = md_edvisits, col = "red", lwd = 3)
points(
  x = 2016:2019,
  target_data[9:12],
  col = "black",
  pch = 16,
  cex = 1.2,
  cex.axis = 1.2
)
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)
abline(h = 0)

par(mfrow = c(1, 1))
legend("bottom",
  legend = c("Target", "Model: mean", "Model: simulation"),
  col = c("black", "red", adjustcolor("indianred1", alpha = 0.2)),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = 3,
  bty = "n",
  pt.cex = 2,
  cex = 1.1,
  text.col = "black",
  horiz = T
)


cal_param_names <- names(calibration_params)

len <- length(cal_param_names)
calib_post <- cal_results[, (dim(cal_results)[2] - len + 1):dim(cal_results)[2]]
# REVIEWED par -> param
par <- rep(colnames(calib_post), 2)
# REVIEWED pe = point estimate, lend = lower end, uend = upper end
case <- c(
  rep("prior", length(cal_param_names)),
  rep("posterior", length(cal_param_names))
)
pe <- rep(0, length(cal_param_names) * 2)
lend <- rep(0, length(cal_param_names) * 2)
uend <- rep(0, length(cal_param_names) * 2)

plot_data <- data.frame(
  par = par, case = case, pe = pe, lend = lend, uend = uend
)
calib_par <- read.xlsx("Inputs/MasterTable.xlsx", sheet = "CalibPar")
plot_data$pe[plot_data$case == "prior"] <- calib_par$pe
plot_data$lend[plot_data$case == "prior"] <- calib_par$pe - calib_par$lower
plot_data$uend[plot_data$case == "prior"] <- calib_par$upper - calib_par$pe

plot_data$pe[plot_data$case == "posterior"] <- apply(calib_post, 2, median)

lend <- apply(calib_post, 2, median) -
  apply(calib_post, 2, quantile, probs = 0.025)
plot_data$lend[plot_data$case == "posterior"] <- lend

uend <- apply(calib_post, 2, quantile, probs = 0.975) -
  apply(calib_post, 2, median)
plot_data$uend[plot_data$case == "posterior"] <- uend

plot_data$case <- factor(plot_data$case, levels = c("prior", "posterior"))

p <- ggplot(plot_data, aes(x = case, y = pe, color = case)) +
  geom_point() +
  geom_errorbar(aes(
    ymin = pe - lend,
    ymax = pe + uend),
    width = .2,
    position = position_dodge(0.05)
  ) +
  facet_wrap(~par, scales = "free_y", ncol = 5) +
  labs(y = "Value", x = "") +
  theme_bw()

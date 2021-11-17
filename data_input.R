#' Data import for PROFOUND
#'
#' @description
#' `data_input()` reads in empirical data to inform the microsimulation.
#'
#' @param main_table The file location for the main empirical data table
#'    (.xlsx format)
#'
#' @returns
#' empirical data for the microsimulation
#'

## Define input excel file for empirical data##
library(openxlsx)

data_input <- function(inputs) {
  # hashed environment to hold params
  params <- new.env(hash = TRUE)
  workbook <- loadWorkbook(inputs$main_table)
  
  # Parameters for initial cohort ---------------------------------------------
  ppl_data <- openxlsx::read.xlsx(workbook, sheet = "InitialPop")
  params$ppl_size <- round(with(
    ppl_data,
    pe[par == "pop.size"]) * with(ppl_data,
    pe[par == "prop.12older"]),
    0
  )

  params$oud_prev <- with(ppl_data, pe[par == "prev.oud"])
  # params$nodu_m_prev <- with(ppl_data, pe[par == "prev.NODU" & sex == "m"])
  # params$nodu_f_prev <- with(ppl_data, pe[par == "prev.NODU" & sex == "f"])
  
  params$demographic <- read.xlsx(workbook, sheet = "Demographic")
  params$demo.mx <- data.matrix(params$demographic[, 4:ncol(params$demographic)])
  
  regions <- colnames(params$demographic)[-c(1:3)] # region names (city/town)
  # change . to space in region names
  params$regions <- gsub("\\.", " ", regions)
  # p_region <- read.xlsx(workbook, sheet = "region_prob")
  # rownames(p_region) <- p_region$
  params$races <- params$demographic$race
  params$ages <- params$demographic$age
  
  params$oud_demo <- openxlsx::read.xlsx(workbook, sheet = "OUDPrevNSDUH")$pe
  params$stim_demo <- openxlsx::read.xlsx(workbook, sheet = "StimPrevNSDUH")$pe
  
  # Opioid use patterns -------------------------------------------------------
  opioid_pattern <- read.xlsx(workbook, sheet = "OpioidPattern")
  params$init_il_m <- with(opioid_pattern, pe[par == "ini.il" & sex == "m"])
  params$init_il_f <- with(opioid_pattern, pe[par == "ini.il" & sex == "f"])
  # % of high-risk among illicit opioid ppl
  params$init_il_hr_m <- with(opioid_pattern, pe[par == "ini.il.hr" & sex == "m"])
  params$init_il_hr_f <- with(opioid_pattern, pe[par == "ini.il.hr" & sex == "f"])
  # % inactive
  params$init_inactive <- with(opioid_pattern, pe[par == "ini.inactive"])
  # TO_REVIEW why is this params instead of initials?
  params$init_oud_fx <- with(opioid_pattern, pe[par == "ini.oud.fx"])
  params$fx_growth <- with(opioid_pattern, pe[par == "fx_growth"])
  # initial probabilities that agent has previously overdosed
  params$init_everod_rx <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "preb"]
  )
  params$init_everod_il_lr <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "il.lr"]
  )
  params$init_everod_il_hr <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "il.hr"]
  )
  
  stimulant_pattern <- read.xlsx(workbook, sheet = "StimulantPattern")
  params$init_noud_fx <- with(stimulant_pattern, pe[par == "ini.NOUD.fx"])
  params$init_everod_sti <- with(stimulant_pattern, pe[par == "ini.everod"])
  
  #   init_everod_rx = init_everod_rx, init_everod_il_lr = init_everod_il_lr,
  #   init_everod_il_hr = init_everod_il_hr,
  #   ini_everod_sti = ini_everod_sti
  # )
  
  # Microsimulation parameters ------------------------------------------------
  # Life table
  annual_mortality_base <- read.xlsx(workbook, sheet = "LifeTable")$pe
  params$mortality_base <- 1 - (1 - annual_mortality_base / 1000000) ^ (1 / 12)
  annual_mortality_drug <- read.xlsx(workbook, sheet = "LifeTable")$drug
  params$mortality_drug <- 1 - (1 - annual_mortality_drug / 1000000) ^ (1 / 12)
  
  mor_gp <- read.xlsx(workbook, sheet = "LifeTable")$age
  rm(list = c("annual_mortality_base", "annual_mortality_drug"))
  
  # Overdose risk
  overdose_risk <- read.xlsx(workbook, sheet = "OverdoseRisk")
  # TO_REVIEW: unused outside calibration/validated city?
  params$od_rx_sub <- with(overdose_risk, pe[par == "od.preb.sub"])
  params$od_il_lr_sub <- with(overdose_risk, pe[par == "od.il.lr.sub"])
  params$od_nodu_sub <- with(overdose_risk, pe[par == "od.NODU.sub"])
  params$multi_hr <- with(overdose_risk, pe[par == "multi.hr"])
  params$multi_fx <- with(overdose_risk, pe[par == "multi.fx"])
  params$multi_relap <- with(overdose_risk, pe[par == "multi.relap"])
  params$multi_sub <- with(overdose_risk, pe[par == "multi.sub"])
  params$multi_NODU_fx <- with(overdose_risk, pe[par == "multi.NODU.fx"])
  
  # Drug use state transition probability
  p_transition <- read.xlsx(workbook, sheet = "TransProb")
  params$p_rx2il_lr <- with(p_transition, pe[par == "p.preb2il.lr"])
  params$p_rx2inact <- with(p_transition, pe[par == "p.preb2inact"])
  params$p_il_lr2il_hr <- with(p_transition, pe[par == "p.il.lr2il.hr"])
  params$p_il_lr2inact <- with(p_transition, pe[par == "p.il.lr2inact"])
  params$p_il_hr2il_lr <- with(p_transition, pe[par == "p.il.hr2il.lr"])
  params$p_il_hr2inact <- with(p_transition, pe[par == "p.il.hr2inact"])
  params$p_inact2relap <- with(p_transition, pe[par == "p.inact2relap"])
  
  
  # Decision tree parameters -------------------------------------------------
  if (inputs$strat == "regional") {
    od_loc_priv <- read.xlsx(workbook, sheet = "ODSettingEMS(sp)")$private
    od_loc_pub <- read.xlsx(workbook, sheet = "ODSettingEMS(sp)")$public
    od_loc <- rbind(od_loc_priv, od_loc_pub)
  } else if (inputs$strat == "overall") {
    od_loc_priv <- read.xlsx(workbook, sheet = "ODSettingEMS")$private
    od_loc_pub <- read.xlsx(workbook, sheet = "ODSettingEMS")$public
    od_loc <- rbind(
      rep(od_loc_priv, length(params$regions)),
      rep(od_loc_pub, length(params$regions))
    )
  } else {
    stop("Please choose a valid strategy from 'regional' or 'overall'")
  }
  
  rownames(od_loc) <- c("priv", "pub")
  colnames(od_loc) <- params$regions
  params$od_loc <- od_loc
  
  decision_tree <- read.xlsx(workbook, sheet = "DecisionTree")
  params$od_loc_pub <- with(decision_tree, pe[par == "OD_loc_pub"])
  params$od_wit_pub <- with(decision_tree, pe[par == "OD_wit_pub"])
  params$rr_od_wit_priv <- with(decision_tree, pe[par == "rr_OD_wit_priv"])
  params$od_wit_priv <- params$od_wit_pub * params$rr_od_wit_priv
  params$od_911_pub <- with(decision_tree, pe[par == "OD_911_pub"])
  params$rr_od_911_priv <- with(decision_tree, pe[par == "rr_OD_911_priv"])
  params$od_911_priv <- params$OD_911_pub * params$rr_OD_911_priv
  params$od_hosp <- with(decision_tree, pe[par == "OD_hosp"])
  params$od_cess <- with(decision_tree, pe[par == "OD_cess"])
  # 
  mortality <- read.xlsx(workbook, sheet = "Mortality")
  params$mor_bl <- with(mortality, pe[par == "mor_bl"])
  params$mortality_nx <- with(mortality, pe[par == "mortality_nx"])
  params$rr_mor_ems <- with(mortality, pe[par == "rr_mor_EMS"])
  
  # Naloxone kit parameters ---------------------------------------------------
  nlx_kit <- read.xlsx(workbook, sheet = "NxKit")
  params$r_loss_exp <- 1 / with(nlx_kit, pe[par == "LossExp"])
  params$low2priv <- with(nlx_kit, pe[par == "Low2Priv"])
  params$nlx_adj <- with(nlx_kit, pe[par == "nlx.adj"])
  params$cap <- with(nlx_kit, pe[par == "cap"])
  nlx_data_oend <- read.xlsx(workbook, sheet = "NxDataOEND")
  nx_oend <- array(
    0,
    dim = c(
      length(unique(nlx_data_oend$year)),
      length(unique(nlx_data_oend$risk)), length(params$regions))
  )
  dimnames(nx_oend)[[1]] <- unique(nlx_data_oend$year)
  dimnames(nx_oend)[[2]] <- unique(nlx_data_oend$risk)
  dimnames(nx_oend)[[3]] <- params$regions
  for (i in seq_len(length(unique(nlx_data_oend$year)))) {
    nx_oend[i, , ] <- data.matrix(
      subset(nlx_data_oend, year == unique(nlx_data_oend$year)[i])[-c(1, 2)]
    )
  }
  params$nx_oend <- nx_oend
  params$nlx_data_pharm <- read.xlsx(workbook, sheet = "NxDataPharm")
  
  # Cost parameters ------------------------------------------------------------
  cost <- read.xlsx(workbook, sheet = "Cost")
  params$c_rx <- with(cost, pe[par == "c.preb"])
  params$c_il_lr <- with(cost, pe[par == "c.il.lr"])
  params$c_il_hr <- with(cost, pe[par == "c.il.hr"])
  params$c_inact <- with(cost, pe[par == "c.inact"])
  params$c.NODU <- with(cost, pe[par == "c.NODU"])
  params$c_nlx_kit <- with(cost, pe[par == "c.nlx.kit"])
  params$c_nlx_dtb <- with(cost, pe[par == "c_nlx_dtb"])
  
  c_relap <- numeric(0)
  c_relap["rx"] <- (params$c_rx + params$c_inact) / 2
  c_relap["il_lr"] <- (params$c_il_lr + params$c_inact) / 2
  c_relap["il_hr"] <- (params$c_il_hr + params$c_inact) / 2
  params$c_relap <- c_relap
  
  params$ems_cost <- with(cost, pe[par == "c.EMS"])
  params$hosp_cost <- with(cost, pe[par == "c.hospcare"])
  
  # Overdose probability matrix (per month) ------------------------------
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["rx", "subs"] <- params$od_rx_sub
  overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
  overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
  overdose_probs["NODU", "subs"] <- params$od_nodu_sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
  params$overdose_probs <- overdose_probs
  
  # Baseline mortality excluding overdose (per month) ---------------
  params$mortality_probs <- matrix(0, nrow = 2, ncol = length(mor_gp))
  rownames(params$mortality_probs) <- c("bg", "drug")
  colnames(params$mortality_probs) <- mor_gp
  params$mortality_probs["bg", ] <- params$mortality_base
  params$mortality_probs["drug", ] <- params$mortality_drug
  
  return(params)
}

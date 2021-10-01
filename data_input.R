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

library(openxlsx)

data_input <- function(main_table) {
  params <- new.env(hash = TRUE)
  workbook <- loadWorkbook(main_table)  # load data table

  # Parameters for initial cohort --------------------------
  ppl_data <- read.xlsx(workbook, sheet = "InitialPop")
  ppl_size <- round(with(
    ppl_data,
    pe[par == "ppl_size"]) * with(ppl_data,
    pe[par == "prop.12older"]),
    0
  )
  oud_prev <- with(ppl_data, pe[par == "prev.oud"])
  nodu_m_prev <- with(ppl_data, pe[par == "prev.NODU" & sex == "m"])
  nodu_f_prev <- with(ppl_data, pe[par == "prev.NODU" & sex == "f"])

  params$demographic <- read.xlsx(workbook, sheet = "Demographic")

  regions <- colnames(params$demographic)[-c(1:3)] # region names (city/town)
  regions <- gsub("\\.", " ", regions)
  p_region <- read.xlsx(workbook, sheet = "region_prob")
  rownames(p_region) <- p_region$region

  params$regions <- regions
  params$races <- params$demographic$race
  params$ages <- params$demographic$age

  oud_demo <- read.xlsx(workbook, sheet = "OUDPrevNSDUH")
  stim_demo <- read.xlsx(workbook, sheet = "StimPrevNSDUH")

  opioid_pattern <- read.xlsx(workbook, sheet = "OpioidPattern")
  # % of illicite opioid use among OUD ppl
  init_il_m <- with(opioid_pattern, pe[par == "ini.il" & sex == "m"])
  init_il_f <- with(opioid_pattern, pe[par == "ini.il" & sex == "f"])
  # % of high-risk among illicit opioid ppl
  init_il_hr_m <- with(opioid_pattern, pe[par == "ini.il.hr" & sex == "m"])
  init_il_hr_f <- with(opioid_pattern, pe[par == "ini.il.hr" & sex == "f"])
  # % inactive
  init_inactive <- with(opioid_pattern, pe[par == "init_inactive"])
  # TO_REVIEW what is this %?
  params$init_oud_fx <- with(opioid_pattern, pe[par == "init_oud_fx"])
  # TO_REVEW what does gw stand for
  fx_growth <- with(opioid_pattern, pe[par == "fx_growth"])
  # initial probabilities that agent has previously overdosed
  init_everod_rx <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "rx"]
  )
  init_everod_il_lr <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "il.lr"]
  )
  init_everod_il_hr <- with(
    opioid_pattern,
    pe[par == "ini.everod" & group == "il.hr"]
  )
  
  stimulant_pattern <- read.xlsx(workbook, sheet = "StimulantPattern")
  params$init_noud_fx <- with(stimulant_pattern, pe[par == "ini.NOUD.fx"])
  ini_everod_sti <- with(stimulant_pattern, pe[par == "ini.everod"])

  params$initials <- list(
    ppl_size = ppl_size, oud_prev = oud_prev, nodu_m_prev = nodu_m_prev,
    nodu_f_prev = nodu_f_prev, regions = regions,
    oud_demo = oud_demo, stim_demo = stim_demo, init_il_m = init_il_m,
    init_il_f = init_il_f, init_il_hr_m = init_il_hr_m,
    init_il_hr_f = init_il_hr_f, init_inactive = init_inactive,
    init_everod_rx = init_everod_rx, init_everod_il_lr = init_everod_il_lr,
    init_everod_il_hr = init_everod_il_hr,
    ini_everod_sti = ini_everod_sti, p_region = p_region
  )

  params$fx_growth <- fx_growth

  ## Parameters for microsimulation ##
  # life table: for mortality
  annual_mortality_base <- read.xlsx(workbook, sheet = "LifeTable")$pe
  params$mortality_base <- 1 - (1 - annual_mortality_base / 1000000)^ (1 / 12)
  annual_mortality_drug <- read.xlsx(workbook, sheet = "LifeTable")$drug
  params$mortality_drug <- 1 - (1 - annual_mortality_drug / 1000000)^ (1 / 12)

  mor.gp <- read.xlsx(workbook, sheet = "LifeTable")$age
  rm(list = c("annual_mortality_base", "annual_mortality_drug"))
  # risk of overdose
  overdose_risk <- read.xlsx(workbook, sheet = "OverdoseRisk")
  params$od_rx_sub <- with(overdose_risk, pe[par == "od.rx.sub"])
  params$od_il_lr_sub <- with(overdose_risk, pe[par == "od.il.lr.sub"])
  params$od_NODU_sub <- with(overdose_risk, pe[par == "od.NODU.sub"])
  params$multi_hr <- with(overdose_risk, pe[par == "multi.hr"])
  params$multi_fx <- with(overdose_risk, pe[par == "multi.fx"])
  params$multi_relap <- with(overdose_risk, pe[par == "multi.relap"])
  params$multi_sub <- with(overdose_risk, pe[par == "multi.sub"])
  params$multi_NODU_fx <- with(overdose_risk, pe[par == "multi.NODU.fx"])
  # transition probability
  p_transition <- read.xlsx(workbook, sheet = "TransProb")
  params$p_rx2il_lr <- with(p_transition, pe[par == "p.rx2il.lr"])
  params$p.rx2inact <- with(p_transition, pe[par == "p.rx2inact"])
  params$p.il_lr2il_hr <- with(p_transition, pe[par == "p.il.lr2il.hr"])
  params$p.il_lr2inact <- with(p_transition, pe[par == "p.il.lr2inact"])
  params$p.il_hr2il_lr <- with(p_transition, pe[par == "p.il.hr2il.lr"])
  params$p.il_hr2inact <- with(p_transition, pe[par == "p.il.hr2inact"])
  params$p.inact2relap <- with(p_transition, pe[par == "p.inact2relap"])


  ## Parameters for decision tree ##
  if (exists("sw.EMS.ODloc")) {
    if (sw.EMS.ODloc == "sp") {
      OD_loc_priv <- read.xlsx(workbook, sheet = "ODSettingEMS(sp)")$private
      OD_loc_pub <- read.xlsx(workbook, sheet = "ODSettingEMS(sp)")$public
      OD_loc <- rbind(OD_loc_priv, OD_loc_pub)
    } else {
      OD_loc_priv <- read.xlsx(workbook, sheet = "ODSettingEMS")$private
      OD_loc_pub <- read.xlsx(workbook, sheet = "ODSettingEMS")$public
      OD_loc <- rbind(
        rep(OD_loc_priv, length(regions)),
        rep(OD_loc_pub, length(regions))
      )
    }
  } else {
    OD_loc_priv <- read.xlsx(workbook, sheet = "ODSettingEMS")$private
    OD_loc_pub <- read.xlsx(workbook, sheet = "ODSettingEMS")$public
    OD_loc <- rbind(
      rep(OD_loc_priv, length(regions)),
      rep(OD_loc_pub, length(regions))
    )
  }
  rownames(OD_loc) <- c("priv", "pub")
  colnames(OD_loc) <- regions
  params$OD_loc <- OD_loc

  DecisionTree <- read.xlsx(workbook, sheet = "DecisionTree")
  params$OD_wit_priv <- with(
    DecisionTree,
    pe[par == "OD_wit" & group == "priv"]
  )
  params$OD_wit_pub <- with(DecisionTree, pe[par == "OD_wit" & group == "pub"])
  params$OD_911_priv <- with(DecisionTree, pe[par == "OD_911_priv"])
  params$OD_911_pub_mul <- with(DecisionTree, pe[par == "OD_911_pub_mul"])
  params$OD_911_pub <- params$OD_911_priv * params$OD_911_pub_mul
  params$OD_hosp <- with(DecisionTree, pe[par == "OD_hosp"])
  params$OD_cess <- with(DecisionTree, pe[par == "OD_cess"])

  Mortality <- read.xlsx(workbook, sheet = "Mortality")
  params$mor_bl <- with(Mortality, pe[par == "mor_bl"])
  params$mortality_nx <- with(Mortality, pe[par == "mortality_nx"])
  params$rr_mor_EMS <- with(Mortality, pe[par == "rr_mor_EMS"])

  ## Parameters for naloxone kits ##
  NxKit <- read.xlsx(workbook, sheet = "NxKit")
  params$r.LossExp <- 1 / with(NxKit, pe[par == "LossExp"])
  params$Low2Priv <- with(NxKit, pe[par == "Low2Priv"])
  params$nlx.adj <- with(NxKit, pe[par == "nlx.adj"])
  params$cap <- with(NxKit, pe[par == "cap"])
  NxDataOEND <- read.xlsx(workbook, sheet = "NxDataOEND")
  NxOEND <- array(
    0,
    dim = c(
      length(unique(NxDataOEND$year)),
      length(unique(NxDataOEND$risk)), length(regions))
  )
  dimnames(NxOEND)[[1]] <- unique(NxDataOEND$year)
  dimnames(NxOEND)[[2]] <- unique(NxDataOEND$risk)
  dimnames(NxOEND)[[3]] <- regions
  for (i in 1:length(unique(NxDataOEND$year))) {
    NxOEND[i, , ] <- data.matrix(
      subset(NxDataOEND, year == unique(NxDataOEND$year)[i])[-c(1, 2)]
    )
  }
  params$NxOEND <- NxOEND
  params$NxDataPharm <- read.xlsx(workbook, sheet = "NxDataPharm")
  NxMvt <- data.matrix(read.xlsx(workbook, sheet = "NxMvt")[, -1])
  row.names(NxMvt) <- regions
  params$NxMvt <- NxMvt

  ## Parameters for cost ##
  Cost <- read.xlsx(workbook, sheet = "Cost")
  params$c.rx <- with(Cost, pe[par == "c.rx"])
  params$c.il_lr <- with(Cost, pe[par == "c.il.lr"])
  params$c.il_hr <- with(Cost, pe[par == "c.il.hr"])
  params$c.inact <- with(Cost, pe[par == "c.inact"])
  params$c.NODU <- with(Cost, pe[par == "c.NODU"])
  params$c.nlx.kit <- with(Cost, pe[par == "c.nlx.kit"])
  # REVIEWED database
  params$c.nlx.dtb <- with(Cost, pe[par == "c.nlx.dtb"])

  c.relap.v <- numeric(0)
  c.relap.v["rx"] <- (params$c.rx + params$c.inact) / 2
  c.relap.v["il_lr"] <- (params$c.il_lr + params$c.inact) / 2
  c.relap.v["il_hr"] <- (params$c.il_hr + params$c.inact) / 2
  params$c.relap.v <- c.relap.v

  params$c.EMS <- with(Cost, pe[par == "c.EMS"])
  params$c.hospcare <- with(Cost, pe[par == "c.hospcare"])

  # Overdose probability matrix (per month)
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["rx", "subs"] <- params$od_rx_sub
  overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
  overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
  overdose_probs["NODU", "subs"] <- params$od_NODU_sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
  params$overdose_probs <- overdose_probs

  # Baseline mortality excluding overdose (per month)
  params$mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
  rownames(params$mortality_probs) <- c("bg", "drug")
  colnames(params$mortality_probs) <- mor.gp
  params$mortality_probs["bg", ] <- params$mortality_base
  params$mortality_probs["drug", ] <- params$mortality_drug

  return(params)
}

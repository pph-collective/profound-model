########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Module for Data Input of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: April 12, 2021
#
#############################################################################

## Define input excel file for empirical data##
library(openxlsx)
params <- list()
WB <- loadWorkbook("Inputs/MasterTable.xlsx")

## Parameters for initial cohort ##
NxDataOEND <- read.xlsx(WB, sheet = "NxDataOEND")
v.region <- unique(NxDataOEND$catchment) 
params$v.region <- v.region
InitialPop <- read.xlsx(WB, sheet = "InitialPop")
ppl.size <- round(with(InitialPop, pe[par == "pop.size"]) * with(InitialPop, pe[par == "prop.12older"]), 0) # size of initial population 12 and older
prev.oud <- with(InitialPop, pe[par == "prev.oud"]) # prevalence of OUD (or at risk for OUD)
prev.NODU.m <- with(InitialPop, pe[par == "prev.NODU" & sex == "m"]) # prevalence of non-opioid drug use among males, in addition to OUD
prev.NODU.f <- with(InitialPop, pe[par == "prev.NODU" & sex == "f"]) # prevalence of non-opioid drug use among females, in addition to OUD

Demographic <- read.xlsx(WB, sheet = "Demographic")

demo.mx <- data.matrix(Demographic[, 4:ncol(Demographic)])
v.demo.race <- Demographic$race
v.demo.age <- Demographic$age

OUDDemo <- read.xlsx(WB, sheet = "OUDPrevNSDUH")$pe
StimDemo <- read.xlsx(WB, sheet = "StimPrevNSDUH")$pe

# REVIEWED init = initial, lr/hr lowrisk/highrisk, il = illegal, inact = inactive, gw = annual growth rate of fx exposure, preb = prescription, opioid_use_patterns
# why create all these variables and then put them into a dataframe? Wastes memory
OpioidPattern <- read.xlsx(WB, sheet = "OpioidPattern")
ini.il.m <- with(OpioidPattern, pe[par == "ini.il" & sex == "m"]) # % of illicite opioid use among OUD ppl
ini.il.f <- with(OpioidPattern, pe[par == "ini.il" & sex == "f"]) # % of illicite opioid use among OUD ppl
ini.il.hr.m <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "m"]) # % of high-risk among illicit opioid ppl
ini.il.hr.f <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "f"]) # % of high-risk among illicit opioid ppl
ini.inactive <- setNames(numeric(3), c("white", "black", "hisp"))  ##NEWLY MODIFIED##
ini.inactive['white'] <- with(OpioidPattern, pe[par == "ini.inactive"])
ini.inactive['black'] <- ini.inactive['white'] * with(OpioidPattern, pe[par == "ini.inactive.b2w"])
ini.inactive['hisp'] <- ini.inactive['white'] * with(OpioidPattern, pe[par == "ini.inactive.h2w"])
ini.oud.fx <- setNames(numeric(3), c("white", "black", "hisp")) ##NEWLY MODIFIED##
ini.oud.fx['white'] <- with(OpioidPattern, pe[par == "ini.oud.fx" & group == "white"]) 
ini.oud.fx['black'] <- with(OpioidPattern, pe[par == "ini.oud.fx" & group == "black"])
ini.oud.fx['hisp'] <- with(OpioidPattern, pe[par == "ini.oud.fx" & group == "hisp"]) ##

gw.fx <- setNames(numeric(3), c("white", "black", "hisp"))
gw.fx['white'] <- with(OpioidPattern, pe[par == "gw.fx" & group == "white"])
gw.fx['black'] <- with(OpioidPattern, pe[par == "gw.fx" & group == "black"])
gw.fx['hisp']  <- with(OpioidPattern, pe[par == "gw.fx" & group == "hisp"])
ini.everod.preb <- with(OpioidPattern, pe[par == "ini.everod" & group == "preb"])
ini.everod.il.lr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.lr"])
ini.everod.il.hr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.hr"])
out.prebopioid <- with(OpioidPattern, pe[par == "out.prebopioid"])

# REVIEWED stimulant_use_patterns
StimulantPattern <- read.xlsx(WB, sheet = "StimulantPattern")
ini.NODU.fx <- setNames(numeric(3), c("white", "black", "hisp"))
ini.NODU.fx[] <- with(StimulantPattern, pe[par == "ini.NODU.fx"])
gw.NODU.fx.ab.yr <- with(StimulantPattern, pe[par == "gw.NODU.fx.ab.yr"])  ##NEWLY ADDED##
covid.NODU.fx <- setNames(numeric(3), c("white", "black", "hisp"))
covid.NODU.fx['white'] <- with(StimulantPattern, pe[par == "covid.NODU.fx" & group == "white"])
covid.NODU.fx['black'] <- with(StimulantPattern, pe[par == "covid.NODU.fx" & group == "black"])
covid.NODU.fx['hisp'] <- with(StimulantPattern, pe[par == "covid.NODU.fx" & group == "hisp"])  ##
ini.everod.sti <- with(StimulantPattern, pe[par == "ini.everod"])

# REVIEWED things used in initilization functions \ see if i can add this without above initials = initial values?
initials <- list(
  ppl.size = ppl.size, prev.oud = prev.oud, prev.NODU.m = prev.NODU.m, prev.NODU.f = prev.NODU.f, demo.mx = demo.mx, v.region = v.region, OUDDemo = OUDDemo, StimDemo = StimDemo,
  ini.il.m = ini.il.m, ini.il.f = ini.il.f, ini.il.hr.m = ini.il.hr.m, ini.il.hr.f = ini.il.hr.f, ini.inactive = ini.inactive,
  # init_oud_fx = init_oud_fx, ini.NODU.fx = ini.NODU.fx,
  ini.everod.preb = ini.everod.preb, ini.everod.il.lr = ini.everod.il.lr, ini.everod.il.hr = ini.everod.il.hr
)

# REVIEWED add these straight to list, no intermediary
params$ini.oud.fx <- ini.oud.fx
params$gw.fx <- gw.fx
params$ini.NODU.fx <- ini.NODU.fx
# params$gw.NODU.fx.ab.yr <- gw.NODU.fx.ab.yr
params$covid.NODU.fx <- covid.NODU.fx
params$out.prebopioid <- out.prebopioid

## Parameters for microsimulation ##
# life table: for mortality
mor.bg.y <- cbind(read.xlsx(WB, sheet = "LifeTable")$white, read.xlsx(WB, sheet = "LifeTable")$black, read.xlsx(WB, sheet = "LifeTable")$hispanic)
rownames(mor.bg.y) <- read.xlsx(WB, sheet = "LifeTable")$age
colnames(mor.bg.y) <- c("white", "black", "hisp")
params$mor.bg <- 1 - (1 - mor.bg.y / 100000)^(1 / 12)
mor.drug.y <- cbind(read.xlsx(WB, sheet = "LifeTable")$drug_white, read.xlsx(WB, sheet = "LifeTable")$drug_black, read.xlsx(WB, sheet = "LifeTable")$drug_hispanic)
rownames(mor.drug.y) <- rownames(mor.bg.y)
colnames(mor.drug.y) <- colnames(mor.bg.y)
params$mor.drug <- 1 - (1 - mor.drug.y / 100000)^(1 / 12)
# REVIEWED gp = general population
mor.gp <- read.xlsx(WB, sheet = "LifeTable")$age
rm(list = c("mor.bg.y"))
# risk of overdose
OverdoseRisk <- read.xlsx(WB, sheet = "OverdoseRisk")
params$od.preb.sub <- with(OverdoseRisk, pe[par == "od.preb.sub"])
params$od.il.lr.sub <- with(OverdoseRisk, pe[par == "od.il.lr.sub"])
params$od.NODU.sub <- with(OverdoseRisk, pe[par == "od.NODU.sub"])
params$multi.hr <- with(OverdoseRisk, pe[par == "multi.hr"])
params$multi.fx <- with(OverdoseRisk, pe[par == "multi.fx"])
params$multi.relap <- with(OverdoseRisk, pe[par == "multi.relap"])
params$multi.sub <- with(OverdoseRisk, pe[par == "multi.sub"])
params$multi.NODU.fx <- with(OverdoseRisk, pe[par == "multi.NODU.fx"])
# transition probability
TransProb <- read.xlsx(WB, sheet = "TransProb")
params$p.preb2il.lr <- with(TransProb, pe[par == "p.preb2il.lr"])
params$p.il.lr2il.hr <- with(TransProb, pe[par == "p.il.lr2il.hr"])
params$p.il.hr2il.lr <- with(TransProb, pe[par == "p.il.hr2il.lr"])
params$p.inact2relap <- with(TransProb, pe[par == "p.inact2relap"])
params$gw.m.2inact <- with(TransProb, pe[par == "gw.m.2inact"])
p.preb2inact.ini <- p.il.lr2inact.ini <- p.il.hr2inact.ini <- setNames(numeric(3), c("white", "black", "hisp"))  ## NEWLY MODIFIED
cov.trmt.b2w <- with(TransProb, pe[par == "cov.trmt.b2w"])
cov.trmt.h2w <- with(TransProb, pe[par == "cov.trmt.h2w"])
p.preb2inact.ini['white']  <- with(TransProb, pe[par == "p.preb2inact.ini"]) 
p.preb2inact.ini['black']  <- cov.trmt.b2w * p.preb2inact.ini['white'] * params$p.inact2relap /(p.preb2inact.ini['white'] + params$p.inact2relap - cov.trmt.b2w * p.preb2inact.ini['white'])
p.preb2inact.ini['hisp']   <- cov.trmt.h2w * p.preb2inact.ini['white'] * params$p.inact2relap /(p.preb2inact.ini['white'] + params$p.inact2relap - cov.trmt.h2w * p.preb2inact.ini['white'])
p.il.lr2inact.ini['white'] <- with(TransProb, pe[par == "p.il.lr2inact.ini"])
p.il.lr2inact.ini['black']  <- cov.trmt.b2w * p.il.lr2inact.ini['white'] * params$p.inact2relap /(p.il.lr2inact.ini['white'] + params$p.inact2relap - cov.trmt.b2w * p.il.lr2inact.ini['white'])
p.il.lr2inact.ini['hisp']   <- cov.trmt.h2w * p.il.lr2inact.ini['white'] * params$p.inact2relap /(p.il.lr2inact.ini['white'] + params$p.inact2relap - cov.trmt.h2w * p.il.lr2inact.ini['white'])
p.il.hr2inact.ini['white'] <- with(TransProb, pe[par == "p.il.hr2inact.ini"])
p.il.hr2inact.ini['black']  <- cov.trmt.b2w * p.il.hr2inact.ini['white'] * params$p.inact2relap /(p.il.hr2inact.ini['white'] + params$p.inact2relap - cov.trmt.b2w * p.il.hr2inact.ini['white'])
p.il.hr2inact.ini['hisp']   <- cov.trmt.h2w * p.il.hr2inact.ini['white'] * params$p.inact2relap /(p.il.hr2inact.ini['white'] + params$p.inact2relap - cov.trmt.h2w * p.il.hr2inact.ini['white'])
params$p.preb2inact.ini  <- p.preb2inact.ini
params$p.il.lr2inact.ini <- p.il.lr2inact.ini
params$p.il.hr2inact.ini <- p.il.hr2inact.ini
params$covid.rd.2inact.black <- with(TransProb, pe[par == "covid.rd.2inact.black"])
params$covid.rd.2inact.hisp <- with(TransProb, pe[par == "covid.rd.2inact.hisp"]) ##

## Parameters for decision tree ##
DecisionTree <- read.xlsx(WB, sheet = "DecisionTree")
params$OD_loc_pub <- with(DecisionTree, pe[par == "OD_loc_pub"])
params$OD_wit_pub <- with(DecisionTree, pe[par == "OD_wit_pub"])
params$rr_OD_wit_priv <- with(DecisionTree, pe[par == "rr_OD_wit_priv"])
params$OD_wit_priv <- params$OD_wit_pub * params$rr_OD_wit_priv
params$OD_911_pub <- with(DecisionTree, pe[par == "OD_911_pub"])
params$rr_OD_911_priv <- with(DecisionTree, pe[par == "rr_OD_911_priv"])
params$OD_911_priv <- params$OD_911_pub * params$rr_OD_911_priv
params$OD_hosp <- with(DecisionTree, pe[par == "OD_hosp"])
params$OD_cess <- setNames(numeric(3), c("white", "black", "hisp"))  ## NEWLY MODIFIED
params$OD_cess['white'] <- with(DecisionTree, pe[par == "OD_cess"])
params$OD_cess['black'] <- params$OD_cess['white'] * with(DecisionTree, pe[par == "rr_OD_cess_black"])
params$OD_cess['hisp']  <- params$OD_cess['white'] * with(DecisionTree, pe[par == "rr_OD_cess_hisp"])

Mortality <- read.xlsx(WB, sheet = "Mortality")
params$mor_bl <- with(Mortality, pe[par == "mor_bl"])
params$mor_nx <- with(Mortality, pe[par == "mor_nx"])
params$rr_mor_EMS <- with(Mortality, pe[par == "rr_mor_EMS"])

## Parameters for naloxone kits ##
NxKit <- read.xlsx(WB, sheet = "NxKit")
params$r.LossExp <- 1 / with(NxKit, pe[par == "LossExp"])
params$nlx.adj <- with(NxKit, pe[par == "nlx.adj"])
params$cap <- with(NxKit, pe[par == "cap"])
params$eff.pharNlx <- with(NxKit, pe[par == "eff.pharNlx"])
NxDataOEND <- read.xlsx(WB, sheet = "NxDataOEND")
NxOENDtotal <- read.xlsx(WB, sheet = "NxOENDtotal")
v.region <- unique(NxDataOEND$catchment) 
params$v.region <- v.region
NxOEND.array <- array(data = 0, dim = c(length(unique(NxDataOEND$year)), length(unique(NxDataOEND$catchment)), ncol(NxDataOEND)-2))
dimnames(NxOEND.array)[[1]] <- unique(NxDataOEND$year)
dimnames(NxOEND.array)[[2]] <- v.region
dimnames(NxOEND.array)[[3]] <- c("white", "black", "hisp", "other")
for (i in 1:length(unique(NxDataOEND$year))){
  NxOEND.array[i,,] <- as.matrix(NxDataOEND[NxDataOEND$year== unique(NxDataOEND$year)[i], -c(1,2)])
  NxOEND.array[i,,] <- round(NxOEND.array[i,,] * (NxOENDtotal$pe[NxOENDtotal$year == unique(NxDataOEND$year)[i]] / sum(NxOEND.array[i,,])), 0)
}
NxOEND.array <- abind(round(NxOEND.array["2016",,]*(NxOENDtotal$pe[NxOENDtotal$year == 2015] / NxOENDtotal$pe[NxOENDtotal$year == 2016]), 0), NxOEND.array, along =1)
dimnames(NxOEND.array)[[1]][1] <- 2015
params$NxOEND.array <- NxOEND.array
params$NxDataPharm <- read.xlsx(WB, sheet = "NxDataPharm")

## Parameters for cost ##
Cost <- read.xlsx(WB, sheet = "Cost")
params$c.preb <- with(Cost, pe[par == "c.preb"])
params$c.il.lr <- with(Cost, pe[par == "c.il.lr"])
params$c.il.hr <- with(Cost, pe[par == "c.il.hr"])
params$c.inact <- with(Cost, pe[par == "c.inact"])
params$c.NODU <- with(Cost, pe[par == "c.NODU"])
params$c.nlx.kit <- with(Cost, pe[par == "c.nlx.kit"])
params$c.nlx.dtb <- with(Cost, pe[par == "c.nlx.dtb"])
c.relap.v <- numeric(0) # cost of remaining one cycle: relapsed, as the average of inactive and prior state
c.relap.v["preb"] <- (params$c.preb + params$c.inact) / 2
c.relap.v["il.lr"] <- (params$c.il.lr + params$c.inact) / 2
c.relap.v["il.hr"] <- (params$c.il.hr + params$c.inact) / 2
params$c.relap.v <- c.relap.v
params$c.EMS <- with(Cost, pe[par == "c.EMS"])
params$c.hospcare <- with(Cost, pe[par == "c.hospcare"])

## Parameters for QALYs ##
QALY <- read.xlsx(WB, sheet = "QALY")
params$q.il.lr <- with(QALY, pe[par == "q.il.lr"])
params$q.il.hr <- with(QALY, pe[par == "q.il.hr"])
params$q.preb <- with(QALY, pe[par == "q.preb"])
params$q.NODU <- with(QALY, pe[par == "q.NODU"])
params$q.inact <- with(QALY, pe[par == "q.inact"])
params$q.death <- with(QALY, pe[par == "q.death"])
q.relap.v <- numeric(0) # cost of remaining one cycle: relapsed, as the average of inactive and prior state
q.relap.v["preb"] <- (params$q.preb + params$q.inact) / 2
q.relap.v["il.lr"] <- (params$q.il.lr + params$q.inact) / 2
q.relap.v["il.hr"] <- (params$q.il.hr + params$q.inact) / 2
params$q.relap.v <- q.relap.v

# Overdose probability matrix (per month)
overdose_probs <- matrix(0, nrow = 4, ncol = 2)
rownames(overdose_probs) <- c("preb", "il.lr", "il.hr", "NODU")
colnames(overdose_probs) <- c("first", "subs")
overdose_probs["preb", "subs"] <- params$od.preb.sub
overdose_probs["il.lr", "subs"] <- params$od.il.lr.sub
overdose_probs["il.hr", "subs"] <- params$od.il.lr.sub * params$multi.hr
overdose_probs["NODU", "subs"] <- params$od.NODU.sub
overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi.sub
params$overdose_probs <- overdose_probs

# Baseline mortality excluding overdose (per month)
mortality_probs <- list(mor.bg = params$mor.bg, mor.drug = params$mor.drug)
params$mortality_probs <- mortality_probs

params$ODdeath2020 <- subset(read.xlsx(WB, sheet = "Target"), par == "ODdeaths" & year == 2020)[,-c(1:2)]
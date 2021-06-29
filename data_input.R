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
###################        DATA INPUT function      #########################
#############################################################################
## Define input excel file ##
vparameters   <- list()
WB            <- loadWorkbook("Inputs/MasterTable.xlsx")

## Parameters for initial cohort ##
InitialPop    <- read.xlsx(WB, sheet="InitialPop")
pop.size      <- round(with(InitialPop, pe[par == "pop.size"]) * with(InitialPop, pe[par == "prop.12older"]), 0)   #size of initial population 12 and older
prev.oud      <- with(InitialPop, pe[par == "prev.oud"])                 #prevalence of OUD (or at risk for OUD)
prev.NODU.m   <- with(InitialPop, pe[par == "prev.NODU" & sex == "m"])   #prevalence of non-opioid drug use among males, in addition to OUD
prev.NODU.f   <- with(InitialPop, pe[par == "prev.NODU" & sex == "f"])   #prevalence of non-opioid drug use among females, in addition to OUD

Demographic   <- read.xlsx(WB, sheet="Demographic")
demo.mx       <- data.matrix(Demographic[ , 4:ncol(Demographic)])
v.rgn         <- colnames(Demographic)[-c(1:3)]                          #region names (city/town)
vparameters$v.rgn <- v.rgn
v.demo.sex    <- Demographic$sex
v.demo.race   <- Demographic$race
v.demo.age    <- Demographic$age

OUDDemo       <- read.xlsx(WB, sheet="OUDPrevNSDUH")$pe
StimDemo      <- read.xlsx(WB, sheet="StimPrevNSDUH")$pe

OpioidPattern <- read.xlsx(WB, sheet="OpioidPattern")
ini.il.m      <- with(OpioidPattern, pe[par == "ini.il" & sex == "m"])   # % of illicite opioid use among OUD pop
ini.il.f      <- with(OpioidPattern, pe[par == "ini.il" & sex == "f"])   # % of illicite opioid use among OUD pop
ini.il.hr.m   <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "m"])# % of high-risk among illicit opioid pop
ini.il.hr.f   <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "f"])# % of high-risk among illicit opioid pop
ini.inact     <- with(OpioidPattern, pe[par == "ini.inact"])
ini.OUD.fx    <- with(OpioidPattern, pe[par == "ini.OUD.fx"])
gw.fx         <- with(OpioidPattern, pe[par == "gw.fx"])
ini.everod.preb  <- with(OpioidPattern, pe[par == "ini.everod" & group == "preb"])
ini.everod.il.lr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.lr"])
ini.everod.il.hr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.hr"])

StimulantPattern <- read.xlsx(WB, sheet="StimulantPattern")
ini.NOUD.fx      <- with(StimulantPattern, pe[par == "ini.NOUD.fx"])
ini.everod.sti   <- with(StimulantPattern, pe[par == "ini.everod"])

initials      <- list(pop.size = pop.size, prev.oud = prev.oud, prev.NODU.m = prev.NODU.m, prev.NODU.f = prev.NODU.f, demo.mx = demo.mx, v.rgn = v.rgn, OUDDemo = OUDDemo, StimDemo = StimDemo, 
                      ini.il.m = ini.il.m, ini.il.f = ini.il.f, ini.il.hr.m = ini.il.hr.m, ini.il.hr.f = ini.il.hr.f, ini.inact = ini.inact,
                      # ini.OUD.fx = ini.OUD.fx, ini.NOUD.fx = ini.NOUD.fx,
                      ini.everod.preb = ini.everod.preb, ini.everod.il.lr = ini.everod.il.lr, ini.everod.il.hr = ini.everod.il.hr)

vparameters$ini.OUD.fx  <- ini.OUD.fx
vparameters$gw.fx       <- gw.fx
vparameters$ini.NOUD.fx <- ini.NOUD.fx

## Parameters for microsimulation ##
#life table: for mortality
mor.bg.y                  <- read.xlsx(WB, sheet="LifeTable")$pe
vparameters$mor.bg        <- 1-(1-mor.bg.y/1000000)^(1/12)
mor.drug.y                <- read.xlsx(WB, sheet="LifeTable")$drug
vparameters$mor.drug      <- 1-(1-mor.drug.y/1000000)^(1/12)
mor.gp                    <- read.xlsx(WB, sheet="LifeTable")$age
#risk of overdose
OverdoseRisk  <- read.xlsx(WB, sheet="OverdoseRisk")
vparameters$od.preb.sub   <- with(OverdoseRisk, pe[par == "od.preb.sub"])
vparameters$od.il.lr.sub  <- with(OverdoseRisk, pe[par == "od.il.lr.sub"])
vparameters$od.NODU.sub   <- with(OverdoseRisk, pe[par == "od.NODU.sub"])
vparameters$multi.hr      <- with(OverdoseRisk, pe[par == "multi.hr"])
vparameters$multi.fx      <- with(OverdoseRisk, pe[par == "multi.fx"])
vparameters$multi.relap   <- with(OverdoseRisk, pe[par == "multi.relap"])
vparameters$multi.sub     <- with(OverdoseRisk, pe[par == "multi.sub"])
vparameters$multi.NODU.fx <- with(OverdoseRisk, pe[par == "multi.NODU.fx"])
#transition probability
TransProb     <- read.xlsx(WB, sheet="TransProb")
vparameters$p.preb2il.lr  <- with(TransProb, pe[par == "p.preb2il.lr"])
vparameters$p.preb2inact  <- with(TransProb, pe[par == "p.preb2inact"])
vparameters$p.il.lr2il.hr <- with(TransProb, pe[par == "p.il.lr2il.hr"])
vparameters$p.il.lr2inact <- with(TransProb, pe[par == "p.il.lr2inact"])
vparameters$p.il.hr2il.lr <- with(TransProb, pe[par == "p.il.hr2il.lr"])
vparameters$p.il.hr2inact <- with(TransProb, pe[par == "p.il.hr2inact"])
vparameters$p.inact2relap <- with(TransProb, pe[par == "p.inact2relap"])


## Parameters for decision tree ##
if (exists("sw.EMS.ODloc")){
  if (sw.EMS.ODloc == "sp"){
    OD_loc_priv         <- read.xlsx(WB, sheet="ODSettingEMS(sp)")$private
    OD_loc_pub          <- read.xlsx(WB, sheet="ODSettingEMS(sp)")$public
    OD_loc  <- rbind(OD_loc_priv, OD_loc_pub)
  } else {
    OD_loc_priv         <- read.xlsx(WB, sheet="ODSettingEMS")$private
    OD_loc_pub          <- read.xlsx(WB, sheet="ODSettingEMS")$public
    OD_loc  <- rbind(rep(OD_loc_priv, length(v.rgn)), rep(OD_loc_pub, length(v.rgn)))
  }
} else {
  OD_loc_priv         <- read.xlsx(WB, sheet="ODSettingEMS")$private
  OD_loc_pub          <- read.xlsx(WB, sheet="ODSettingEMS")$public
  OD_loc  <- rbind(rep(OD_loc_priv, length(v.rgn)), rep(OD_loc_pub, length(v.rgn)))
}
rownames(OD_loc)      <- c("priv", "pub")
colnames(OD_loc)      <- v.rgn
vparameters$OD_loc    <- OD_loc

DecisionTree     <- read.xlsx(WB, sheet="DecisionTree")
vparameters$OD_wit_priv   <- with(DecisionTree, pe[par == "OD_wit" & group == "priv"])
vparameters$OD_wit_pub    <- with(DecisionTree, pe[par == "OD_wit" & group == "pub"])
vparameters$OD_911_priv   <- with(DecisionTree, pe[par == "OD_911_priv"])
vparameters$OD_911_pub_mul<- with(DecisionTree, pe[par == "OD_911_pub_mul"])
vparameters$OD_911_pub    <- vparameters$OD_911_priv * vparameters$OD_911_pub_mul
vparameters$OD_hosp       <- with(DecisionTree, pe[par == "OD_hosp"])
vparameters$OD_cess       <- with(DecisionTree, pe[par == "OD_cess"])

Mortality     <- read.xlsx(WB, sheet="Mortality")
vparameters$mor_bl        <- with(Mortality, pe[par == "mor_bl"])
vparameters$mor_Nx        <- with(Mortality, pe[par == "mor_Nx"])
vparameters$rr_mor_EMS    <- with(Mortality, pe[par == "rr_mor_EMS"])

## Parameters for naloxone kits ##
NxKit         <- read.xlsx(WB, sheet="NxKit")
vparameters$r.LossExp     <- 1/with(NxKit, pe[par == "LossExp"])
vparameters$Low2Priv      <- with(NxKit, pe[par == "Low2Priv"])
vparameters$nlx.adj       <- with(NxKit, pe[par == "nlx.adj"])
vparameters$cap           <- with(NxKit, pe[par == "cap"])
NxDataOEND    <- read.xlsx(WB, sheet="NxDataOEND")
NxOEND.array  <- array(0, dim = c(length(unique(NxDataOEND$year)), length(unique(NxDataOEND$risk)), length(v.rgn)))
dimnames(NxOEND.array)[[1]] <- unique(NxDataOEND$year)
dimnames(NxOEND.array)[[2]] <- unique(NxDataOEND$risk)
dimnames(NxOEND.array)[[3]] <- v.rgn
for (i in 1:length(unique(NxDataOEND$year))){
  NxOEND.array[i, , ] <- data.matrix(subset(NxDataOEND, year == unique(NxDataOEND$year)[i])[-c(1,2)])
}
vparameters$NxOEND.array <- NxOEND.array
vparameters$NxDataPharm  <- read.xlsx(WB, sheet="NxDataPharm")
NxMvt                    <- data.matrix(read.xlsx(WB, sheet="NxMvt")[,-1])
row.names(NxMvt)         <- v.rgn
vparameters$NxMvt        <- NxMvt

## Parameters for cost ##
Cost          <- read.xlsx(WB, sheet="Cost")
vparameters$c.preb        <- with(Cost, pe[par == "c.preb"])
vparameters$c.il.lr       <- with(Cost, pe[par == "c.il.lr"])
vparameters$c.il.hr       <- with(Cost, pe[par == "c.il.hr"])
vparameters$c.inact       <- with(Cost, pe[par == "c.inact"])
vparameters$c.NODU        <- with(Cost, pe[par == "c.NODU"])
vparameters$c.nlx.kit     <- with(Cost, pe[par == "c.nlx.kit"])
vparameters$c.nlx.dtb     <- with(Cost, pe[par == "c.nlx.dtb"])
c.relap.v          <- numeric(0)                      # cost of remaining one cycle: relapsed, as the average of inactive and prior state
c.relap.v["preb"]  <- (vparameters$c.preb  + vparameters$c.inact)/2
c.relap.v["il.lr"] <- (vparameters$c.il.lr + vparameters$c.inact)/2
c.relap.v["il.hr"] <- (vparameters$c.il.hr + vparameters$c.inact)/2
vparameters$c.relap.v     <- c.relap.v

vparameters$c.EMS         <- with(Cost, pe[par == "c.EMS"])
vparameters$c.hospcare    <- with(Cost, pe[par == "c.hospcare"])

# Overdose probability matrix (per month)
od.matrix             <- matrix(0, nrow = 4, ncol = 2)
rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
colnames(od.matrix)   <- c("first", "subs")
od.matrix["preb", "subs"]   <- vparameters$od.preb.sub
od.matrix["il.lr", "subs"]  <- vparameters$od.il.lr.sub
od.matrix["il.hr", "subs"]  <- vparameters$od.il.lr.sub * vparameters$multi.hr
od.matrix["NODU", "subs"]   <- vparameters$od.NODU.sub
od.matrix[ , "first"]       <- od.matrix[ , "subs"] / vparameters$multi.sub
vparameters$od.matrix       <- od.matrix

# Baseline mortality excluding overdose (per month)
mor.matrix                  <- matrix(0, nrow = 2, ncol = length(mor.gp))
rownames(mor.matrix)        <- c("bg", "drug")
colnames(mor.matrix)        <- mor.gp
mor.matrix["bg", ]          <- vparameters$mor.bg
mor.matrix["drug", ]        <- vparameters$mor.drug
vparameters$mor.matrix      <- mor.matrix